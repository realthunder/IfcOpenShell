/********************************************************************************
 *                                                                              *
 * This file is part of IfcOpenShell.                                           *
 *                                                                              *
 * IfcOpenShell is free software: you can redistribute it and/or modify         *
 * it under the terms of the Lesser GNU General Public License as published by  *
 * the Free Software Foundation, either version 3.0 of the License, or          *
 * (at your option) any later version.                                          *
 *                                                                              *
 * IfcOpenShell is distributed in the hope that it will be useful,              *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of               *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                 *
 * Lesser GNU General Public License for more details.                          *
 *                                                                              *
 * You should have received a copy of the Lesser GNU General Public License     *
 * along with this program. If not, see <http://www.gnu.org/licenses/>.         *
 *                                                                              *
 ********************************************************************************/

// This file was generated with the assistance of an AI coding tool.

// Slicing a product into its material layers, on the OpenCASCADE kernel.
//
// Both ends of this have been in the tree all along. The mapping builds one
// boundary per layer -- a plane, for anything whose body is a single
// extrusion, and for a wall the Axis curve pushed sideways by the accumulated
// layer offset -- and util::apply_layerset splits the body by those
// boundaries and gives each slice its layer's style. What was missing was the
// middle: nothing turned a boundary into a surface. abstract_kernel threw
// not_implemented_error, the hybrid kernel swallowed it, and every product
// with a material layer set came out as one undivided body. That is the state
// IfcOpenShell has been in since the taxonomy rewrite -- see issue #6607,
// where the same setting also lost the products outright, which is the null
// Axis mapping fixed separately in mapping.cpp.
//
// A wall boundary is a planar offset, and a planar offset here goes through
// libarea rather than BRepOffsetAPI. Not for the joins: for stability.
// BRepOffsetAPI_MakeOffset offsets locally and then repairs the topology,
// which is exactly what it is worst at; libarea's route rebuilds the result
// from scratch. The side an offset falls on is then libarea's answer rather
// than a normal derived by hand here, which is the same answer for a straight
// axis and a curved one.

#include "area_2d.h"
#include "layerset.h"
#include "opencascade_conversion_result.h"
#include "opencascade_kernel.h"

#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepGProp.hxx>
#include <BRep_Tool.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Edge.hxx>
#include <Geom_Curve.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <string>
#include <limits>
#include <GProp_GProps.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Plane.hxx>
#include <Precision.hxx>
#include <TopoDS_Wire.hxx>
#include <gp.hxx>
#include <gp_Ax3.hxx>
#include <gp_Dir.hxx>
#include <gp_Vec.hxx>

#include <cmath>
#include <iterator>
#include <map>
#include <list>
#include <mutex>

using namespace ifcopenshell::geom;
using namespace ifcopenshell::geom::kernels;

namespace {

	// A wall's Axis is a Curve2D in the product's own coordinates, in the
	// plane the wall extrudes out of, which is z = 0 there. The body being
	// split is in those same coordinates -- the product placement is applied
	// after this, by the iterator -- so nothing is placed here. The offset
	// curve the mapping hands over carries the product placement in its
	// matrix, and that is why it is not read.
	const gp_Ax3& axis_frame() {
		static const gp_Ax3 f(gp::Origin(), gp::DZ(), gp::DX());
		return f;
	}

	TopoDS_Wire wire_of(open_cascade_kernel* kernel, const taxonomy::ptr& curve) {
		open_cascade_kernel::curve_creation_visitor_result_type r;
		try {
			r = kernel->convert_curve(curve);
		} catch (const std::exception&) {
			return TopoDS_Wire();
		}

		if (r.index() == 1) {
			BRepBuilderAPI_MakeWire mw(BRepBuilderAPI_MakeEdge(std::get<Handle(Geom_Curve)>(r)).Edge());
			return mw.IsDone() ? mw.Wire() : TopoDS_Wire();
		} else if (r.index() == 2) {
			return std::get<TopoDS_Wire>(r);
		}

		return TopoDS_Wire();
	}

	// The surface a layer boundary lies on, unbounded: apply_layerset bounds it
	// to whatever it is about to cut, so an infinite surface is what it wants
	// and is also what makes a single span enough to describe the boundary.
	Handle(Geom_Surface) boundary_surface(open_cascade_kernel* kernel, const taxonomy::ptr& layer, double eps) {
		if (auto pln = taxonomy::dcast<taxonomy::plane>(layer)) {
			// A body that is one extrusion: its layers are planes already, and
			// no offset is involved.
			return kernel->convert_surface(pln);
		}

		taxonomy::ptr basis = layer;
		double offset = 0.;
		if (auto oc = taxonomy::dcast<taxonomy::offset_curve>(layer)) {
			basis = oc->basis;
			offset = oc->offset;
		}

		if (!basis) {
			return Handle(Geom_Surface)();
		}

		const TopoDS_Wire w = wire_of(kernel, basis);
		if (w.IsNull()) {
			return Handle(Geom_Surface)();
		}

#ifdef IFOPSH_WITH_LIBAREA
		const gp_Ax3& f = axis_frame();

		CCurve axis;
		{
			// libarea keeps its tolerances and its Clipper scratch state in
			// statics, so one caller at a time, shared with the 2D boolean.
			std::lock_guard<std::mutex> guard(util::area2d::lock());
			util::area2d::configure(eps);

			if (!util::area2d::wire_to_curve(w, f, false, axis)) {
				return Handle(Geom_Surface)();
			}
			// An offset curve's offset is measured along tangent x reference,
			// which is to the right of the curve; libarea offsets to the left.
			// Measured on walls whose layers are not palindromic, the stack
			// comes out mirrored about the reference line without this.
			if (offset != 0. && !axis.Offset(-offset)) {
				return Handle(Geom_Surface)();
			}
		}

		// An unbounded surface can carry one line or one circle, and one of
		// those is every wall axis seen so far: a two point polyline, or a
		// trimmed circle. An axis with joints in it would need a boundary made
		// of several surfaces, which apply_layerset has no way to take.
		const std::list<CVertex>& vs = axis.m_vertices;
		if (vs.size() < 2) {
			return Handle(Geom_Surface)();
		}

		if (vs.size() == 2 && vs.back().m_type == 0) {
			const gp_Pnt a = util::area2d::to_3d(vs.front().m_p, f);
			const gp_Pnt b = util::area2d::to_3d(vs.back().m_p, f);
			const gp_Vec d(a, b);
			if (d.Magnitude() <= Precision::Confusion()) {
				return Handle(Geom_Surface)();
			}
			// Across the wall rather than along it, and vertical: the wall
			// stands on this plane's line.
			return new Geom_Plane(a, gp_Dir(d) ^ f.Direction());
		}

		// A circle wider than half a turn arrives as two spans, because a
		// CVertex arc is an end point and a centre and cannot say which way
		// round it goes otherwise. Both spans are the same circle, so the
		// boundary is still one cylinder.
		Point centre(0., 0.);
		double r = 0.;
		Point prev = vs.front().m_p;
		for (auto it = std::next(vs.begin()); it != vs.end(); ++it) {
			if (it->m_type == 0) {
				return Handle(Geom_Surface)();
			}

			Point c = it->m_c;
			double ri = 0.;
			util::area2d::reconcile_arc(prev, it->m_p, c, ri);

			if (r == 0.) {
				centre = c;
				r = ri;
			} else if (c.dist(centre) > eps || std::fabs(ri - r) > eps) {
				// not one circle after all
				return Handle(Geom_Surface)();
			}

			prev = it->m_p;
		}

		if (r <= Precision::Confusion()) {
			return Handle(Geom_Surface)();
		}

		return new Geom_CylindricalSurface(gp_Ax3(util::area2d::to_3d(centre, f), f.Direction()), r);
#else
		(void) offset;
		(void) eps;
		return Handle(Geom_Surface)();
#endif
	}

	// Where a piece sits across the layer stack, as a distance from the first
	// boundary. Every boundary is an offset of the same basis -- parallel
	// planes for a slab, coaxial cylinders for a curved wall -- so this single
	// coordinate places a piece in the stack whatever the boundaries are.
	bool extent_across_layers(const TopoDS_Shape& shape, const Handle(Geom_Surface)& first, double& lo, double& hi, std::string* worst = nullptr) {
		lo = std::numeric_limits<double>::max();
		hi = -std::numeric_limits<double>::max();
		int n = 0;

		double worst_d = -1.;
		auto sample = [&](const gp_Pnt& p, const char* what, double detail) {
			GeomAPI_ProjectPointOnSurf proj(p, first);
			if (!proj.IsDone() || proj.NbPoints() == 0) {
				return false;
			}
			const double d = proj.LowerDistance();
			lo = std::min(lo, d);
			hi = std::max(hi, d);
			if (worst && d > worst_d) {
				worst_d = d;
				*worst = std::string(what) + " " + std::to_string(detail) + " at " + std::to_string(d) +
					" (" + std::to_string(p.X()) + ", " + std::to_string(p.Y()) + ", " + std::to_string(p.Z()) + ")";
			}
			return true;
		};

		for (TopExp_Explorer exp(shape, TopAbs_VERTEX); exp.More(); exp.Next(), ++n) {
			if (!sample(BRep_Tool::Pnt(TopoDS::Vertex(exp.Current())), "vertex", 0.)) {
				return false;
			}
		}

		// Along the edges as well, not only at their ends. A split by a curved
		// boundary can leave a piece carrying a face whose edge runs off as a
		// straight chord while its two vertices stay where they belong -- a fin
		// of no volume, so neither the vertices nor the volume of the slices
		// show it, and only a point taken between the ends does.
		for (TopExp_Explorer exp(shape, TopAbs_EDGE); exp.More(); exp.Next()) {
			double a, b;
			opencascade::handle<Geom_Curve> crv = BRep_Tool::Curve(TopoDS::Edge(exp.Current()), a, b);
			if (crv.IsNull()) {
				continue;
			}
			GeomAdaptor_Curve gac(crv);
			for (int k = 1; k < 4; ++k) {
				gp_Pnt p;
				crv->D0(a + (b - a) * k / 4., p);
				if (!sample(p, "edge of curve kind", (double)gac.GetType() * 10 + k)) {
					return false;
				}
			}
		}

		return n > 0;
	}

	// The curve a layer boundary is an offset of. Every boundary of a wall is
	// an offset of the one axis, which is what lets a single coordinate place
	// anything across the stack.
	taxonomy::ptr boundary_basis(const taxonomy::ptr& layer) {
		if (auto oc = taxonomy::dcast<taxonomy::offset_curve>(layer)) {
			return oc->basis;
		}
		return layer;
	}

	// The product placement the mapping left on every boundary it built. It is
	// deliberately not read when a body is sliced on its own -- body and
	// boundaries are in the same coordinates there -- but a neighbour's
	// boundaries arrive in the neighbour's coordinates and have to be carried
	// into ours.
	bool placement_of(const taxonomy::ptr& layer, Eigen::Matrix4d& m) {
		auto g = taxonomy::dcast<taxonomy::geom_item>(layer);
		if (!g) {
			return false;
		}
		m = g->matrix ? g->matrix->ccomponents() : Eigen::Matrix4d::Identity();
		return true;
	}

	// Both are placements, so the step between them is a rigid motion. A
	// matrix that is not one is refused rather than approximated: a scaled or
	// sheared fold surface would cut the body in the wrong place, and cutting
	// in the wrong place is worse than not folding.
	bool relative_placement(const layerset_information& ours, const layerset_information& theirs, gp_Trsf& trsf) {
		Eigen::Matrix4d a, b;
		if (ours.layers.empty() || theirs.layers.empty()) {
			return false;
		}
		if (!placement_of(ours.layers.front(), a) || !placement_of(theirs.layers.front(), b)) {
			return false;
		}

		const Eigen::Matrix4d rel = a.inverse() * b;
		const Eigen::Matrix3d r = rel.block<3, 3>(0, 0);
		if ((r.transpose() * r - Eigen::Matrix3d::Identity()).norm() > 1.e-6) {
			return false;
		}

		trsf.SetValues(
			rel(0, 0), rel(0, 1), rel(0, 2), rel(0, 3),
			rel(1, 0), rel(1, 1), rel(1, 2), rel(1, 3),
			rel(2, 0), rel(2, 1), rel(2, 2), rel(2, 3));
		return true;
	}

	// The ends of the wall axis, in the wall's own coordinates.
	bool axis_end_points(open_cascade_kernel* kernel, const layerset_information& info, gp_Pnt& a, gp_Pnt& b) {
		if (info.layers.empty()) {
			return false;
		}
		const taxonomy::ptr basis = boundary_basis(info.layers.front());
		if (!basis) {
			return false;
		}
		const TopoDS_Wire w = wire_of(kernel, basis);
		if (w.IsNull()) {
			return false;
		}
		TopoDS_Vertex v1, v2;
		TopExp::Vertices(w, v1, v2);
		if (v1.IsNull() || v2.IsNull()) {
			return false;
		}
		a = BRep_Tool::Pnt(v1);
		b = BRep_Tool::Pnt(v2);
		return a.Distance(b) > Precision::Confusion();
	}

	// The point of a surface nearest a given point.
	bool foot_on(const Handle(Geom_Surface)& s, const gp_Pnt& from, gp_Pnt& foot) {
		GeomAPI_ProjectPointOnSurf proj(from, s);
		if (!proj.IsDone() || proj.NbPoints() == 0) {
			return false;
		}
		foot = proj.NearestPoint();
		return true;
	}

	// Where a boundary sits when measured from a point along a direction. This
	// asks the surface for its nearest point rather than for a normal, so a
	// plane and a cylinder answer the same way.
	bool score_along(const Handle(Geom_Surface)& s, const gp_Pnt& from, const gp_Vec& dir, double& score) {
		GeomAPI_ProjectPointOnSurf proj(from, s);
		if (!proj.IsDone() || proj.NbPoints() == 0) {
			return false;
		}
		score = gp_Vec(from, proj.NearestPoint()).Dot(dir);
		return true;
	}

	// A surface with its normal made to point a chosen way.
	//
	// This matters because apply_folded_layerset decides which part of each
	// surface to keep from the surface's own normal: it keeps the boundary on
	// the far side of the fold and the fold on the far side of the boundary,
	// and those two half-space choices are what make the pair turn a corner
	// one way rather than another. The normal a boundary happens to arrive
	// with is whichever way its axis ran, so left alone the fold turns at
	// random -- three of the four combinations cut the wall in a way its
	// layers never asked for, and one of them drops the boundary along the
	// whole length of the wall and keeps only the corner.
	Handle(Geom_Surface) oriented_surface(const Handle(Geom_Surface)& s, const gp_Pnt& near, const gp_Vec& want) {
		GeomAPI_ProjectPointOnSurf proj(near, s);
		if (!proj.IsDone() || proj.NbPoints() == 0) {
			return s;
		}
		double u, v;
		proj.LowerDistanceParameters(u, v);
		gp_Pnt p;
		gp_Vec du, dv;
		s->D1(u, v, p, du, dv);
		const gp_Vec n = du ^ dv;
		if (n.Magnitude() <= Precision::Confusion() || n.Dot(want) >= 0.) {
			return s;
		}
		return s->UReversed();
	}

	// Give each interior boundary the neighbours' boundaries that continue it
	// around a corner. Returns the number of folds made; nothing to fold is
	// not a failure, it just means the wall is sliced as a wall on its own.
	//
	// The pairing is by thickness from the outside of the corner. Two walls
	// meeting at an end have their layers reconciled there: the cladding on
	// the outside of the corner runs on into the neighbour's cladding, and so
	// on inwards, so our boundary a given distance in from the outer face
	// continues as the neighbour's boundary the same distance in from theirs.
	// Which end of a stack is the outer one is not in the file -- it is
	// geometry -- and the neighbour's own direction answers it: the neighbour
	// extends away from the junction, so the outside of the corner is the end
	// of our stack that lies against that direction.
	int fold_boundaries(open_cascade_kernel* kernel, const layerset_information& info,
		const std::vector<Handle(Geom_Surface)>& surfaces,
		const std::map<express::base, layerset_information>& neighbours, double eps,
		std::vector<std::vector<Handle(Geom_Surface)>>& folded) {

		folded.assign(surfaces.size(), std::vector<Handle(Geom_Surface)>());
		for (size_t i = 0; i < surfaces.size(); ++i) {
			folded[i].push_back(surfaces[i]);
		}

		std::vector<double> cum(1, 0.);
		for (auto& t : info.thicknesses) {
			cum.push_back(cum.back() + t);
		}
		if (cum.size() != surfaces.size() || cum.back() <= Precision::Confusion()) {
			return 0;
		}

		gp_Pnt our_1, our_2;
		if (!axis_end_points(kernel, info, our_1, our_2)) {
			return 0;
		}

		int made = 0;
		bool chosen = false;
		gp_Vec chosen_inner;
		gp_Pnt chosen_at;

		for (auto& kv : neighbours) {
			const layerset_information& other = kv.second;
			if (other.layers.size() < 3 || other.layers.size() != other.thicknesses.size() + 1) {
				continue;
			}

			gp_Trsf to_ours;
			if (!relative_placement(info, other, to_ours)) {
				continue;
			}

			gp_Pnt their_1, their_2;
			if (!axis_end_points(kernel, other, their_1, their_2)) {
				continue;
			}
			their_1.Transform(to_ours);
			their_2.Transform(to_ours);

			// Which ends meet is measured rather than taken from the
			// connection type: a file that says ATSTART where the geometry
			// says otherwise would fold the wrong end of the wall.
			const gp_Pnt* ours[2] = { &our_1, &our_2 };
			const gp_Pnt* theirs[2] = { &their_1, &their_2 };
			int oi = 0, ti = 0;
			double best = std::numeric_limits<double>::max();
			for (int a = 0; a < 2; ++a) {
				for (int b = 0; b < 2; ++b) {
					const double d = ours[a]->Distance(*theirs[b]);
					if (d < best) {
						best = d;
						oi = a;
						ti = b;
					}
				}
			}

			// The junction has to be a junction. Walls whose axes stop short of
			// one another by more than the wall is thick are not meeting at a
			// corner, whatever the file says about them.
			if (best > cum.back() * 2.) {
				continue;
			}

			const gp_Pnt& junction = *ours[oi];
			gp_Vec our_away(*ours[oi], *ours[1 - oi]);
			gp_Vec their_away(*theirs[ti], *theirs[1 - ti]);
			if (our_away.Magnitude() <= Precision::Confusion() || their_away.Magnitude() <= Precision::Confusion()) {
				continue;
			}
			our_away.Normalize();
			their_away.Normalize();

			std::vector<double> our_score(surfaces.size());
			bool scored = true;
			for (size_t i = 0; i < surfaces.size() && scored; ++i) {
				scored = score_along(surfaces[i], junction, their_away, our_score[i]);
			}
			if (!scored) {
				continue;
			}

			// Two walls running on into one another have no corner: their
			// directions are the same line, our stack does not spread out
			// along theirs, and there is no outer end to pair from.
			if (std::fabs(our_score.back() - our_score.front()) < cum.back() * .5) {
				continue;
			}

			std::vector<Handle(Geom_Surface)> their_surfaces;
			for (auto& layer : other.layers) {
				Handle(Geom_Surface) s = boundary_surface(kernel, layer, eps);
				if (s.IsNull()) {
					their_surfaces.clear();
					break;
				}
				their_surfaces.push_back(Handle(Geom_Surface)::DownCast(s->Transformed(to_ours)));
			}
			if (their_surfaces.size() != other.layers.size()) {
				continue;
			}

			std::vector<double> their_cum(1, 0.);
			for (auto& t : other.thicknesses) {
				their_cum.push_back(their_cum.back() + t);
			}

			std::vector<double> their_score(their_surfaces.size());
			for (size_t j = 0; j < their_surfaces.size() && scored; ++j) {
				scored = score_along(their_surfaces[j], junction, our_away, their_score[j]);
			}
			if (!scored || std::fabs(their_score.back() - their_score.front()) < their_cum.back() * .5) {
				continue;
			}

			const bool our_front_outer = our_score.front() < our_score.back();
			const bool their_front_outer = their_score.front() < their_score.back();

			// Across our own stack, from the face on the outside of this
			// corner towards the far one. Every fold at this junction turns
			// the same way about it.
			gp_Pnt foot_first, foot_last;
			if (!foot_on(surfaces.front(), junction, foot_first) || !foot_on(surfaces.back(), junction, foot_last)) {
				continue;
			}
			gp_Vec inner(foot_first, foot_last);
			if (!our_front_outer) {
				inner.Reverse();
			}
			if (inner.Magnitude() <= Precision::Confusion()) {
				continue;
			}
			inner.Normalize();

			// A wall can be on the outside of the corner at one end and on the
			// inside at the other -- a Z of three walls is the everyday case.
			// Its boundaries would then have to turn both ways at once, which
			// one surface cannot do, so the second junction is left unfolded
			// rather than folded backwards.
			if (chosen) {
				if (inner.Dot(chosen_inner) <= 0.) {
					continue;
				}
			} else {
				chosen = true;
				chosen_inner = inner;
				chosen_at = junction;
			}

			const double tol = std::max(eps * 10., cum.back() * 1.e-3);

			for (size_t i = 1; i + 1 < surfaces.size(); ++i) {
				const double target = our_front_outer ? cum[i] : cum.back() - cum[i];

				for (size_t j = 1; j + 1 < their_surfaces.size(); ++j) {
					const double d = their_front_outer ? their_cum[j] : their_cum.back() - their_cum[j];
					if (std::fabs(d - target) <= tol) {
						// The fold keeps the boundary on its far side from the
						// junction, so its normal has to point down the wall
						// away from the corner.
						folded[i].push_back(oriented_surface(their_surfaces[j], junction, our_away));
						++made;
						break;
					}
				}
			}
		}

		// The boundary keeps the fold on its inner side, so where anything was
		// folded the boundary's own normal has to point across the stack.
		if (made > 0) {
			for (size_t i = 1; i + 1 < surfaces.size(); ++i) {
				if (folded[i].size() > 1) {
					folded[i][0] = oriented_surface(folded[i][0], chosen_at, chosen_inner);
				}
			}
		}

		return made;
	}

	double total_volume(const std::vector<conversion_result>& items) {
		double v = 0.;
		for (auto& item : items) {
			GProp_GProps props;
			BRepGProp::VolumeProperties(std::static_pointer_cast<open_cascade_shape>(item.shape())->shape(), props);
			v += props.Mass();
		}
		return v;
	}

}

bool open_cascade_kernel::apply_layerset(std::vector<conversion_result>& items, const layerset_information& info) {
	// Two boundaries and one layer between them is the whole body: there is
	// nothing to divide, and saying so is not a failure.
	if (info.layers.size() < 3 || info.layers.size() != info.styles.size() + 1) {
		return false;
	}

	std::vector<Handle(Geom_Surface)> surfaces;
	for (auto& layer : info.layers) {
		auto s = boundary_surface(this, layer, precision_);
		if (s.IsNull()) {
			logger_.message(ifcopenshell::logger::LOG_WARNING, "GEO", 330,
				"Unable to build a surface for a material layer boundary of kind "
					+ taxonomy::kind_to_string(layer->kind()) + ", the body is left unsliced");
			return false;
		}
		surfaces.push_back(s);
	}

	std::vector<taxonomy::style::ptr> styles;
	for (auto& style : info.styles) {
		styles.push_back(taxonomy::make<taxonomy::style>(style));
	}

	// The body itself has to be the thickness of its own layer set before there
	// is any point dividing it. Bodies that are not are common -- King has a
	// wall 48 units thick carrying a 0.4167 layer set -- and slicing one hands
	// back an outermost piece that silently absorbs everything the layer set
	// never described.
	std::vector<double> declared(1, 0.);
	for (auto& t : info.thicknesses) {
		declared.push_back(declared.back() + t);
	}
	const double band_tol = std::max(precision_ * 10., declared.back() * 1.e-3);

	if (declared.size() == surfaces.size()) {
		for (auto& item : items) {
			double lo, hi;
			std::string worst;
			const TopoDS_Shape& shape = std::static_pointer_cast<open_cascade_shape>(item.shape())->shape();
			if (!extent_across_layers(shape, surfaces.front(), lo, hi, &worst) ||
				std::fabs(lo) > band_tol || std::fabs(hi - declared.back()) > band_tol) {
				logger_.message(ifcopenshell::logger::LOG_WARNING, "GEO", 337,
					"The body measures " + std::to_string(lo) + ".." + std::to_string(hi) +
					" across its own layer boundaries, not 0.." + std::to_string(declared.back()) +
					", so its material layer set does not describe it and it is left unsliced"
					" [furthest point: " + worst + "]");
				return false;
			}
		}
	}

	std::vector<conversion_result> sliced;
	if (!util::apply_layerset(items, surfaces, styles, sliced, precision_)) {
		return false;
	}

	// The slices are the same material as the body, rearranged, so they have
	// to add up to it. A split by a curved boundary can hand back pieces that
	// overlap -- measured at up to 62 per cent too much on a curved wall --
	// and a body that has silently gained material is worse than one that was
	// never divided, because nothing downstream can tell.
	const double before = total_volume(items);
	const double after = total_volume(sliced);
	if (before <= 0. || std::fabs(after - before) > 1.e-3 * before) {
		logger_.message(ifcopenshell::logger::LOG_WARNING, "GEO", 331,
			"Material layer slices do not add up to the body they came from, leaving it unsliced");
		return false;
	}

	// Every piece has to span exactly one declared layer, measured across the
	// boundaries. Without this a body that is not the thickness of its own
	// layer set still divides, and the outermost piece silently absorbs all
	// the material the layer set never described -- King has curved walls
	// whose bodies are twenty-eight times their layer set, which came out as
	// two right layers and one 11-unit one against a declared 0.0625. The
	// volume guard above does not see this: the pieces do add up to the body.
	const std::vector<double>& cumulative = declared;

	if (cumulative.size() == surfaces.size()) {
		for (auto& piece : sliced) {
			double lo, hi;
			const TopoDS_Shape& shape = std::static_pointer_cast<open_cascade_shape>(piece.shape())->shape();
			bool spans_a_layer = false;
			std::string worst;
			if (extent_across_layers(shape, surfaces.front(), lo, hi, &worst)) {
				for (size_t i = 0; i + 1 < cumulative.size(); ++i) {
					if (std::fabs(lo - cumulative[i]) <= band_tol && std::fabs(hi - cumulative[i + 1]) <= band_tol) {
						spans_a_layer = true;
						break;
					}
				}
			}
			if (!spans_a_layer) {
				logger_.message(ifcopenshell::logger::LOG_WARNING, "GEO", 335,
					"A material layer slice measures " + std::to_string(lo) + ".." + std::to_string(hi) +
					" across the boundaries, which is not the thickness of any declared layer, leaving the body unsliced"
					" [furthest point: " + worst + "]");
				return false;
			}
		}
	}

	items.swap(sliced);
	return true;
}

bool open_cascade_kernel::apply_folded_layerset(std::vector<conversion_result>& items, const layerset_information& info, const std::map<express::base, layerset_information>& neighbours) {
	if (info.layers.size() < 3 || info.layers.size() != info.styles.size() + 1) {
		return false;
	}

	std::vector<Handle(Geom_Surface)> surfaces;
	for (auto& layer : info.layers) {
		auto s = boundary_surface(this, layer, precision_);
		if (s.IsNull()) {
			// Reported by the unfolded path, which cannot get past this either.
			return apply_layerset(items, info);
		}
		surfaces.push_back(s);
	}

	std::vector<std::vector<Handle(Geom_Surface)>> folded;
	const int made = fold_boundaries(this, info, surfaces, neighbours, precision_, folded);
	if (made <= 0) {
		// A wall whose neighbours turn out to fold nothing -- they meet it end
		// on, or their layers line up with none of its own -- is a wall on its
		// own as far as the geometry goes.
		return apply_layerset(items, info);
	}

	// The same guard the unfolded path applies: a body that is not the
	// thickness of its own layer set is not described by it, and dividing it
	// hands back an outer piece that silently absorbs the rest.
	std::vector<double> declared(1, 0.);
	for (auto& t : info.thicknesses) {
		declared.push_back(declared.back() + t);
	}
	const double band_tol = std::max(precision_ * 10., declared.back() * 1.e-3);

	if (declared.size() == surfaces.size()) {
		for (auto& item : items) {
			double lo, hi;
			std::string worst;
			const TopoDS_Shape& shape = std::static_pointer_cast<open_cascade_shape>(item.shape())->shape();
			if (!extent_across_layers(shape, surfaces.front(), lo, hi, &worst) ||
				std::fabs(lo) > band_tol || std::fabs(hi - declared.back()) > band_tol) {
				logger_.message(ifcopenshell::logger::LOG_WARNING, "GEO", 337,
					"The body measures " + std::to_string(lo) + ".." + std::to_string(hi) +
					" across its own layer boundaries, not 0.." + std::to_string(declared.back()) +
					", so its material layer set does not describe it and it is left unsliced"
					" [furthest point: " + worst + "]");
				return false;
			}
		}
	}

	std::vector<taxonomy::style::ptr> styles;
	for (auto& style : info.styles) {
		styles.push_back(taxonomy::make<taxonomy::style>(style));
	}

	std::vector<conversion_result> sliced;
	if (!util::apply_folded_layerset(items, folded, styles, sliced, precision_)) {
		logger_.message(ifcopenshell::logger::LOG_WARNING, "GEO", 338,
			"Unable to divide the body by " + std::to_string(made) +
			" folded material layer boundaries, slicing it as a wall on its own instead");
		return apply_layerset(items, info);
	}

	// The pieces are the body rearranged, so they still have to add up to it.
	const double before = total_volume(items);
	const double after = total_volume(sliced);
	if (before <= 0. || std::fabs(after - before) > 1.e-3 * before) {
		logger_.message(ifcopenshell::logger::LOG_WARNING, "GEO", 339,
			"Folded material layer slices do not add up to the body they came from, slicing it as a wall on its own instead");
		return apply_layerset(items, info);
	}

	// A folded piece is deliberately not one layer band -- that is what the
	// fold is -- so it cannot be checked against a declared thickness the way
	// an unfolded slice is. What still holds is that no piece may reach
	// outside the layer set: the body was measured as 0..total above, and a
	// piece that leaves that range is material the layer set never described.
	for (auto& piece : sliced) {
		double lo, hi;
		std::string worst;
		const TopoDS_Shape& shape = std::static_pointer_cast<open_cascade_shape>(piece.shape())->shape();
		if (!extent_across_layers(shape, surfaces.front(), lo, hi, &worst) ||
			lo < -band_tol || hi > declared.back() + band_tol) {
			logger_.message(ifcopenshell::logger::LOG_WARNING, "GEO", 340,
				"A folded material layer slice measures " + std::to_string(lo) + ".." + std::to_string(hi) +
				" across the boundaries, which reaches outside the declared 0.." + std::to_string(declared.back()) +
				", slicing the body as a wall on its own instead [furthest point: " + worst + "]");
			return apply_layerset(items, info);
		}
	}

	logger_.message(ifcopenshell::logger::LOG_NOTICE, "GEO", 341,
		"The material layers were folded around " + std::to_string(made) + " boundary continuation(s) where this wall meets its neighbours");

	items.swap(sliced);
	return true;
}
