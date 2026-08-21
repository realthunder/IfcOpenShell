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
			logger_.message(ifcopenshell::logger::LOG_WARNING, "GEO", 323,
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
		logger_.message(ifcopenshell::logger::LOG_WARNING, "GEO", 327,
			"Material layer slices do not add up to the body they came from, leaving it unsliced");
		return false;
	}

	items.swap(sliced);
	return true;
}

bool open_cascade_kernel::apply_folded_layerset(std::vector<conversion_result>& items, const layerset_information& info, const std::map<express::base, layerset_information>&) {
	// Folding the layers around the corner where two walls meet takes the
	// neighbours' boundaries as well and is a feature of its own. Until it is
	// written, a wall that has neighbours is sliced as though it had none,
	// which differs only within a wall thickness of the junction -- and is
	// what every such wall got before, minus the slicing.
	return apply_layerset(items, info);
}
