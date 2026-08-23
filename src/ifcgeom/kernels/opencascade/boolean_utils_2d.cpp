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

// A planar subtraction done on an integer lattice instead of by the boolean
// kernel. The kernel's face boolean is the same BOPAlgo pipeline as the 3D
// one, so a profile it dislikes costs the same tens of seconds and fails the
// same way: a second operand that starts at the first operand's own boundary
// is a notch in the outer contour rather than a hole in the face, and deriving
// that over a couple of hundred operands at once is what produces
// BRepCheck_InvalidImbricationOfWires. Clipper rebuilds contour nesting from
// scratch under exact integer predicates, so the notch case is not special to
// it.
//
// The lattice is reached through libarea rather than through Clipper directly.
// Clipper cannot return a curve -- everything that goes onto the lattice comes
// back as the chords it was walked along -- and a wall profile that comes back
// faceted is not merely less accurate: on the model this was written for it
// went from 96 faces to 2426, and every operation downstream of it then paid
// for that, which cost more than the boolean being replaced. libarea carries
// an arc as a centre and an end point, only walks it at the last moment, and
// runs FitArcs over the result to recognise the arcs again. That is the whole
// reason for the dependency.

#include "boolean_utils.h"

#include "../../../ifcparse/logger.h"

#ifndef IFOPSH_WITH_LIBAREA

// Built without libarea, so every caller falls through to the boolean kernel
// exactly as it did before this path existed.
bool ifcopenshell::geom::util::boolean_subtraction_2d_using_area(const TopoDS_Shape&, const NCollection_List<TopoDS_Shape>&, TopoDS_Shape&, double, bool, ifcopenshell::logger&) {
	return false;
}

bool ifcopenshell::geom::util::boolean_2d_area_supports(const TopoDS_Shape&) {
	return false;
}

#else

#include "area_2d.h"

#include <libarea/Area.h>
#include <clipper2/clipper.h>

#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBndLib.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepTools.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box.hxx>
#include <Geom_Plane.hxx>
#include <Precision.hxx>
#include <Standard_Failure.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>
#include <gp_Circ.hxx>
#include <gp_Pln.hxx>
#include <gp_Vec.hxx>

#include <cmath>
#include <list>
#include <mutex>
#include <string>
#include <vector>

using namespace ifcopenshell::geom::util::area2d;

namespace {
	// A right-handed frame on the face's plane. gp_Ax3 is allowed to be
	// left-handed, and the direction an arc turns is read off the sign of its
	// axis against this frame's normal, so taking the plane's own X and Y
	// without checking would silently reverse arcs on some faces.
	gp_Ax3 frame_of(const Handle(Geom_Plane)& plane) {
		const gp_Ax3& p = plane->Position();
		return gp_Ax3(p.Location(), p.Direction(), p.XDirection());
	}

	// Under the non-zero rule it is the winding that says which side of a
	// contour is solid, and nothing about the order OCCT hands its wires over
	// guarantees that -- a face whose own orientation is reversed hands them
	// over backwards. So say it explicitly: outer contours one way, holes the
	// other. Only the relative sense matters, since the fill rule tests the
	// winding for zero rather than for a sign.
	void wind(CCurve& curve, bool hole) {
		if (curve.IsClockwise() != hole) {
			curve.Reverse();
		}
	}

	bool face_to_area(const TopoDS_Face& face, const gp_Ax3& f, CArea& area) {
		const TopoDS_Wire outer = BRepTools::OuterWire(face);
		if (outer.IsNull()) {
			return false;
		}

		bool any = false;
		for (TopoDS_Iterator it(face); it.More(); it.Next()) {
			if (it.Value().ShapeType() != TopAbs_WIRE) {
				continue;
			}
			CCurve curve;
			if (!wire_to_curve(TopoDS::Wire(it.Value()), f, true, curve)) {
				// Dropping a boundary here would quietly hand back a shape
				// with a hole missing, which is worse than not answering.
				return false;
			}
			wind(curve, !it.Value().IsSame(outer));
			area.append(curve);
			any = true;
		}
		return any;
	}

	// Null on failure. A wire that cannot be rebuilt is not patched up with
	// straight edges: that would hand back a shape that is quietly the wrong
	// shape, and the kernel is still there to be asked instead.
	TopoDS_Wire curve_to_wire(const CCurve& c, const gp_Ax3& f) {
		BRepBuilderAPI_MakeWire mw;

		bool first = true;
		Point prev(0, 0);
		gp_Pnt prev3d;

		for (const CVertex& v : c.m_vertices) {
			if (first) {
				first = false;
				prev = v.m_p;
				prev3d = to_3d(prev, f);
				continue;
			}

			const gp_Pnt next3d = to_3d(v.m_p, f);
			if (next3d.SquareDistance(prev3d) <= Precision::SquareConfusion()) {
				// An arc that ends where it began is the whole circle -- what
				// libarea says when it has put one back together -- and a line
				// that ends where it began is nothing.
				if (v.m_type == 0) {
					continue;
				}

				const gp_Pnt centre3d = to_3d(v.m_c, f);
				const double r = centre3d.Distance(prev3d);
				if (r <= Precision::Confusion()) {
					continue;
				}

				gp_Dir axis = f.Direction();
				if (v.m_type < 0) {
					axis.Reverse();
				}

				try {
					mw.Add(BRepBuilderAPI_MakeEdge(gp_Circ(gp_Ax2(centre3d, axis), r)).Edge());
				} catch (const Standard_Failure&) {
					return TopoDS_Wire();
				}

				if (!mw.IsDone()) {
					return TopoDS_Wire();
				}

				continue;
			}

			try {
				if (v.m_type == 0) {
					mw.Add(BRepBuilderAPI_MakeEdge(prev3d, next3d).Edge());
				} else {
					Point centre = v.m_c;
					double r = 0.;
					reconcile_arc(prev, v.m_p, centre, r);

					gp_Dir axis = f.Direction();
					if (v.m_type < 0) {
						axis.Reverse();
					}

					mw.Add(BRepBuilderAPI_MakeEdge(gp_Circ(gp_Ax2(to_3d(centre, f), axis), r), prev3d, next3d).Edge());
				}
			} catch (const Standard_Failure&) {
				return TopoDS_Wire();
			}

			if (!mw.IsDone()) {
				return TopoDS_Wire();
			}

			prev = v.m_p;
			prev3d = next3d;
		}

		return mw.IsDone() ? mw.Wire() : TopoDS_Wire();
	}

	// Summed signed and made positive at the end, not per contour: a hole is
	// wound against its outer and so comes out negative, which is what makes
	// this the net area. Taking each contour's magnitude instead would count
	// holes as material, and since a cut is in the business of making holes
	// the total would grow with every one of them -- which reads as a cut that
	// added area, and rejects a correct result.
	double lattice_area(const Clipper2Lib::Paths64& paths) {
		double a = 0.;
		for (const auto& p : paths) {
			a += Clipper2Lib::Area(p);
		}
		return std::fabs(a);
	}

	TopoDS_Wire path_to_wire(const Clipper2Lib::Path64& contour, const gp_Ax3& f, bool hole) {
		// Back through libarea rather than straight off the lattice: this is
		// where SimplifyPath and FitArcs run, and a hole that came in as an
		// arc leaves as an arc rather than as the several hundred chords it
		// was walked along.
		CCurve curve;
		CurveFromClipperPath(curve, const_cast<Clipper2Lib::Path64&>(contour));

		// Clipper winds its output the way its own rules require, and the
		// conversion reverses on the way in and again on the way out, so what
		// arrives here says nothing about which side is solid. Wind it to say
		// so: MakeFace takes a face's sense from its outer wire, and an outer
		// contour that arrives clockwise builds a face pointing the wrong way,
		// which a prism raised on it turns into volume where there should be
		// none. That is not a small error -- it took the wall to more volume
		// than it had before anything was cut out of it.
		wind(curve, hole);
		return curve_to_wire(curve, f);
	}

	// A Clipper solution is a forest: solid contours, their holes, and any
	// islands standing inside those holes, which is a nesting that the flat
	// list of contours cannot express and that would otherwise have to be
	// worked out again from the geometry.
	bool faces_from_node(const Clipper2Lib::PolyPath64& node, const gp_Ax3& f, double eps, NCollection_List<TopoDS_Shape>& faces, const std::vector<TopoDS_Wire>& carry_through) {
		const TopoDS_Wire outer = path_to_wire(node.Polygon(), f, false);
		if (outer.IsNull()) {
			return false;
		}

		// On the plane of A rather than one deduced from the contour:
		// MakeFace(wire) takes its normal from the way the wire winds, so a
		// contour wound the other way would give an inside-out face, and a
		// prism raised on that has negative volume.
		BRepBuilderAPI_MakeFace mf(gp_Pln(f), outer);

		for (const auto& hole : node) {
			// Not reversed. Clipper already hands back a hole wound against
			// its outer contour, which is how a face wants an inner boundary.
			// The other 2D path in this file does reverse, because there the
			// wires come off an OCCT face, where an inner wire is already
			// oriented as a hole of that face -- the same wire from the two
			// sources arrives with opposite conventions.
			const TopoDS_Wire hw = path_to_wire(hole->Polygon(), f, true);
			if (hw.IsNull()) {
				return false;
			}
			mf.Add(hw);
		}

		if (!mf.IsDone()) {
			return false;
		}

		// A carried boundary belongs to the piece of the result that still
		// surrounds it. Nothing came near it, so exactly one piece does.
		if (!carry_through.empty()) {
			const TopoDS_Face probe = mf.Face();
			for (const TopoDS_Wire& w : carry_through) {
				BRepTools_WireExplorer exp(w);
				if (!exp.More()) {
					continue;
				}
				BRepClass_FaceClassifier fc(probe, BRep_Tool::Pnt(exp.CurrentVertex()), eps);
				if (fc.State() == TopAbs_IN) {
					mf.Add(w);
				}
			}

			if (!mf.IsDone()) {
				return false;
			}
		}

		faces.Append(mf.Face());

		// An island inside a hole is solid again, and a face of its own.
		for (const auto& hole : node) {
			for (const auto& island : *hole) {
				if (!faces_from_node(*island, f, eps, faces, carry_through)) {
					return false;
				}
			}
		}

		return true;
	}

	// Whether any operand comes near enough to this boundary to change it.
	// Bounding boxes only prune: a reveal is a long thin sliver whose box
	// covers most of a wall, so on boxes alone every window in that wall reads
	// as reachable. What settles it is the distance from the boundary itself,
	// and an operand lying wholly inside a hole or wholly outside it leaves
	// that hole's boundary alone either way.
	bool operand_reaches(const TopoDS_Wire& w, const NCollection_List<TopoDS_Shape>& b, const std::vector<Bnd_Box>& boxes, double eps) {
		Bnd_Box wb;
		BRepBndLib::Add(w, wb);
		wb.Enlarge(eps);

		size_t i = 0;
		for (NCollection_List<TopoDS_Shape>::Iterator it(b); it.More(); it.Next(), ++i) {
			if (i < boxes.size() && boxes[i].IsOut(wb)) {
				continue;
			}

			BRepExtrema_DistShapeShape d(w, it.Value(), Extrema_ExtFlag_MIN);
			if (d.IsDone() && d.Value() > eps) {
				continue;
			}
			return true;
		}

		return false;
	}

}

// Whether this operand can go onto the lattice without changing shape.
bool ifcopenshell::geom::util::boolean_2d_area_supports(const TopoDS_Shape& s) {
	if (s.ShapeType() != TopAbs_FACE) {
		return false;
	}
	if (Handle(Geom_Plane)::DownCast(BRep_Tool::Surface(TopoDS::Face(s))).IsNull()) {
		return false;
	}
	for (TopExp_Explorer exp(s, TopAbs_EDGE); exp.More(); exp.Next()) {
		if (!edge_is_representable(TopoDS::Edge(exp.Current()))) {
			return false;
		}
	}
	return true;
}

bool ifcopenshell::geom::util::boolean_subtraction_2d_using_area(const TopoDS_Shape& a_input, const NCollection_List<TopoDS_Shape>& b_input, TopoDS_Shape& result, double eps, bool fit_circles, ifcopenshell::logger& logger) {
	if (a_input.ShapeType() != TopAbs_FACE) {
		return false;
	}

	const TopoDS_Face& a = TopoDS::Face(a_input);
	auto plane = Handle(Geom_Plane)::DownCast(BRep_Tool::Surface(a));
	if (plane.IsNull()) {
		return false;
	}

	const TopoDS_Wire a_outer = BRepTools::OuterWire(a);
	if (a_outer.IsNull()) {
		return false;
	}

	const gp_Ax3 f = frame_of(plane);

	configure(eps);
	std::lock_guard<std::mutex> guard(lock());
	// Whether the arcs the fit finds are joined back into whole circles is a
	// caller's choice rather than a tolerance, so it is set here, under the
	// lock, and not in the once-only configure.
	CArea::m_fit_circles = fit_circles;

	PERF("2d area: total");

	// All the operands go into one area rather than being united first. Under
	// the non-zero rule two overlapping tools wind to two, which is still
	// inside, so the overlap is cut once and correctly; under even-odd it
	// would cancel and silently not be cut at all, and the only way to avoid
	// that would be a union per operand.
	CArea tools;
	std::vector<Bnd_Box> operand_boxes;
	bool any_tool = false;

	{
	PERF("2d area: convert operands");
	for (NCollection_List<TopoDS_Shape>::Iterator it(b_input); it.More(); it.Next()) {
		if (it.Value().ShapeType() != TopAbs_FACE || !face_to_area(TopoDS::Face(it.Value()), f, tools)) {
			logger.notice("GEO", 406, "An operand is outside of the 2D engine's domain. Retrying with the kernel.");
			return false;
		}
		any_tool = true;

		operand_boxes.push_back(Bnd_Box());
		BRepBndLib::Add(it.Value(), operand_boxes.back());
		operand_boxes.back().Enlarge(eps);
	}
	}

	if (!any_tool) {
		return false;
	}

	CArea subject;
	{
		PERF("2d area: convert subject");
		CCurve outer;
		if (!wire_to_curve(a_outer, f, true, outer)) {
			logger.notice("GEO", 406, "Operand A's outer boundary is outside of the 2D engine's domain. Retrying with the kernel.");
			return false;
		}
		wind(outer, false);
		subject.append(outer);
	}

	// An inner boundary the lattice cannot represent does not have to cost the
	// whole attempt: if no operand comes near it, the subtraction cannot
	// change it, so it is set aside and put back on the result as it was.
	//
	// Conversion is tried first and this is only the fallback, which matters
	// for more than tidiness. Deciding whether an operand reaches a boundary
	// costs a distance query per operand, and A here is a wall that already
	// has its windows cut, so asking it of every boundary is that many
	// queries times two hundred operands -- which measured at 378 of the
	// wall's 390 seconds, more than the boolean it is replacing.
	std::vector<TopoDS_Wire> carry_through;

	for (TopoDS_Iterator it(a); it.More(); it.Next()) {
		if (it.Value().ShapeType() != TopAbs_WIRE || it.Value().IsSame(a_outer)) {
			continue;
		}

		const TopoDS_Wire& w = TopoDS::Wire(it.Value());

		CCurve inner;
		if (wire_to_curve(w, f, true, inner)) {
			wind(inner, true);
			subject.append(inner);
			continue;
		}

		if (operand_reaches(w, b_input, operand_boxes, eps)) {
			logger.notice("GEO", 406, "Operand A has a boundary an operand reaches that the 2D engine cannot represent. Retrying with the kernel.");
			return false;
		}

		carry_through.push_back(w);
	}

	// Driven here rather than through CArea::Subtract so the result comes back
	// as a PolyTree. A flat list of contours does not say which lies inside
	// which, and re-deriving that afterwards means comparing every contour
	// with every other -- which measured at 30 of the wall's 48 seconds, and
	// at 383 before libarea learned to reject a pair whose bounding boxes are
	// apart. Clipper knows the nesting while it is building the result and
	// hands it over for nothing.
	Clipper2Lib::Paths64 subject_paths, clip_paths;
	for (const CCurve& c : subject.m_curves) {
		subject_paths.push_back(Clipper2Lib::Path64());
		CurveToClipperPath(c, subject_paths.back());
	}
	for (const CCurve& c : tools.m_curves) {
		clip_paths.push_back(Clipper2Lib::Path64());
		CurveToClipperPath(c, clip_paths.back());
	}

	Clipper2Lib::PolyTree64 solution;
	{
		PERF("2d area: subtract");
		// AddSubject and AddClip do not report, unlike Clipper1's AddPaths;
		// Execute is the only thing left that can say no.
		Clipper2Lib::Clipper64 clipper;
		clipper.AddSubject(subject_paths);
		clipper.AddClip(clip_paths);
		if (!clipper.Execute(Clipper2Lib::ClipType::Difference, Clipper2Lib::FillRule::NonZero, solution)) {
			return false;
		}
	}

	// The same reasoning as the material-conservation check the 3D path makes,
	// one dimension down: a cut cannot add area. Neither engine always says
	// when it has broken that -- a contour that comes back read as solid where
	// it should have been read as a hole adds area instead of taking it away,
	// and the total is what shows it. In lattice units on both sides, which is
	// the only place the two are comparable.
	{
		const Clipper2Lib::Paths64 flat = Clipper2Lib::PolyTreeToPaths64(solution);
		if (lattice_area(flat) > lattice_area(subject_paths) * (1. + 1e-9) + 1.) {
			logger.notice("GEO", 405, "2D subtraction returned more area than it started with. Retrying with the kernel.");
			return false;
		}
	}

	PERF("2d area: rebuild");

	NCollection_List<TopoDS_Shape> faces;
	for (const auto& node : solution) {
		if (!faces_from_node(*node, f, eps, faces, carry_through)) {
			logger.notice("GEO", 405, "The 2D result could not be made into a face. Retrying with the kernel.");
			return false;
		}
	}

	if (faces.IsEmpty()) {
		logger.notice("GEO", 404, "2D subtraction left nothing of the first operand. Retrying with the kernel.");
		return false;
	}

	if (faces.Extent() == 1) {
		result = faces.First();
	} else {
		// More than one piece is a legitimate answer: a cut can sever the
		// profile.
		TopoDS_Compound comp;
		BRep_Builder bb;
		bb.MakeCompound(comp);
		for (NCollection_List<TopoDS_Shape>::Iterator it(faces); it.More(); it.Next()) {
			bb.Add(comp, it.Value());
		}
		result = comp;
	}

	return true;
}

#endif
