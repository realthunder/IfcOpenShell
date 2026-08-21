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

// Lifted out of boolean_utils_2d.cpp when the layer-set path became a second
// caller. Nothing here changed in the move except that wire_to_curve now takes
// whether to close the curve, which used to be unconditional because a contour
// was all it was ever given.

#include "area_2d.h"

#ifdef IFOPSH_WITH_LIBAREA

#include <BRepAdaptor_Curve.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRep_Tool.hxx>
#include <Precision.hxx>
#include <TopoDS.hxx>
#include <gp_Circ.hxx>
#include <gp_Vec.hxx>

#include <cmath>

namespace {
	bool is_left(const Point& a, const Point& b, const Point& p) {
		return (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x) > 0.;
	}

}

namespace ifcopenshell::geom::util::area2d {

	// libarea's Clipper bridge keeps the points of the polygon it is currently
	// building in a file-static list, and its tolerances are static members of
	// CArea. Neither is per-call state, so the whole of the 2D path is taken
	// one at a time. It costs nothing worth measuring: the operations being
	// serialised are the microsecond ones, and the reason this path exists is
	// that the alternative takes seconds.
	std::mutex& lock() {
		static std::mutex m;
		return m;
	}

	// CArea's tolerances are static, so they are set once rather than per
	// call: a second thread in here would otherwise be writing them while the
	// first reads them. They are the same values every time in any case --
	// eps is the model's precision, which does not vary within a run.
	//
	// m_accuracy is the sagitta an arc is allowed when it is walked onto the
	// lattice, and also how far from a circle a run of segments may be and
	// still be recognised as an arc coming back. The two have to agree or
	// nothing that goes in as an arc comes back as one, which is why it is one
	// number and not two.
	void configure(double eps) {
		static std::once_flag once;
		std::call_once(once, [eps]() {
			CArea::m_units = 1.0;
			CArea::m_fit_arcs = true;
			CArea::m_accuracy = eps > 0. ? eps : 1e-5;
			CArea::m_clipper_clean_distance = 0.;
			CArea::m_clipper_simple = false;
		});
	}

	Point to_2d(const gp_Pnt& p, const gp_Ax3& f) {
		const gp_Vec d(f.Location(), p);
		return Point(d.Dot(gp_Vec(f.XDirection())), d.Dot(gp_Vec(f.YDirection())));
	}

	gp_Pnt to_3d(const Point& p, const gp_Ax3& f) {
		return f.Location().Translated(gp_Vec(f.XDirection()) * p.x + gp_Vec(f.YDirection()) * p.y);
	}

	// A line or an arc goes onto the lattice without losing its identity.
	// Anything else -- a B-spline, an ellipse -- would have to be walked into
	// segments and would come back as segments, which is the outcome this
	// whole path exists to avoid, so it is left to the kernel instead.
	bool edge_is_representable(const TopoDS_Edge& e) {
		const GeomAbs_CurveType t = BRepAdaptor_Curve(e).GetType();
		return t == GeomAbs_Line || t == GeomAbs_Circle;
	}

	bool wire_to_curve(const TopoDS_Wire& w, const gp_Ax3& f, bool close, CCurve& curve) {
		BRepTools_WireExplorer xp(w);
		if (!xp.More()) {
			return false;
		}

		curve.append(CVertex(to_2d(BRep_Tool::Pnt(xp.CurrentVertex()), f)));

		for (; xp.More(); xp.Next()) {
			const TopoDS_Edge& edge = xp.Current();
			BRepAdaptor_Curve c(edge);
			const bool reversed = edge.Orientation() == TopAbs_REVERSED;
			const double first = c.FirstParameter();
			const double last = c.LastParameter();
			const gp_Pnt end = c.Value(reversed ? first : last);

			switch (c.GetType()) {
			case GeomAbs_Line:
				curve.append(CVertex(to_2d(end, f)));
				break;
			case GeomAbs_Circle: {
				const gp_Circ circle = c.Circle();
				// Which way the arc turns in this frame, not in the frame the
				// circle was defined in.
				int type = circle.Axis().Direction().Dot(f.Direction()) < 0. ? -1 : 1;
				if (reversed) {
					type = -type;
				}
				const Point centre = to_2d(circle.Location(), f);
				// A CVertex arc is defined by its two end points and its
				// centre, which does not distinguish an arc from its
				// complement once it is more than half a circle. Split it.
				if (std::fabs(last - first) > M_PI) {
					curve.append(CVertex(type, to_2d(c.Value(first + (last - first) * 0.5), f), centre));
				}
				curve.append(CVertex(type, to_2d(end, f), centre));
				break;
			}
			default:
				return false;
			}
		}

		if (close && !curve.IsClosed()) {
			curve.append(curve.m_vertices.front());
		}

		return curve.m_vertices.size() >= (close ? 3u : 2u);
	}

	// FitArcs picks a centre that fits a run of lattice points, so the two end
	// points of the arc it produces are not exactly the same distance from it.
	// Move the centre onto the perpendicular bisector of the chord, on the
	// side it was already on, at the averaged radius.
	void reconcile_arc(const Point& p0, const Point& p1, Point& centre, double& r) {
		const double r0 = p0.dist(centre);
		const double r1 = p1.dist(centre);
		r = (r0 + r1) * 0.5;

		if (std::fabs(r0 - r1) <= Precision::Confusion()) {
			return;
		}

		const double d = p0.dist(p1);
		if (d <= Precision::Confusion()) {
			return;
		}

		double rr = r * r;
		const double dd = d * d * 0.25;
		if (rr < dd) {
			r = d * 0.5;
			rr = dd;
		}

		const double q = std::sqrt(rr - dd);
		const double x = (p0.x + p1.x) * 0.5;
		const double y = (p0.y + p1.y) * 0.5;
		const double dx = q * (p0.y - p1.y) / d;
		const double dy = q * (p1.x - p0.x) / d;

		Point moved(x + dx, y + dy);
		if (is_left(p0, p1, centre) != is_left(p0, p1, moved)) {
			moved = Point(x - dx, y - dy);
		}
		centre = moved;
	}

}

#endif
