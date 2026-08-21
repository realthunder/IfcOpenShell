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

// The plumbing shared by everything in this kernel that reaches libarea. Its
// state -- the point list its Clipper bridge builds into, and the tolerances
// that are static members of CArea -- belongs to the process rather than to a
// call, so there is one lock over it and one place it is configured, however
// many callers there are.

#ifndef AREA_2D_H
#define AREA_2D_H

#ifdef IFOPSH_WITH_LIBAREA

#include <libarea/Area.h>

#include <TopoDS_Edge.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Ax3.hxx>
#include <gp_Pnt.hxx>

#include <mutex>

namespace ifcopenshell::geom::util::area2d {

	std::mutex& lock();
	void configure(double eps);

	Point to_2d(const gp_Pnt& p, const gp_Ax3& f);
	gp_Pnt to_3d(const Point& p, const gp_Ax3& f);

	bool edge_is_representable(const TopoDS_Edge& e);

	// close: append the first vertex again, which a contour needs and an open
	// curve -- a wall axis, say -- must not have.
	bool wire_to_curve(const TopoDS_Wire& w, const gp_Ax3& f, bool close, CCurve& curve);

	void reconcile_arc(const Point& p0, const Point& p1, Point& centre, double& r);

}

#endif

#endif
