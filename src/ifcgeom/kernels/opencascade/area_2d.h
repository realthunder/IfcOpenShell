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

// libarea's kurve/geometry.h, reached through Area.h, defines a block of
// unprefixed macros -- NONE, PI, CW, ACW, LINEAR, ORIGIN, TANGENT and the
// rest -- and an installed header that does that poisons every translation
// unit downstream of it. rocksdb/db.h is the one that found us:
//
//     enum class SizeApproximationFlags : uint8_t { NONE = 0, ... };
//
// becomes "0 = 0" as soon as this header is included first, so
// layerset_offset.cpp stopped compiling the moment the build actually found
// libarea. Nothing here uses any of them (checked, not assumed), so they are
// dropped at the door rather than left for the next includer to trip over.
// The real fix belongs upstream in libarea, where NONE in particular is dead
// code, but this header is the boundary we control.
#undef PI
#undef DegreesToRadians
#undef RadiansToDegrees
#undef NEARLY_ONE
#undef CPTANGENTTOL
#undef TANTO
#undef ANTITANTO
#undef TANGENT
#undef NEARINT
#undef FARINT
#undef LEFTINT
#undef RIGHTINT
#undef CFILLET
#undef CHAMFER
#undef GEOFF_LEFT
#undef NONE
#undef GEOFF_RIGHT
#undef LINEAR
#undef ACW
#undef CW
#undef INVALID_POINT
#undef INVALID_POINT3D
#undef INVALID_CLINE
#undef INVALID_CIRCLE
#undef ORIGIN
#undef NULL_VECTOR
#undef Z_VECTOR
#undef Y_VECTOR
#undef X_VECTOR
#undef HORIZ_CLINE

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
