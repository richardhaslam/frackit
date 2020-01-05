// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief \todo TODO doc me.
 */
#ifndef FRACKIT_INTERSECT_HH
#define FRACKIT_INTERSECT_HH

#include <algorithm>
#include <type_traits>
#include <stdexcept>
#include <cassert>
#include <utility>
#include <variant>
#include <vector>
#include <cmath>

// some algorithms/classes receive handles instead of objects
#include <Standard_Handle.hxx>

// classes from the geometric processors package
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <gp_Dir.hxx>
#include <gp_Ax2.hxx>
#include <gp_Elips.hxx>

// classes from the geometry package
#include <GeomAPI_IntSS.hxx>
#include <GeomAPI_IntCS.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <Geom_Plane.hxx>
#include <Geom_Surface.hxx>
#include <Geom_Line.hxx>
#include <Geom_Curve.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <Geom_Ellipse.hxx>

// builders for TopoDS_Shapes
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRep_Tool.hxx>

// BRep primitives and operations
#include <BRepPrimAPI_MakeCylinder.hxx>

// shapes to be passed to intersection algorithms
#include <TopTools_ListOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Vertex.hxx>

// class to explore the resulting shape of an intersection
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopAbs_Orientation.hxx>

// Algorithm for boolean operations on shapes
#include <BRepAlgoAPI_BuilderAlgo.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Cut.hxx>

// Geometries of intersections
#include <frackit/geometry/point.hh>
#include <frackit/geometry/line.hh>
#include <frackit/geometry/segment.hh>
#include <frackit/geometry/direction.hh>
#include <frackit/geometry/plane.hh>
#include <frackit/geometry/disk.hh>
#include <frackit/geometry/cylindricalsurface.hh>

// base tolerances for floating point comparisons
#include <frackit/geometry/precision.hh>

// utility functionality
#include <frackit/common/utilities.hh>
#include <frackit/common/math.hh>

#include "intersectiontraits.hh"
#include "algorithms/intersection_segment_segment.hh"
#include "algorithms/intersection_plane_plane.hh"
#include "algorithms/intersection_plane_line.hh"
#include "algorithms/intersection_disk_line.hh"
#include "algorithms/intersection_disk_disk.hh"
#include "algorithms/intersection_cylsurface_disk.hh"

namespace Frackit {

/*!
 * \brief Interface for intersecting two geometries.
 * \param ge1 The first geometry
 * \param geo2 The second geometry
 * \param eps Tolerance to be used for floating point comparison (optional)
 */
template<class Geom1, class Geom2, class ctype = double>
EmptyIntersection<0> intersect(const Geom1& geo1, const Geom2& geo2, ctype eps = 0.0)
{
    std::string msg = "Intersection algorithm between \"";
    msg += geo1.name();
    msg += "\" and \"";
    msg += geo2.name();
    msg += "\" not implemented";
    throw std::runtime_error(msg);
}

/*!
 * \brief Intersect two segments.
 * \param segment1 The first segment
 * \param segment2 The second segment
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype, int wd>
Intersection< Segment<ctype, wd>, Segment<ctype, wd> >
intersect(const Segment<ctype, wd>& segment1, const Segment<ctype, wd>& segment2, ctype eps)
{ return IntersectionAlgorithms::intersect_segment_segment(segment1, segment2, eps); }

/*!
 * \brief Intersect two segments.
 * \param segment1 The first segment
 * \param segment2 The second segment
 * \note This overload uses the default tolerance
 */
template<class ctype, int wd>
Intersection< Segment<ctype, wd>, Segment<ctype, wd> >
intersect(const Segment<ctype, wd>& segment1, const Segment<ctype, wd>& segment2)
{
    using std::min;
    const auto eps = Precision<ctype>::confusion()*min(segment1.length(), segment2.length());
    return intersect(segment1, segment2, eps);
}

/*!
 * \brief Intersect two planes.
 * \param plane1 The first plane
 * \param plane2 The second plane
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Plane<ctype, 3>, Plane<ctype, 3> >
intersect(const Plane<ctype, 3>& plane1, const Plane<ctype, 3>& plane2, ctype eps)
{ return IntersectionAlgorithms::intersect_plane_plane(plane1, plane2, eps); }

/*!
 * \brief Intersect two planes.
 * \param plane1 The first plane
 * \param plane2 The second plane
 * \note This overload uses the default tolerance
 */
template<class ctype>
Intersection< Plane<ctype, 3>, Plane<ctype, 3> >
intersect(const Plane<ctype, 3>& plane1, const Plane<ctype, 3>& plane2)
{ return intersect(plane1, plane2, Precision<ctype>::confusion()); }

/*!
 * \brief Intersect a plane and a line.
 * \param plane The plane
 * \param line The line
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Plane<ctype, 3>, Line<ctype, 3> >
intersect(const Plane<ctype, 3>& plane, const Line<ctype, 3>& line, ctype eps)
{ return IntersectionAlgorithms::intersect_plane_line(plane, line, eps); }

/*!
 * \brief Intersect a plane and a line.
 * \param plane The plane
 * \param line The line
 * \note This overload uses the default tolerance
 */
template<class ctype>
Intersection< Plane<ctype, 3>, Line<ctype, 3> >
intersect(const Plane<ctype, 3>& plane, const Line<ctype, 3>& line)
{ return intersect(plane, line, Precision<ctype>::confusion()); }

/*!
 * \brief Intersect a line and a plane.
 * \param line The plane
 * \param plane The line
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Line<ctype, 3>, Plane<ctype, 3> >
intersect(const Line<ctype, 3>& line, const Plane<ctype, 3>& plane, ctype eps)
{ return intersect(plane, line, eps); }

/*!
 * \brief Intersect a line and a plane.
 * \param line The plane
 * \param plane The line
 * \note This overload uses the default tolerance
 */
template<class ctype>
Intersection< Line<ctype, 3>, Plane<ctype, 3> >
intersect(const Line<ctype, 3>& line, const Plane<ctype, 3>& plane)
{ return intersect(plane, line); }

/*!
 * \brief Intersect a disk and a line.
 * \param disk The disk
 * \param plane The line
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Disk<ctype>, Line<ctype, 3> >
intersect(const Disk<ctype>& disk, const Line<ctype, 3>& line, ctype eps)
{ return IntersectionAlgorithms::intersect_disk_line(disk, line, eps); }

/*!
 * \brief Intersect a disk and a line.
 * \param disk The disk
 * \param plane The line
 * \note This overload uses the default tolerance
 */
template<class ctype>
Intersection< Disk<ctype>, Line<ctype, 3> >
intersect(const Disk<ctype>& disk, const Line<ctype, 3>& line)
{ return intersect(disk, line, Precision<ctype>::confusion()*disk.minorAxisLength()); }

/*!
 * \brief Intersect a line and a disk.
 * \param plane The line
 * \param disk The disk
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Line<ctype, 3>, Disk<ctype> >
intersect(const Line<ctype, 3>& line, const Disk<ctype>& disk, ctype eps)
{ return intersect(disk, line, eps); }

/*!
 * \brief Intersect a line and a disk.
 * \param plane The line
 * \param disk The disk
 * \note This overload uses the default tolerance
 */
template<class ctype>
Intersection< Line<ctype, 3>, Disk<ctype> >
intersect(const Line<ctype, 3>& line, const Disk<ctype>& disk)
{ return intersect(disk, line); }

/*!
 * \brief Intersect two disks.
 * \param disk1 The first disk
 * \param disk2 The second disk
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Disk<ctype>, Disk<ctype> >
intersect(const Disk<ctype>& disk1,
          const Disk<ctype>& disk2,
          ctype eps)
{ return IntersectionAlgorithms::intersect_disk_disk(disk1, disk2, eps); }

/*!
 * \brief Intersect two disks.
 * \param disk1 The first disk
 * \param disk2 The second disk
 * \note This overload uses the default tolerance
 */
template<class ctype>
Intersection< Disk<ctype>, Disk<ctype> >
intersect(const Disk<ctype>& disk1, const Disk<ctype>& disk2)
{
    using std::min;
    auto eps = min(disk1.minorAxisLength(), disk2.minorAxisLength());
    eps *= Precision<ctype>::confusion();
    return intersect(disk1, disk2, eps);
}

/*!
 * \brief Intersect a lateral cylinder surface and a disk.
 * \param cylSurface The lateral cylinder surface
 * \param disk The disk
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< CylindricalSurface<ctype>, Disk<ctype> >
intersect(const CylindricalSurface<ctype>& cylSurface, const Disk<ctype>& disk, ctype eps)
{ return IntersectionAlgorithms::intersect_cylinderSurface_disk(cylSurface, disk, eps); }

/*!
 * \brief Intersect a lateral cylinder surface and a disk.
 * \param cylSurface The lateral cylinder surface
 * \param disk The disk
 * \note This overload uses the default tolerance
 */
template<class ctype>
Intersection< CylindricalSurface<ctype>, Disk<ctype> >
intersect(const CylindricalSurface<ctype>& cylSurface, const Disk<ctype>& disk)
{
    using std::min;
    auto eps = min(disk.minorAxisLength(), cylSurface.radius());
    eps = min(eps, cylSurface.height());
    eps *= Precision<ctype>::confusion();
    return intersect(cylSurface, disk, eps);
}

/*!
 * \brief Intersect a disk and a lateral cylinder surface.
 * \param disk The disk
 * \param cylSurface The lateral cylinder surface
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Disk<ctype>, CylindricalSurface<ctype> >
intersect(const Disk<ctype>& disk, const CylindricalSurface<ctype>& cylSurface, ctype eps)
{ return intersect(cylSurface, disk, eps); }

/*!
 * \brief Intersect a disk and a lateral cylinder surface.
 * \param disk The disk
 * \param cylSurface The lateral cylinder surface
 * \note This overload uses the default tolerance
 */
template<class ctype>
Intersection< Disk<ctype>, CylindricalSurface<ctype> >
intersect(const Disk<ctype>& disk, const CylindricalSurface<ctype>& cylSurface)
{ return intersect(cylSurface, disk); }

} // end namespace Frackit

#endif // FRACKIT_INTERSECT_HH
