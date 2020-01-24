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
 * \brief Contains functionality to compute the
 *        intersections between geometries.
 */
#ifndef FRACKIT_INTERSECT_HH
#define FRACKIT_INTERSECT_HH

#include <cmath>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shell.hxx>

// Geometries of intersections
#include <frackit/geometry/line.hh>
#include <frackit/geometry/segment.hh>
#include <frackit/geometry/plane.hh>
#include <frackit/geometry/disk.hh>
#include <frackit/geometry/quadrilateral.hh>
#include <frackit/geometry/cylindersurface.hh>
#include <frackit/precision/defaultepsilon.hh>

#include "intersectiontraits.hh"
#include "algo_segment_segment.hh"
#include "algo_plane_plane.hh"
#include "algo_plane_line.hh"
#include "algo_disk_line.hh"
#include "algo_quadrilateral_line.hh"
#include "algo_quadrilateral_quadrilateral.hh"
#include "algo_quadrilateral_disk.hh"
#include "algo_disk_disk.hh"
#include "algo_cylsurface_disk.hh"
#include "algo_cylsurface_quadrilateral.hh"
#include "algo_shell_disk.hh"
#include "algo_face_disk.hh"

namespace Frackit {

/*!
 * \brief Interface for intersecting two geometries.
 * \param geo1 The first geometry
 * \param geo2 The second geometry
 * \param eps Tolerance to be used for floating point comparison
 */
template<class Geom1, class Geom2>
EmptyIntersection<0> intersect(const Geom1& geo1,
                               const Geom2& geo2,
                               typename Geom1::ctype eps)
{
    std::string msg = "Intersection algorithm between \"";
    msg += geo1.name();
    msg += "\" and \"";
    msg += geo2.name();
    msg += "\" not implemented";
    throw std::runtime_error(msg);
}

/*!
 * \brief Interface for intersecting two geometries.
 * \param geo1 The first geometry
 * \param geo2 The second geometry
 * \note This overload selects a default tolerance
 */
template<class Geom1, class Geom2>
auto intersect(const Geom1& geo1, const Geom2& geo2)
{
    using std::min;
    const auto eps1 = defaultEpsilon(geo1);
    const auto eps2 = defaultEpsilon(geo2);
    return intersect(geo1, geo2, min(eps1, eps2));
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
 * \brief Intersect a disk and a line.
 * \param disk The disk
 * \param line The line
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Disk<ctype>, Line<ctype, 3> >
intersect(const Disk<ctype>& disk, const Line<ctype, 3>& line, ctype eps)
{ return IntersectionAlgorithms::intersect_disk_line(disk, line, eps); }

/*!
 * \brief Intersect a line and a disk.
 * \param line The line
 * \param disk The disk
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Line<ctype, 3>, Disk<ctype> >
intersect(const Line<ctype, 3>& line, const Disk<ctype>& disk, ctype eps)
{ return intersect(disk, line, eps); }

/*!
 * \brief Intersect a three-dimensional quadrilateral and a line.
 * \param quad The quadrilateral
 * \param line The line
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Quadrilateral<ctype, 3>, Line<ctype, 3> >
intersect(const Quadrilateral<ctype, 3>& quad, const Line<ctype, 3>& line, ctype eps)
{ return IntersectionAlgorithms::intersect_quadrilateral_line(quad, line, eps); }


/*!
 * \brief Intersect a line and a three-dimensional quadrilateral.
 * \param line The line
 * \param quad The quadrilateral
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Line<ctype, 3>, Quadrilateral<ctype, 3>>
intersect(const Line<ctype, 3>& line, const Quadrilateral<ctype, 3>& quad, ctype eps)
{ return intersect(quad, line, eps); }

/*!
 * \brief Intersect two disks.
 * \param disk1 The first disk
 * \param disk2 The second disk
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Disk<ctype>, Disk<ctype> >
intersect(const Disk<ctype>& disk1, const Disk<ctype>& disk2, ctype eps)
{ return IntersectionAlgorithms::intersect_disk_disk(disk1, disk2, eps); }

/*!
 * \brief Intersect two quadrilaterals in 3d space.
 * \param quad1 The first quadrilateral
 * \param quad2 The second quadrilateral
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Quadrilateral<ctype, 3>, Quadrilateral<ctype, 3> >
intersect(const Quadrilateral<ctype, 3>& quad1, const Quadrilateral<ctype, 3>& quad2, ctype eps)
{ return IntersectionAlgorithms::intersect_quadrilateral_quadrilateral(quad1, quad2, eps); }

/*!
 * \brief Intersect a quadrilateral and a disk in 3d space.
 * \param quad The quadrilateral
 * \param disk The disk
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Quadrilateral<ctype, 3>, Disk<ctype> >
intersect(const Quadrilateral<ctype, 3>& quad, const Disk<ctype>& disk, ctype eps)
{ return IntersectionAlgorithms::intersect_quadrilateral_disk(quad, disk, eps); }

/*!
 * \brief Intersect a disk and a quadrilateral in 3d space.
 * \param disk The disk
 * \param quad The quadrilateral
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Disk<ctype>, Quadrilateral<ctype, 3> >
intersect(const Disk<ctype>& disk, const Quadrilateral<ctype, 3>& quad, ctype eps)
{ return intersect(quad, disk, eps); }

/*!
 * \brief Intersect a lateral cylinder surface and a disk.
 * \param cylSurface The lateral cylinder surface
 * \param disk The disk
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< CylinderSurface<ctype>, Disk<ctype> >
intersect(const CylinderSurface<ctype>& cylSurface, const Disk<ctype>& disk, ctype eps)
{ return IntersectionAlgorithms::intersect_cylinderSurface_disk(cylSurface, disk, eps); }

/*!
 * \brief Intersect a disk and a lateral cylinder surface.
 * \param disk The disk
 * \param cylSurface The lateral cylinder surface
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Disk<ctype>, CylinderSurface<ctype> >
intersect(const Disk<ctype>& disk, const CylinderSurface<ctype>& cylSurface, ctype eps)
{ return intersect(cylSurface, disk, eps); }

/*!
 * \brief Intersect a lateral cylinder surface and a quadrilateral.
 * \param cylSurface The lateral cylinder surface
 * \param quad The quadrilateral
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< CylinderSurface<ctype>, Quadrilateral<ctype, 3> >
intersect(const CylinderSurface<ctype>& cylSurface, const Quadrilateral<ctype, 3>& quad, ctype eps)
{ return IntersectionAlgorithms::intersect_cylinderSurface_quadrilateral(cylSurface, quad, eps); }

/*!
 * \brief Intersect a quadrilateral and a lateral cylinder surface.
 * \param quad The quadrilateral
 * \param cylSurface The lateral cylinder surface
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Quadrilateral<ctype, 3>, CylinderSurface<ctype> >
intersect(const Quadrilateral<ctype, 3>& quad, const CylinderSurface<ctype>& cylSurface, ctype eps)
{ return intersect(cylSurface, quad, eps); }

/*!
 * \brief Intersect a disk and the boundary (TopoDS_Shell) of a solid.
 * \param disk The disk
 * \param shell The shell of a solid
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Disk<ctype>, TopoDS_Shell >
intersect(const Disk<ctype>& disk, const TopoDS_Shell& shell, ctype eps)
{ return IntersectionAlgorithms::intersect_shell_disk(shell, disk, eps); }

/*!
 * \brief Intersect the boundary (TopoDS_Shell) of a solid and a disk.
 * \param shell The shell of a solid
 * \param disk The disk
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< TopoDS_Shell, Disk<ctype> >
intersect(const TopoDS_Shell& shell, const Disk<ctype>& disk, ctype eps)
{ return intersect(disk, shell, eps); }

/*!
 * \brief Intersect a disk and a face shape.
 * \param disk The disk
 * \param face The face shape
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Disk<ctype>, TopoDS_Face >
intersect(const Disk<ctype>& disk, const TopoDS_Face& face, ctype eps)
{ return IntersectionAlgorithms::intersect_face_disk(face, disk, eps); }

/*!
 * \brief Intersect a face shape and a disk.
 * \param face The face shape
 * \param disk The disk
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< TopoDS_Face, Disk<ctype> >
intersect(const TopoDS_Face& face, const Disk<ctype>& disk, ctype eps)
{ return intersect(disk, face, eps); }

} // end namespace Frackit

#endif // FRACKIT_INTERSECT_HH
