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

// Geometries of intersections
#include <frackit/geometry/point.hh>
#include <frackit/geometry/line.hh>
#include <frackit/geometry/segment.hh>
#include <frackit/geometry/direction.hh>
#include <frackit/geometry/plane.hh>
#include <frackit/geometry/disk.hh>
#include <frackit/geometry/cylindersurface.hh>
#include <frackit/geometry/precision.hh>
#include <frackit/magnitude/length.hh>

#include "intersectiontraits.hh"
#include "algo_segment_segment.hh"
#include "algo_plane_plane.hh"
#include "algo_plane_line.hh"
#include "algo_disk_line.hh"
#include "algo_disk_disk.hh"
#include "algo_cylsurface_disk.hh"

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
    const auto eps = Precision<ctype>::confusion()
                     *min(computeLength(segment1), computeLength(segment2));
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
Intersection< CylinderSurface<ctype>, Disk<ctype> >
intersect(const CylinderSurface<ctype>& cylSurface, const Disk<ctype>& disk, ctype eps)
{ return IntersectionAlgorithms::intersect_cylinderSurface_disk(cylSurface, disk, eps); }

/*!
 * \brief Intersect a lateral cylinder surface and a disk.
 * \param cylSurface The lateral cylinder surface
 * \param disk The disk
 * \note This overload uses the default tolerance
 */
template<class ctype>
Intersection< CylinderSurface<ctype>, Disk<ctype> >
intersect(const CylinderSurface<ctype>& cylSurface, const Disk<ctype>& disk)
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
Intersection< Disk<ctype>, CylinderSurface<ctype> >
intersect(const Disk<ctype>& disk, const CylinderSurface<ctype>& cylSurface, ctype eps)
{ return intersect(cylSurface, disk, eps); }

/*!
 * \brief Intersect a disk and a lateral cylinder surface.
 * \param disk The disk
 * \param cylSurface The lateral cylinder surface
 * \note This overload uses the default tolerance
 */
template<class ctype>
Intersection< Disk<ctype>, CylinderSurface<ctype> >
intersect(const Disk<ctype>& disk, const CylinderSurface<ctype>& cylSurface)
{ return intersect(cylSurface, disk); }

} // end namespace Frackit

#endif // FRACKIT_INTERSECT_HH
