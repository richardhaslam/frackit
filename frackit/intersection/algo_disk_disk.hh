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
 * \brief Contains the intersection algorithm
 *        between two disks.
 */
#ifndef FRACKIT_DISK_DISK_INTERSECTION_HH
#define FRACKIT_DISK_DISK_INTERSECTION_HH

#include <variant>
#include <stdexcept>

#include <frackit/geometry/point.hh>
#include <frackit/geometry/line.hh>
#include <frackit/geometry/segment.hh>
#include <frackit/geometry/plane.hh>
#include <frackit/geometry/disk.hh>
#include <frackit/precision/precision.hh>

#include "intersectiontraits.hh"
#include "emptyintersection.hh"

#include "algo_segment_segment.hh"
#include "algo_plane_plane.hh"
#include "algo_disk_line.hh"

namespace Frackit {
namespace IntersectionAlgorithms {

//! Intersect two disks
//! The result can be:
//! - a surface bounded by segments and/or elliptical arcs
//! - a segment
//! - a point
//! - no intersection
template<class ctype>
Intersection< Disk<ctype>, Disk<ctype> >
intersect_disk_disk(const Disk<ctype>& disk1,
                    const Disk<ctype>& disk2,
                    ctype eps)
{
    using ResultType  = Intersection< Disk<ctype>, Disk<ctype> >;
    static constexpr int worldDim = 3;

    // first intersect the supporting planes
    std::string isGeomName;
    const auto planeIS = intersect_plane_plane(disk1.supportingPlane(), disk2.supportingPlane(), eps);

    // if the planes don't intersect, the disks don't either
    if (std::holds_alternative<EmptyIntersection<worldDim>>(planeIS))
        return ResultType( EmptyIntersection<worldDim>() );

    // if the result is a line, find the possible segment or point intersection
    else if (std::holds_alternative<Line<ctype, worldDim>>(planeIS))
    {
        const auto& isLine = std::get<Line<ctype, worldDim>>(planeIS);
        const auto is1 = intersect_disk_line(disk1, isLine, eps);
        const auto is2 = intersect_disk_line(disk2, isLine, eps);

        // both disks intersect the support plane of the other disks
        if (std::holds_alternative<Segment<ctype, worldDim>>(is1)
            && std::holds_alternative<Segment<ctype, worldDim>>(is2))
        {
            const auto& seg1 = std::get<Segment<ctype, worldDim>>(is1);
            const auto& seg2 = std::get<Segment<ctype, worldDim>>(is2);
            const auto segmentIS = intersect_segment_segment(seg1, seg2, eps);

            if (std::holds_alternative<Segment<ctype, worldDim>>(segmentIS))
                return ResultType( std::get<Segment<ctype, worldDim>>(segmentIS) );
            if (std::holds_alternative<Point<ctype, worldDim>>(segmentIS))
                return ResultType( std::get<Point<ctype, worldDim>>(segmentIS) );
            else
                return ResultType( EmptyIntersection<worldDim>() );
        }

        // one of the disks might still touch the other disk
        if (std::holds_alternative<Segment<ctype, worldDim>>(is1)
            && std::holds_alternative<Point<ctype, worldDim>>(is2))
        {
            const auto& seg1 = std::get<Segment<ctype, worldDim>>(is1);
            const auto& p2 = std::get<Point<ctype, worldDim>>(is2);
            if (seg1.contains(p2, eps))
                return ResultType( p2 );
        }
        else if (std::holds_alternative<Point<ctype, worldDim>>(is1)
            && std::holds_alternative<Segment<ctype, worldDim>>(is2))
        {
            const auto& p1 = std::get<Point<ctype, worldDim>>(is1);
            const auto& seg2 = std::get<Segment<ctype, worldDim>>(is2);
            if (seg2.contains(p1, eps))
                return ResultType( p1 );
        }

        // intersection is empty
        return ResultType( EmptyIntersection<worldDim>() );
    }

    // if the result is a plane, find the intersection surface
    else if (std::holds_alternative<Plane<ctype, worldDim>>(planeIS))
        throw std::runtime_error(std::string("NotImplemented: planar disk-disk intersections"));

    throw std::runtime_error(std::string("Unexpected plane-plane intersection result"));
}

} // end namespace IntersectionAlgorithms
} // end namespace Frackit

#endif // FRACKIT_DISK_DISK_INTERSECTION_HH
