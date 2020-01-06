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
 *        between a disk and a line.
 */
#ifndef FRACKIT_DISK_LINE_INTERSECTION_HH
#define FRACKIT_DISK_LINE_INTERSECTION_HH

#include <variant>
#include <stdexcept>

#include <frackit/geometry/disk.hh>
#include <frackit/geometry/line.hh>
#include <frackit/geometry/precision.hh>

#include <frackit/occ/breputilities.hh>
#include <frackit/occ/gputilities.hh>

#include "intersectiontraits.hh"
#include "emptyintersection.hh"
#include "algo_plane_line.hh"

namespace Frackit {
namespace IntersectionAlgorithms {

//! Intersect a disk and a line
//! The result can be:
//! - a segment
//! - a point
//! - no intersection
template<class ctype>
Intersection< Disk<ctype>, Line<ctype, 3> >
intersect_disk_line(const Disk<ctype>& disk,
                    const Line<ctype, 3>& line,
                    ctype eps)
{
    using ResultType = Intersection< Disk<ctype>, Line<ctype, 3> >;

    // intersect the line with the plane
    const auto planeLineIS = intersect_plane_line(disk.supportingPlane(), line, eps);

    // empty result
    if (std::holds_alternative<EmptyIntersection<3>>(planeLineIS))
        return ResultType( EmptyIntersection<3>() );

    // potential point intersection
    else if (std::holds_alternative<Point<ctype, 3>>(planeLineIS))
    {
        // check if point is on the disk
        const auto p = std::get<Point<ctype, 3>>(planeLineIS);
        if (disk.contains(p, Precision<ctype>::confusion(), false))
            return ResultType( p );
        else
            return ResultType( EmptyIntersection<3>() );
    }

    // potential segment intersection
    else if (std::holds_alternative<Line<ctype, 3>>(planeLineIS))
    {
        // Intersection might be a segment, a touching point on the rim or empty
        // Use the BRep algorithm package to determine the part of the line on the ellipse.
        // For this, we create a segment of the line in the neighborhood of ellipse which is
        // large enough to cover the entire disk in case they intersect.
        auto source = line.projection(disk.center());
        auto target = source;
        auto dir = Vector<ctype, 3>(line.direction());
        dir *= disk.majorAxisLength()*20.0;
        source += dir;
        target -= dir;

        // find the part of the segment on the disk
        const auto segmentShape = OCCUtilities::makeEdge(source, target);
        const auto diskShape = OCCUtilities::getShape(disk);
        const auto isShape = OCCUtilities::intersect(segmentShape, diskShape, 0.1*eps);
        const auto isEdges = OCCUtilities::getEdges(isShape);

        assert(isEdges.size() <= 1);
        if (isEdges.size() == 1)
        {
            const auto vertices = OCCUtilities::getVertices(isEdges[0]);
            assert(vertices.size() == 2);
            return ResultType( Segment<ctype, 3>(OCCUtilities::point(vertices[0]),
                                                 OCCUtilities::point(vertices[1])) );
        }

        // the line can still touch the ellipse on the rim
        const auto cutShape = OCCUtilities::cut(segmentShape, diskShape, 0.1*eps);
        const auto vertices = OCCUtilities::getVertices(cutShape);

        assert(vertices.size() == 2 || vertices.size() == 4);
        if (vertices.size() == 2)
            return ResultType( EmptyIntersection<3>() );
        else
        {
            // find the new point
            for (const auto& v : vertices)
            {
                const auto curPoint = OCCUtilities::point(v);
                if (!curPoint.isEqual(source, eps) && !curPoint.isEqual(target, eps))
                    return ResultType( curPoint );
            }
        }

        throw std::runtime_error(std::string("Unexpected code behaviour"));
    }
    else
        throw std::runtime_error(std::string("Unexpected Plane-Line intersection result"));
}

} // end namespace IntersectionAlgorithms
} // end namespace Frackit

#endif // FRACKIT_INTERSECT_HH
