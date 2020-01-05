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
 *        between two segments.
 */
#ifndef FRACKIT_SEGMENT_SEGMENT_INTERSECTION_HH
#define FRACKIT_SEGMENT_SEGMENT_INTERSECTION_HH

#include <stdexcept>
#include <vector>

#include <gp_Vec.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopExp_Explorer.hxx>

#include <frackit/geometry/segment.hh>
#include <frackit/geometry/precision.hh>
#include <frackit/common/utilities.hh>

#include <frackit/intersection/intersectiontraits.hh>
#include <frackit/intersection/emptyintersection.hh>

namespace Frackit {
namespace IntersectionAlgorithms {

//! Intersect two segments
//! The result can be:
//! - a segment
//! - a point
//! - no intersection
template<class ctype, int wd>
Intersection< Segment<ctype, wd>, Segment<ctype, wd> >
intersect_segment_segment(const Segment<ctype, wd>& segment1,
                          const Segment<ctype, wd>& segment2,
                          ctype eps)
{
    using ResultType = Intersection<Segment<ctype, wd>, Segment<ctype, wd>>;
    using Segment = Frackit::Segment<ctype, wd>;
    using Direction = typename Segment::Direction;

    const auto& s1 = segment1.source(); const auto& t1 = segment1.target();
    const auto& s2 = segment2.source(); const auto& t2 = segment2.target();

    // if the segments are parallel, there might be an overlap
    const auto dir1 = Direction(Vector<ctype, wd>(s1, t1));
    const auto dir2 = Direction(Vector<ctype, wd>(s2, t2));
    gp_Vec v1(dir1.x(), dir1.y(), dir1.z());
    gp_Vec v2(dir2.x(), dir2.y(), dir2.z());

    if (v1.IsParallel(v2, Precision<ctype>::angular()))
    {
        // check if the lines overlap
        if (!segment1.supportingLine().contains(s2, eps))
            return ResultType( EmptyIntersection<wd>() );

        // find overlap
        const bool s2On1 = segment1.contains(s2, eps, false);
        const bool t2On1 = segment1.contains(t2, eps, false);
        if (s2On1 && t2On1) return ResultType( segment2 );

        const bool s1On2 = segment2.contains(s1, eps, false);
        const bool t1On2 = segment2.contains(t1, eps, false);
        if (s1On2 && t1On2) return ResultType( segment1 );

        if (s2On1 && s1On2) return ResultType( Segment(s2, s1) );
        if (s2On1 && t1On2) return ResultType( Segment(s2, t1) );

        if (t2On1 && s1On2) return ResultType( Segment(t2, s1) );
        if (t2On1 && t1On2) return ResultType( Segment(t2, t1) );
    }

    // Fragment the lines to find the intersection point
    std::vector<TopoDS_Shape> segmentEdges;
    segmentEdges.emplace_back(OCCUtilities::makeEdge(s1, t1));
    segmentEdges.emplace_back(OCCUtilities::makeEdge(s2, t2));

    const auto fragmentShape = OCCUtilities::fragment(segmentEdges, 0.1*eps);
    const auto edges = OCCUtilities::getEdges(fragmentShape);

    // If the result has 3 edges, a corner of one segment splits the other segment.
    // If it has 4 edges, both segments are split by the intersection point.
    // In either case, find the newly created point which is the intersection.
    const auto numEdges = edges.size();
    if (numEdges == 3 || numEdges == 4)
    {
        for (const auto& edge : edges)
        {
            for (TopExp_Explorer explorer(edge, TopAbs_VERTEX); explorer.More(); explorer.Next())
            {
                auto curPoint = OCCUtilities::point(TopoDS::Vertex(explorer.Current()));
                if (!curPoint.isEqual(s1, eps) && !curPoint.isEqual(t1, eps)
                    && !curPoint.isEqual(s2, eps) && !curPoint.isEqual(t2, eps))
                    return ResultType( std::move(curPoint) );
            }
        }

        throw std::runtime_error(std::string("Could not find intersection point"));
    }

    // a pair of corners being equal is the last possible option
    if (s1.isEqual(s2, eps)) return ResultType( s1 );
    else if (s1.isEqual(t2, eps)) return ResultType( s1 );
    else if (t1.isEqual(s2, eps)) return ResultType( t1 );
    else if (t1.isEqual(t2, eps)) return ResultType( t1 );

    return ResultType( EmptyIntersection<wd>() );
}

} // end namespace IntersectionAlgorithms
} // end namespace Frackit

#endif // FRACKIT_SEGMENT_SEGMENT_INTERSECTION_HH
