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
 *        between a planar two-dimensional geometry
 *        and a line in three-dimensional space.
 */
#ifndef FRACKIT_PLANAR_GEOMETRY_LINE_INTERSECTION_HH
#define FRACKIT_PLANAR_GEOMETRY_LINE_INTERSECTION_HH

#include <variant>
#include <stdexcept>

#include <frackit/geometry/line.hh>
#include <frackit/common/extractdimension.hh>
#include <frackit/occ/breputilities.hh>
#include <frackit/precision/precision.hh>

#include "intersectiontraits.hh"
#include "emptyintersection.hh"
#include "algo_plane_line.hh"

namespace Frackit {
namespace IntersectionAlgorithms {

/*!
 * \brief Intersect a planar face and a line
 *        The result can be:
 *        - a segment
 *        - a point
 *        - no intersection
 * \param faceGeom The planar geometry
 * \param line The line
 * \param charLength A characteristic length scale of the face
 * \param containsEps Tolerance to be used for contains() queries on the planar geometry
 * \param eps Tolerance to be used for boolean operations
 */
template<class PlanarGeometry, class ctype>
Intersection< PlanarGeometry, Line<ctype, 3> >
intersect_planarGeometry_line(const PlanarGeometry& faceGeom,
                              const Line<ctype, 3>& line,
                              ctype charLength,
                              ctype containsEps,
                              ctype eps)
{
    static_assert( (DimensionalityTraits<PlanarGeometry>::geometryDimension() == 2 &&
                    DimensionalityTraits<PlanarGeometry>::worldDimension() == 3),
                  "This algorithm expects planar, two-dimensional geometries in 3d space");

    using ResultType = Intersection< PlanarGeometry, Line<ctype, 3> >;

    // intersect the line with the plane
    const auto planeLineIS = intersect_plane_line(faceGeom.supportingPlane(), line, eps);

    // empty result
    if (std::holds_alternative<EmptyIntersection<3>>(planeLineIS))
        return ResultType( EmptyIntersection<3>() );

    // potential point intersection
    else if (std::holds_alternative<Point<ctype, 3>>(planeLineIS))
    {
        // check if point is on the face
        const auto p = std::get<Point<ctype, 3>>(planeLineIS);
        if (faceGeom.contains(p, containsEps, false))
            return ResultType( p );
        else
            return ResultType( EmptyIntersection<3>() );
    }

    // potential segment intersection
    else if (std::holds_alternative<Line<ctype, 3>>(planeLineIS))
    {
        // Intersection might be a segment, a touching point on the rim of the face, or empty.
        // Use the BRep algorithm package to determine the part of the line on the face.
        // For this, we create a segment of the line in the neighborhood of the face which is
        // large enough to cover the entire face in case they intersect.
        auto source = line.projection(faceGeom.center());
        auto target = source;
        auto dir = Vector<ctype, 3>(line.direction());
        dir *= charLength*20.0;
        source += dir;
        target -= dir;

        // find the part of the segment on the disk
        const auto segmentShape = OCCUtilities::makeEdge(source, target);
        const auto faceShape = OCCUtilities::getShape(faceGeom);
        const auto isShape = OCCUtilities::intersect(segmentShape, faceShape, 0.1*eps);
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
        const auto cutShape = OCCUtilities::cut(segmentShape, faceShape, 0.1*eps);
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

        throw std::runtime_error("Unexpected code behaviour");
    }

    throw std::runtime_error("Unexpected Plane-Line intersection result");
}

} // end namespace IntersectionAlgorithms
} // end namespace Frackit

#endif // FRACKIT_PLANAR_GEOMETRY_LINE_INTERSECTION_HH
