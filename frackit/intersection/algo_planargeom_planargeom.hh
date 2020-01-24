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
 *        between two planar two-dimensional geometries
 *        in three-dimensional space.
 */
#ifndef FRACKIT_PLANARGEOM_PLANARGEOM_INTERSECTION_HH
#define FRACKIT_PLANARGEOM_PLANARGEOM_INTERSECTION_HH

#include <variant>
#include <stdexcept>

#include <frackit/geometry/point.hh>
#include <frackit/geometry/line.hh>
#include <frackit/geometry/segment.hh>
#include <frackit/geometry/plane.hh>

#include "intersectiontraits.hh"
#include "emptyintersection.hh"

#include "algo_segment_segment.hh"
#include "algo_plane_plane.hh"
#include "algo_planargeom_line.hh"

namespace Frackit {
namespace IntersectionAlgorithms {

/*!
 * \brief Intersect two planar faces
 *        The result can be:
 *        - a face
 *        - a segment
 *        - a point
 *        - no intersection
 * \param faceGeom1 The first planar geometry
 * \param faceGeom2 The second planar geometry
 * \param charLength A characteristic length scale to be used
 * \param containsEps1 Tolerance to be used for contains() queries on the first planar geometry
 * \param containsEps2 Tolerance to be used for contains() queries on the second planar geometry
 * \param eps Tolerance to be used for boolean operations
 */
template<class PlanarGeom1, class PlanarGeom2, class ctype>
Intersection< PlanarGeom1, PlanarGeom2 >
intersect_planargeometry_planargeometry(const PlanarGeom1& faceGeom1,
                                        const PlanarGeom2& faceGeom2,
                                        ctype charLength,
                                        ctype containsEps1,
                                        ctype containsEps2,
                                        ctype eps)
{
    static_assert( (DimensionalityTraits<PlanarGeom1>::geometryDimension() == 2 &&
                    DimensionalityTraits<PlanarGeom1>::worldDimension() == 3),
                   "This algorithm expects planar, two-dimensional geometries in 3d space");
    static_assert( (DimensionalityTraits<PlanarGeom2>::geometryDimension() == 2 &&
                    DimensionalityTraits<PlanarGeom2>::worldDimension() == 3),
                   "This algorithm expects planar, two-dimensional geometries in 3d space");

    using ResultType  = Intersection< PlanarGeom1, PlanarGeom2 >;
    static constexpr int worldDim = 3;

    // first intersect the supporting planes
    const auto planeIS = intersect_plane_plane(faceGeom1.supportingPlane(), faceGeom2.supportingPlane(), eps);

    // if the planes don't intersect, the geometries don't either
    if (std::holds_alternative<EmptyIntersection<worldDim>>(planeIS))
        return ResultType( EmptyIntersection<worldDim>() );

    // if the result is a line, find the possible segment or point intersection
    else if (std::holds_alternative<Line<ctype, worldDim>>(planeIS))
    {
        const auto& isLine = std::get<Line<ctype, worldDim>>(planeIS);
        const auto is1 = intersect_planargeometry_line(faceGeom1, isLine, charLength, containsEps1, eps);
        const auto is2 = intersect_planargeometry_line(faceGeom2, isLine, charLength, containsEps2, eps);

        // each faces intersect the support plane of the other face
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

        // one of the faces might still touch the other face
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
        throw std::runtime_error("NotImplemented: two-dimensional intersections between two faces");

    throw std::runtime_error("Unexpected plane-plane intersection result");
}

} // end namespace IntersectionAlgorithms
} // end namespace Frackit

#endif // FRACKIT_PLANARGEOM_PLANARGEOM_INTERSECTION_HH
