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
 *        between a face shape and a planar
 *        two-dimensional geometry in 3 space.
 */
#ifndef FRACKIT_FACE_PLANAR_GEOM_INTERSECTION_HH
#define FRACKIT_FACE_PLANAR_GEOM_INTERSECTION_HH

#include <algorithm>
#include <stdexcept>

#include <TopoDS_Vertex.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>

#include <frackit/geometry/point.hh>
#include <frackit/occ/breputilities.hh>
#include <frackit/common/extractctype.hh>

#include <frackit/intersection/intersectiontraits.hh>
#include <frackit/intersection/emptyintersection.hh>
#include "algo_find_touching_points.hh"

namespace Frackit {
namespace IntersectionAlgorithms {

/*!
 * \brief Intersect a face shape and a planar 2d geometry in 3d space.
 *        The result can be:
 *        - points
 *        - edge shapes
 *        - face shapes
 *        Multiple of the above are possible since the face shape might be curved.
 * \param face The (possible curved) face shape
 * \param faceGeom The planar geometry (planar bounded 2d surface in 3d space)
 * \param eps Tolerance to be used for boolean operations
 */
template<class PlanarGeometry>
Intersection< TopoDS_Face, PlanarGeometry >
intersect_face_planarGeometry(const TopoDS_Face& face,
                              const PlanarGeometry& faceGeom,
                              typename CoordinateTypeTraits<PlanarGeometry>::type eps)
{
    using ResultType = Intersection< TopoDS_Face, PlanarGeometry >;
    using ctype = typename CoordinateTypeTraits<PlanarGeometry>::type;

    std::vector<Point<ctype, 3>> pointIntersections;
    std::vector<TopoDS_Edge> edgeIntersections;
    std::vector<TopoDS_Face> faceIntersections;

    // find possible face intersections
    const auto faceGeomShape = OCCUtilities::getShape(faceGeom);
    const auto faceGeomIS = OCCUtilities::intersect(faceGeomShape, face, eps);
    const auto faceGeomISFaces = OCCUtilities::getFaces(faceGeomIS);
    for (auto& f : faceGeomISFaces)
        faceIntersections.emplace_back(std::move(f));

    // find possible edge intersections
    const auto faceGeomCut = OCCUtilities::cut(faceGeomShape, face, eps);
    const auto faceGeomCutFaces = OCCUtilities::getFaces(faceGeomCut);
    for (unsigned int i = 0; i < faceGeomCutFaces.size(); ++i)
    {
        // edges of the fragments coinciding with the face
        const auto cutFaceWires = OCCUtilities::getWires(faceGeomCutFaces[i]);
        for (const auto& wire : cutFaceWires)
        {
            const auto wireIS = OCCUtilities::intersect(wire, face, eps);
            const auto edgesOnFace = OCCUtilities::getEdges(wireIS);
            for (auto& edge : edgesOnFace)
            {
                // only use this edge if it is not on a previously found face
                bool isContained = false;
                for (const auto& faceIS : faceIntersections)
                {
                    const auto edgeFaceIS = OCCUtilities::intersect(edge, faceIS, eps);
                    const auto edgeFaceISEdges = OCCUtilities::getEdges(edgeFaceIS);
                    if (!edgeFaceISEdges.empty()) { isContained = true; break; }
                }

                if (isContained)
                    continue;

                if (std::none_of(edgeIntersections.begin(),
                                 edgeIntersections.end(),
                                 [&] (const auto& e) { return edge.IsSame(e); }))
                    edgeIntersections.emplace_back(edge);
            }
        }
    }

    // Find possible touching points. We might find points
    // again that resulted from intersections in faces or edges.
    // Thus, we have to check for duplicates here and avoid all
    // points that have been detected so far. The newly detected
    // ones are true touching points. Note that this would not be
    // necessary if we restricted the face to be planar. But, on
    // nonplanar and periodic faces you can have all kinds of
    // intersections. To this end, collect all points found so far.
    std::vector<Point<ctype, 3>> curIsPoints;

    for (const auto& faceIS : faceIntersections)
        for (const auto& v : OCCUtilities::getVertices(faceIS))
        {
            const auto p = OCCUtilities::point(v);
            if (std::none_of(curIsPoints.begin(),
                             curIsPoints.end(),
                             [&] (const auto& isP) { return isP.isEqual(p, eps); }))
                curIsPoints.emplace_back(p);
        }

    for (const auto& edgeIS : edgeIntersections)
        for (const auto& v : OCCUtilities::getVertices(edgeIS))
        {
            const auto p = OCCUtilities::point(v);
            if (std::none_of(curIsPoints.begin(),
                             curIsPoints.end(),
                             [&] (const auto& isP) { return isP.isEqual(p, eps); }))
                curIsPoints.emplace_back(p);
        }

    const auto touchPoints = find_touching_points(faceGeomShape, face, eps);
    for (auto&& p : touchPoints)
        if (std::none_of(curIsPoints.begin(),
                         curIsPoints.end(),
                         [&] (const auto& isP) { return isP.isEqual(p, eps); })
            && std::none_of(pointIntersections.begin(),
                            pointIntersections.end(),
                            [&] (const auto& isP) { return isP.isEqual(p, eps); }))
            pointIntersections.emplace_back(std::move(p));

    ResultType result;
    for (auto& p : pointIntersections) result.emplace_back(std::move(p));
    for (auto& e : edgeIntersections) result.emplace_back(std::move(e));
    for (auto& f : faceIntersections) result.emplace_back(std::move(f));

    if (result.empty())
        return ResultType({EmptyIntersection<3>()});
    return result;
}

} // end namespace IntersectionAlgorithms
} // end namespace Frackit

#endif // FRACKIT_FACE_PLANAR_GEOM_INTERSECTION_HH
