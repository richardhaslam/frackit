// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \brief Contains the intersection algorithm between two face shapes in 3d space.
 */
#ifndef FRACKIT_FACE_FACE_FACE_3D_INTERSECTION_HH
#define FRACKIT_FACE_FACE_FACE_3D_INTERSECTION_HH

#include <algorithm>
#include <stdexcept>

#include <TopoDS_Vertex.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>

#include <frackit/geometry/ctype.hh>
#include <frackit/geometry/point.hh>
#include <frackit/occ/breputilities.hh>

#include <frackit/intersection/intersectiontraits.hh>
#include <frackit/intersection/emptyintersection.hh>

namespace Frackit {
namespace IntersectionAlgorithms {

/*!
 * \brief Intersect two face shapes in 3d space.
 *        The result can be:
 *        - points
 *        - edge shapes
 *        - face shapes
 *        Multiple of the above are possible since the faces shape might be curved.
 * \param face1 The first face shape
 * \param face2 The seond face shape
 * \param eps Tolerance to be used for boolean operations
 */
Intersection< TopoDS_Face, TopoDS_Face >
intersect_face_face_3d(const TopoDS_Face& face1,
                       const TopoDS_Face& face2,
                       typename CoordinateTypeTraits<TopoDS_Face>::type eps)
{
    using ResultType = Intersection< TopoDS_Face, TopoDS_Face >;
    using ctype = typename CoordinateTypeTraits<TopoDS_Face>::type;
    using Point = Point<ctype, 3>;

    std::vector<Point> pointIntersections;
    std::vector<TopoDS_Edge> edgeIntersections;
    std::vector<TopoDS_Face> faceIntersections;

    // find possible face intersections
    const auto faceIS = OCCUtilities::intersect(face1, face2, eps);
    const auto faceISFaces = OCCUtilities::getFaces(faceIS);
    for (auto& f : faceISFaces)
        faceIntersections.emplace_back(std::move(f));

    // find possible edge intersections
    const auto faceCut = OCCUtilities::cut(face1, face2, eps);
    const auto faceCutFaces = OCCUtilities::getFaces(faceCut);
    for (unsigned int i = 0; i < faceCutFaces.size(); ++i)
    {
        // edges of the fragments coinciding with the face
        const auto cutFaceWires = OCCUtilities::getWires(faceCutFaces[i]);
        for (const auto& wire : cutFaceWires)
        {
            const auto wireIS = OCCUtilities::intersect(wire, face2, eps);
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

    // Find possible touching points by cutting the wires of the faces
    // with the other faces. We will also find points that are part of previously
    // found edge or face intersections. To this end, we collect all points
    // of the previously found intersections to be able to check for duplicates.
    // TODO: There must be an easier way to find touching points than this here!

    // lambdas to add the vertices of shapes to a container
    auto addVertexShape = [] (const auto&v, auto& c) { c.emplace_back(OCCUtilities::point(v)); };
    auto addVertices = [&] (const auto& shape, auto&c)
    {
        const auto verts = OCCUtilities::getVertices(shape);
        std::for_each(verts.begin(), verts.end(), [&] (const auto& v) { addVertexShape(v, c); });
    };

    // collect all intersection points found so far
    std::vector<Point> curIsPoints;
    for (const auto& s : faceIntersections) addVertices(s, curIsPoints);
    for (const auto& s : edgeIntersections) addVertices(s, curIsPoints);

    // collect all points of the wires of the two faces and add those that result
    // from cutting the wires with the other face (potential touching points)
    const auto wires1 = OCCUtilities::getWires(face1);
    const auto wires2 = OCCUtilities::getWires(face2);

    std::vector<Point> pointsWireCut1, pointsWireCut2;
    for (const auto& s : wires1) addVertices(OCCUtilities::cut(s, face2, eps), pointsWireCut1);
    for (const auto& s : wires2) addVertices(OCCUtilities::cut(s, face1, eps), pointsWireCut2);

    // find those points that live on the other faces' boundary. These will also
    // include the points that are part of face or edge intersections previously found
    using namespace OCCUtilities;
    auto isOnShape = [eps] (const auto& p, const auto& shape)
    { return getVertices(intersect(getShape(p), shape, eps)).size() > 0; };

    std::vector<Point> candidates;
    for (const auto& p : pointsWireCut1)
        if ( isOnShape(p, face2) ) candidates.emplace_back( std::move(p) );
    for (const auto& p : pointsWireCut2)
        if ( isOnShape(p, face1) ) candidates.emplace_back( std::move(p) );

    // lambda to check if a point is in a container
    auto hasPoint = [eps] (const auto& p, const auto& c)
    { return std::any_of(c.begin(), c.end(), [&] (const auto& p2) { return p2.isEqual(p, eps); }); };

    // find those candidates that are not in curIsPoints
    for (auto&& p : candidates)
        if (!hasPoint(p, curIsPoints) && !hasPoint(p, pointIntersections))
            pointIntersections.emplace_back( std::move(p) );

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

#endif // FRACKIT_FACE_FACE_FACE_3D_INTERSECTION_HH
