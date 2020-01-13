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
 *        between a face shape and a disk.
 */
#ifndef FRACKIT_FACE_DISK_INTERSECTION_HH
#define FRACKIT_FACE_DISK_INTERSECTION_HH

#include <algorithm>
#include <stdexcept>

#include <TopoDS_Vertex.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>

#include <frackit/geometry/point.hh>
#include <frackit/geometry/disk.hh>
#include <frackit/occ/breputilities.hh>

#include "algo_find_touching_points.hh"
#include "intersectiontraits.hh"
#include "emptyintersection.hh"

namespace Frackit {
namespace IntersectionAlgorithms {

//! Intersect a face shape and a disk
//! The result can be composed of:
//! - points
//! - edges
//! - faces
//! Multiple of the above are possible
//! since the face shape might be curved.
template<class ctype>
Intersection< TopoDS_Face, Disk<ctype> >
intersect_face_disk(const TopoDS_Face& face,
                    const Disk<ctype>& disk,
                    ctype eps)
{
    using ResultType = Intersection< TopoDS_Face, Disk<ctype> >;

    std::vector<Point<ctype, 3>> pointIntersections;
    std::vector<TopoDS_Edge> edgeIntersections;
    std::vector<TopoDS_Face> faceIntersections;

    // find possible face intersections
    const auto diskShape = OCCUtilities::getShape(disk);
    const auto diskIS = OCCUtilities::intersect(diskShape, face, eps);
    const auto diskISFaces = OCCUtilities::getFaces(diskIS);
    for (auto& f : diskISFaces)
        faceIntersections.emplace_back(std::move(f));

    // find possible edge intersections
    const auto diskCut = OCCUtilities::cut(diskShape, face, eps);
    const auto diskCutFaces = OCCUtilities::getFaces(diskCut);
    for (unsigned int i = 0; i < diskCutFaces.size(); ++i)
    {
        // edges of the fragments coinciding with the face
        const auto cutFaceWires = OCCUtilities::getWires(diskCutFaces[i]);
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

    const auto touchPoints = find_touching_points(diskShape, face, eps);
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

#endif // FRACKIT_FACE_DISK_INTERSECTION_HH
