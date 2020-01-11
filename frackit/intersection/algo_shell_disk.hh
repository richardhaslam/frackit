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
 *        between the shell of a solid and a disk.
 */
#ifndef FRACKIT_SHELL_DISK_INTERSECTION_HH
#define FRACKIT_SHELL_DISK_INTERSECTION_HH

#include <algorithm>
#include <stdexcept>

#include <TopoDS.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>

#include <frackit/geometry/disk.hh>
#include <frackit/precision/precision.hh>
#include <frackit/occ/breputilities.hh>

#include "intersectiontraits.hh"
#include "emptyintersection.hh"

namespace Frackit {
namespace IntersectionAlgorithms {

//! Intersect two planes
//! The result can be:
//! - a point
//! - one or more egdges
//! - a face
template<class ctype>
Intersection< TopoDS_Shell, Disk<ctype> >
intersect_shell_disk(const TopoDS_Shell& shell,
                     const Disk<ctype>& disk,
                     ctype eps)
{
    using ResultType = Intersection< TopoDS_Shell, Disk<ctype> >;

    // first intersect the disk with the solid inside the shell
    BRepBuilderAPI_MakeSolid makeSolid(shell);
    makeSolid.Build();
    if (!makeSolid.IsDone())
        throw std::runtime_error(std::string("Could not build solid"));

    const auto& solid = TopoDS::Solid(makeSolid.Shape());
    const auto& diskFace = OCCUtilities::getShape(disk);
    const auto& containedShape = OCCUtilities::intersect(diskFace, solid, eps);
    const auto containedFaces = OCCUtilities::getFaces(containedShape);
    assert(containedFaces.size() <= 1);

    // lambda to detect if a point is in a vector of points
    auto hasPoint = [&] (const auto& p, const auto& pointVec)
    { return std::any_of(pointVec.begin(),
                         pointVec.end(),
                         [&] (const auto& pw) { return p.isEqual(pw, eps); }); };

    // lambda to find touching points of the disk with the shell
    auto getTouchingPoints = [&] () -> ResultType
    {
        const auto diskWires = OCCUtilities::getWires(diskFace);

        assert(diskWires.size() == 1);
        const auto wireVertices = OCCUtilities::getVertices(diskWires[0]);
        const auto wireEdges = OCCUtilities::getEdges(diskWires[0]);

        const auto wireCut = OCCUtilities::cut(diskWires[0], shell, eps);
        const auto wireCutEdges = OCCUtilities::getEdges(wireCut);

        // there are touching points
        if (wireCutEdges.size() > wireEdges.size())
        {
            std::vector<Point<ctype, 3>> isPoints;
            std::vector<Point<ctype, 3>> wirePoints;
            for (const auto& v : wireVertices)
            {
                // if this point is contained in the shell,
                // directly add it to the intersection points
                const auto pointIS = OCCUtilities::intersect(v, shell, eps);
                const auto pointISPoints = OCCUtilities::getVertices(pointIS);

                if (pointISPoints.empty())
                    wirePoints.push_back(OCCUtilities::point(v));
                else
                {
                    assert(pointISPoints.size() == 1);
                    const auto p = OCCUtilities::point(v);
                    if (!hasPoint(p, isPoints))
                        isPoints.push_back(p);
                }
            }

            // find newly created poitns
            const auto wireCutVertices = OCCUtilities::getVertices(wireCut);
            for (const auto& v : wireCutVertices)
            {
                const auto p = OCCUtilities::point(v);
                if (!hasPoint(p, wirePoints))
                    if (!hasPoint(p, isPoints))
                        isPoints.push_back(p);
            }

            ResultType result;
            for (auto& p : isPoints)
                result.emplace_back(std::move(p));
            return result;
        }

        // no intersection
        return ResultType({ EmptyIntersection<3>() });
    };

    // only touching points are possible
    if (containedFaces.size() == 0)
        return getTouchingPoints();

    // detect intersection face
    const auto& containedFace = containedFaces[0];
    const auto faceIntersection = OCCUtilities::intersect(containedFace, shell, eps);
    const auto faceIsFaces = OCCUtilities::getFaces(faceIntersection);
    if (faceIsFaces.size() > 0)
    {
        assert(faceIsFaces.size() == 1);
        return ResultType({ faceIsFaces[0] });
    }

    // there might be intersection edges
    const auto containedWires = OCCUtilities::getWires(containedFaces[0]);
    assert(containedWires.size() == 1);

    // intersect the wire with the shell
    const auto& containedWire = containedWires[0];
    const auto shellIntersection = OCCUtilities::intersect(containedWire, shell, eps);
    const auto shellIsEdges = OCCUtilities::getEdges(shellIntersection);

    // no intersection -> there might still be touching points
    if (shellIsEdges.empty())
        return getTouchingPoints();

    // intersection edges
    ResultType result;
    for (auto&& edge : shellIsEdges)
        result.emplace_back( std::move(edge) );
    return result;
}

} // end namespace IntersectionAlgorithms
} // end namespace Frackit

#endif // FRACKIT_SHELL_DISK_INTERSECTION_HH
