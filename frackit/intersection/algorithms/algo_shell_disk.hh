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
#include <frackit/occ/breputilities.hh>

#include <frackit/intersection/intersectiontraits.hh>
#include <frackit/intersection/emptyintersection.hh>
#include "algo_find_touching_points.hh"

namespace Frackit {
namespace IntersectionAlgorithms {

//! Intersect a shell and a disk
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
        throw std::runtime_error("Could not build solid");

    const auto& solid = TopoDS::Solid(makeSolid.Shape());
    const auto& diskFace = OCCUtilities::getShape(disk);
    const auto& containedShape = OCCUtilities::intersect(diskFace, solid, eps);
    const auto containedFaces = OCCUtilities::getFaces(containedShape);
    assert(containedFaces.size() <= 1);

    // only touching points are possible
    if (containedFaces.size() == 0)
    {
        auto isPoints = find_touching_points(diskFace, shell, eps);
        if (!isPoints.empty())
        {
            ResultType result;
            for (auto&& p : isPoints)
                result.emplace_back(std::move(p));
            return result;
        }
        return ResultType({EmptyIntersection<3>()});
    }

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
    {
        auto isPoints = find_touching_points(diskFace, shell, eps);
        if (!isPoints.empty())
        {
            ResultType result;
            for (auto&& p : isPoints)
                result.emplace_back(std::move(p));
            return result;
        }
        return ResultType({EmptyIntersection<3>()});
    }

    // intersection edges
    ResultType result;
    for (auto&& edge : shellIsEdges)
        result.emplace_back( std::move(edge) );
    return result;
}

} // end namespace IntersectionAlgorithms
} // end namespace Frackit

#endif // FRACKIT_SHELL_DISK_INTERSECTION_HH
