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
 * \brief Contains the algorithm to find touching points
 *        between a face and a shape.
 */
#ifndef FRACKIT_FIND_TOUCHING_POINTS_HH
#define FRACKIT_FIND_TOUCHING_POINTS_HH

#include <algorithm>
#include <stdexcept>

#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>

#include <frackit/geometry/point.hh>
#include <frackit/occ/breputilities.hh>

namespace Frackit {
namespace IntersectionAlgorithms {

//! Find touching points between a face and a shape.
//! This is to be used in other intersection algorithms
//! after verifying that the result can only be touching
//! points. This must ONLY be used when it is certain that
//! the result will be touching points!
template<class ctype>
std::vector<Point<ctype, 3>>
find_touching_points(const TopoDS_Face& face,
                     const TopoDS_Shape& shape,
                     ctype eps)
{
    std::vector<Point<ctype, 3>> result;

    // lambda to detect if a point is in a vector of points
    auto hasPoint = [&] (const auto& p, const auto& pointVec)
    { return std::any_of(pointVec.begin(),
                         pointVec.end(),
                         [&] (const auto& pw) { return p.isEqual(pw, eps); }); };

    // cut the wire by the shape and check if new edges were created
    const auto faceWires = OCCUtilities::getWires(face);
    if (faceWires.size() != 1)
        throw std::runtime_error(std::string("Algorithm expects face bounded by a single wire"));

    const auto wireCut = OCCUtilities::cut(faceWires[0], shape, eps);

    // detect all points on the wire cut that are on the shape
    const auto wireCutVertices = OCCUtilities::getVertices(wireCut);
    for (const auto& v : wireCutVertices)
    {
        const auto pointIS = OCCUtilities::intersect(v, shape, eps);
        const auto pointISPoints = OCCUtilities::getVertices(pointIS);

        assert(pointISPoints.size() <= 1);
        if (!pointISPoints.empty())
        {
            auto p = OCCUtilities::point(v);
            if (!hasPoint(p, result))
                result.emplace_back(std::move(p));
        }
    }

    return result;
}

} // end namespace IntersectionAlgorithms
} // end namespace Frackit

#endif // FRACKIT_FIND_TOUCHING_POINTS_HH
