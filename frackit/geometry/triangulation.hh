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
 * \ingroup Geometry
 * \brief Functionality to compute triangulations of geometric objects.
 */
#ifndef FRACKIT_GEOMETRY_TRIANGULATION_HH
#define FRACKIT_GEOMETRY_TRIANGULATION_HH

#include <cassert>

#include "point.hh"
#include "triangle.hh"

namespace Frackit {

/*!
 * \ingroup Geometry
 * \brief Compute the triangulation of a convex hull,
 *        given as a vector of points.
 * \param corners The points of the convex hull
 * \note this expects the provided container is sorted, that is,
 *       corners[i] and corners[i+1] are actually neighboring
 *       corners of the convex hull.
 */
template<class CT, int worldDim>
std::vector<Triangle<CT, worldDim>>
triangulate(const std::vector<Point<CT, worldDim>>& corners)
{
    const auto numCorners = corners.size();
    assert(numCorners >= 3);

    using Triangle = Frackit::Triangle<CT, worldDim>;
    std::vector<Triangle> result;
    result.reserve(numCorners-2);
    result.emplace_back(corners[0], corners[1], corners[2]);

    for (std::size_t i = 2; i < numCorners-1; ++i)
        result.emplace_back(corners[0], corners[i], corners[i+1]);

    assert(result.size() == numCorners-2);
    return result;
}

} // end namespace Frackit

#endif // FRACKIT_GEOMETRY_TRIANGULATION_HH
