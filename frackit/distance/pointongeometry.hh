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
 * \ingroup Distance
 * \brief Contains functionality to evaluate if points are contained
 *        on geometries, i.e. if the distance is below a numerical threshold.
 */
#ifndef FRACKIT_POINT_ON_GEOMETRY_HH
#define FRACKIT_POINT_ON_GEOMETRY_HH

#include <TopoDS_Shape.hxx>

#include <frackit/common/extractdimension.hh>
#include <frackit/precision/defaultepsilon.hh>
#include <frackit/geometry/point.hh>

#include "distance.hh"

namespace Frackit {

/*!
 * \ingroup Distance
 * \brief Evaluate if a point lies on a geometry.
 * \param p The point
 * \param geo The geometry
 * \param eps Epsilon value to be used for the check
 * \note In the general case we forward to the internal function of the geometry
 */
template<class ctype1, int wd, class Geo, class ctype2>
bool pointOnGeometry(const Point<ctype1, wd>& p, const Geo& geo, ctype2 eps)
{
    static_assert(wd == DimensionalityTraits<Geo>::geomDim, "World dimension mismatch");
    return geo.contains(p, eps);
}

/*!
 * \ingroup Distance
 * \brief Evaluate if a point lies on a shape.
 * \param p The point
 * \param shape The shape
 * \param eps Epsilon value to be used for the check
 */
template<class ctype1, class Shape, class ctype2>
bool pointOnGeometry(const Point<ctype1, 3>& p, const Shape& shape, ctype2 eps)
{
    const auto distance = computeDistance(p, shape);
    return distance < eps;
}

/*!
 * \ingroup Distance
 * \brief Evaluate if a point lies on a geometry.
 * \param p The point
 * \param geo The geometry
 * \note This overload selects the default epsilon
 */
template<class ctype1, int wd, class Geo>
bool pointOnGeometry(const Point<ctype1, wd>& p, const Geo& geo)
{ return pointOnGeometry(p, geo, defaultEpsilon(geo)); }

} // end namespace Frackit

#endif // FRACKIT_POINT_ON_GEOMETRY_HH
