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
 *        between two polygons in 3d space.
 */
#ifndef FRACKIT_POLYGON_POLYGON_INTERSECTION_HH
#define FRACKIT_POLYGON_POLYGON_INTERSECTION_HH

#include <cmath>
#include <stdexcept>

#include <frackit/geometry/polygon.hh>
#include <frackit/intersection/intersectiontraits.hh>
#include "algo_planargeom_planargeom.hh"

namespace Frackit {
namespace IntersectionAlgorithms {

//! Intersect two polygons in 3d space
//! The result can be:
//! - a polygon bounded by segments
//! - a segment
//! - a point
//! - no intersection
template<class ctype>
Intersection< Polygon<ctype, 3>, Polygon<ctype, 3> >
intersect_polygon_polygon(const Polygon<ctype, 3>& polygon1,
                          const Polygon<ctype, 3>& polygon2,
                          ctype eps)
{
    if (!polygon1.isConvex() || !polygon2.isConvex())
        throw std::runtime_error("Polygon-Polygon algorithm only works for convex polygons");

    using std::max;
    using std::sqrt;
    using Segment = typename Polygon<ctype, 3>::Segment;

    ctype charLength = 0.0;
    for (unsigned int cIdx = 0; cIdx < polygon1.numCorners(); ++cIdx)
        charLength = max(charLength, 2.0*Segment(polygon1.center(), polygon1.corner(cIdx)).squaredLength());
    for (unsigned int cIdx = 0; cIdx < polygon2.numCorners(); ++cIdx)
        charLength = max(charLength, 2.0*Segment(polygon2.center(), polygon2.corner(cIdx)).squaredLength());

    return intersect_planarGeometry_planarGeometry(polygon1, polygon2,
                                                   sqrt(charLength),
                                                   eps, eps, eps);
}

} // end namespace IntersectionAlgorithms
} // end namespace Frackit

#endif // FRACKIT_POLYGON_POLYGON_INTERSECTION_HH
