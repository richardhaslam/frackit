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
 *        between a lateral cylinder surface and
 *        a polygon in 3d space.
 */
#ifndef FRACKIT_CYLINDERSURFACE_POLYGON_INTERSECTION_HH
#define FRACKIT_CYLINDERSURFACE_POLYGON_INTERSECTION_HH

#include <cmath>
#include <stdexcept>

#include <frackit/geometry/polygon.hh>
#include <frackit/geometry/cylindersurface.hh>

#include <frackit/intersection/intersectiontraits.hh>
#include "algo_cylsurface_planargeom.hh"

namespace Frackit {
namespace IntersectionAlgorithms {

//! Intersect a lateral cylinder surface and a polygon
//! The result can be:
//! - an ellipse
//! - ellipse arc(s)
//! - segment(s)
//! - touching points
template<class ctype>
Intersection< CylinderSurface<ctype>, Quadrilateral<ctype, 3> >
intersect_cylinderSurface_polygon(const CylinderSurface<ctype>& cylSurface,
                                  const Polygon<ctype, 3>& polygon,
                                  ctype eps)
{
    if (!polygon.isConvex())
        throw std::runtime_error("CylinderSurface-Polygon algorithm only works for convex polygons");

    using std::max;
    using Segment = typename Polygon<ctype, 3>::Segment;

    ctype charLength = 0.0;
    for (unsigned int cIdx = 0; cIdx < polygon.numCorners(); ++cIdx)
        charLength = max(charLength, 2.0*Segment(polygon.center(), polygon.corner(cIdx)).squaredLength());

    using std::sqrt;
    return intersect_cylinderSurface_planarGeometry(cylSurface, polygon, sqrt(charLength), eps);
}

} // end namespace IntersectionAlgorithms
} // end namespace Frackit

#endif // FRACKIT_CYLINDERSURFACE_POLYGON_INTERSECTION_HH
