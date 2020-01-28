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
 *        between two quadrilaterals.
 */
#ifndef FRACKIT_QUADRILATERAL_QUADRILATERAL_INTERSECTION_HH
#define FRACKIT_QUADRILATERAL_QUADRILATERAL_INTERSECTION_HH

#include <cmath>

#include <frackit/precision/precision.hh>
#include <frackit/geometry/quadrilateral.hh>

#include <frackit/intersection/intersectiontraits.hh>
#include "algo_planargeom_planargeom.hh"

namespace Frackit {
namespace IntersectionAlgorithms {

//! Intersect two quadrilaterals in 3d space
//! The result can be:
//! - a polygon bounded by segments
//! - a segment
//! - a point
//! - no intersection
template<class ctype>
Intersection< Quadrilateral<ctype, 3>, Quadrilateral<ctype, 3> >
intersect_quadrilateral_quadrilateral(const Quadrilateral<ctype, 3>& quad1,
                                      const Quadrilateral<ctype, 3>& quad2,
                                      ctype eps)
{
    using std::max;
    ctype charLength = 0.0;
    for (unsigned int edgeIdx = 0; edgeIdx < quad1.numEdges(); ++edgeIdx)
        charLength = max(charLength, quad1.edge(edgeIdx).length());
    for (unsigned int edgeIdx = 0; edgeIdx < quad2.numEdges(); ++edgeIdx)
        charLength = max(charLength, quad2.edge(edgeIdx).length());

    return intersect_planarGeometry_planarGeometry(quad1,
                                                   quad2,
                                                   charLength,
                                                   eps,
                                                   eps,
                                                   eps);
}

} // end namespace IntersectionAlgorithms
} // end namespace Frackit

#endif // FRACKIT_QUADRILATERAL_QUADRILATERAL_INTERSECTION_HH
