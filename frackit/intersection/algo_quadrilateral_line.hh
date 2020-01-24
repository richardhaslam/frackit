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
 *        between a quadrilateral and a line in 3d space.
 */
#ifndef FRACKIT_QUADRILATERAL_LINE_INTERSECTION_HH
#define FRACKIT_QUADRILATERAL_LINE_INTERSECTION_HH

#include <cmath>
#include <limits>

#include <frackit/geometry/quadrilateral.hh>
#include <frackit/geometry/line.hh>

#include "intersectiontraits.hh"
#include "algo_planargeom_line.hh"

namespace Frackit {
namespace IntersectionAlgorithms {

//! Intersect a 3d quadrilateral and a line
//! The result can be:
//! - a segment
//! - a point
//! - no intersection
template<class ctype>
Intersection< Quadrilateral<ctype, 3>, Line<ctype, 3> >
intersect_quadrilateral_line(const Quadrilateral<ctype, 3>& quad,
                             const Line<ctype, 3>& line,
                             ctype eps)
{
    // characteristic length
    ctype charLength = std::numeric_limits<ctype>::max();

    using std::min;
    for (unsigned int edgeIdx = 0; edgeIdx < quad.numEdges(); ++edgeIdx)
        charLength = min(charLength, quad.edge(edgeIdx).length());

    return intersect_planargeometry_line(quad, line, charLength, eps, eps);
}

} // end namespace IntersectionAlgorithms
} // end namespace Frackit

#endif // FRACKIT_QUADRILATERAL_LINE_INTERSECTION_HH
