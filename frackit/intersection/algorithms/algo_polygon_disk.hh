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
 *        between a polygon and a disk.
 */
#ifndef FRACKIT_POLYGON_DISK_INTERSECTION_HH
#define FRACKIT_POLYGON_DISK_INTERSECTION_HH

#include <cmath>

#include <frackit/precision/precision.hh>
#include <frackit/geometry/polygon.hh>
#include <frackit/geometry/disk.hh>
#include <frackit/geometry/segment.hh>

#include <frackit/intersection/intersectiontraits.hh>
#include "algo_planargeom_planargeom.hh"

namespace Frackit {
namespace IntersectionAlgorithms {

//! Intersect a polygon and a disk
//! The result can be:
//! - a surface bounded by segments and/or elliptical arcs
//! - a segment
//! - a point
//! - no intersection
template<class ctype>
Intersection< Polygon<ctype, 3>, Disk<ctype> >
intersect_polygon_disk(const Polygon<ctype, 3>& polygon,
                             const Disk<ctype>& disk,
                             ctype eps)
{
    using Segment = Segment<ctype, 3>;
    ctype charLength = 0.0;

    using std::max;
    for (unsigned int i = 0; i < polygon.numCorners(); ++i)
        charLength = max(charLength, Segment(polygon.center(), polygon.corner(i)).squaredLength());

    using std::sqrt;
    return intersect_planarGeometry_planarGeometry(polygon, disk,
                                                   sqrt(charLength),
                                                   eps, eps, eps);
}

} // end namespace IntersectionAlgorithms
} // end namespace Frackit

#endif // FRACKIT_POLYGON_DISK_INTERSECTION_HH
