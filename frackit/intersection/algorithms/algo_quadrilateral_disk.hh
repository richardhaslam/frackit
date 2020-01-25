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
 *        between a quadrilateral and a disk.
 */
#ifndef FRACKIT_QUADRILATERAL_DISK_INTERSECTION_HH
#define FRACKIT_QUADRILATERAL_DISK_INTERSECTION_HH

#include <cmath>

#include <frackit/precision/precision.hh>
#include <frackit/geometry/quadrilateral.hh>
#include <frackit/geometry/disk.hh>

#include <frackit/intersection/intersectiontraits.hh>
#include "algo_planargeom_planargeom.hh"

namespace Frackit {
namespace IntersectionAlgorithms {

//! Intersect two disks
//! The result can be:
//! - a surface bounded by segments and/or elliptical arcs
//! - a segment
//! - a point
//! - no intersection
template<class ctype>
Intersection< Quadrilateral<ctype, 3>, Disk<ctype> >
intersect_quadrilateral_disk(const Quadrilateral<ctype, 3>& quad,
                             const Disk<ctype>& disk,
                             ctype eps)
{
    using std::max;
    ctype charLength = disk.majorAxisLength();
    for (unsigned int edgeIdx = 0; edgeIdx < quad.numEdges(); ++edgeIdx)
        charLength = max(charLength, quad.edge(edgeIdx).length());

    return intersect_planarGeometry_planarGeometry(quad,
                                                   disk,
                                                   charLength,
                                                   eps,
                                                   Precision<ctype>::confusion(),
                                                   eps);
}

} // end namespace IntersectionAlgorithms
} // end namespace Frackit

#endif // FRACKIT_QUADRILATERAL_DISK_INTERSECTION_HH
