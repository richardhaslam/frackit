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
 *        between two disks.
 */
#ifndef FRACKIT_DISK_DISK_INTERSECTION_HH
#define FRACKIT_DISK_DISK_INTERSECTION_HH

#include <cmath>

#include <frackit/precision/precision.hh>
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
Intersection< Disk<ctype>, Disk<ctype> >
intersect_disk_disk(const Disk<ctype>& disk1,
                    const Disk<ctype>& disk2,
                    ctype eps)
{
    using std::max;
    ctype charLength = disk1.majorAxisLength();
    charLength = max(charLength, disk2.majorAxisLength());

    return intersect_planargeometry_planargeometry(disk1,
                                                   disk2,
                                                   charLength,
                                                   Precision<ctype>::confusion(),
                                                   Precision<ctype>::confusion(),
                                                   eps);
}

} // end namespace IntersectionAlgorithms
} // end namespace Frackit

#endif // FRACKIT_DISK_DISK_INTERSECTION_HH
