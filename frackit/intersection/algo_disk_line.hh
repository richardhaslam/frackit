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
 *        between a disk and a line.
 */
#ifndef FRACKIT_DISK_LINE_INTERSECTION_HH
#define FRACKIT_DISK_LINE_INTERSECTION_HH

#include <variant>
#include <stdexcept>

#include <frackit/geometry/disk.hh>
#include <frackit/geometry/line.hh>
#include <frackit/precision/precision.hh>

#include "intersectiontraits.hh"
#include "algo_planargeom_line.hh"

namespace Frackit {
namespace IntersectionAlgorithms {

//! Intersect a disk and a line
//! The result can be:
//! - a segment
//! - a point
//! - no intersection
template<class ctype>
Intersection< Disk<ctype>, Line<ctype, 3> >
intersect_disk_line(const Disk<ctype>& disk,
                    const Line<ctype, 3>& line,
                    ctype eps)
{
    // characteristic length
    const auto charLength = disk.majorAxisLength();
    return intersect_planarGeometry_line(disk, line, charLength, Precision<ctype>::confusion(), eps);
}

} // end namespace IntersectionAlgorithms
} // end namespace Frackit

#endif // FRACKIT_INTERSECT_HH
