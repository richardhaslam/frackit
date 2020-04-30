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
 * \brief Contains the intersection algorithm between the mantle of
 *        a cylinder surface and a disk.
 */
#ifndef FRACKIT_CYLINDERMANTLE_DISK_INTERSECTION_HH
#define FRACKIT_CYLINDERMANTLE_DISK_INTERSECTION_HH

#include <frackit/geometry/disk.hh>
#include <frackit/geometry/cylindermantle.hh>

#include <frackit/intersection/intersectiontraits.hh>
#include "algo_cylmantle_planargeom.hh"

namespace Frackit {
namespace IntersectionAlgorithms {

//! Intersect a lateral cylinder surface and a disk
//! The result can be:
//! - an ellipse
//! - ellipse arc(s)
//! - segment(s)
//! - touching points
template<class ctype>
Intersection< CylinderMantle<ctype>, Disk<ctype> >
intersect_cylinderMantle_disk(const CylinderMantle<ctype>& cylMantle,
                               const Disk<ctype>& disk,
                               ctype eps)
{ return intersect_cylinderMantle_planarGeometry(cylMantle, disk, disk.majorAxisLength(), eps); }

} // end namespace IntersectionAlgorithms
} // end namespace Frackit

#endif // FRACKIT_CYLINDERMANTLE_DISK_INTERSECTION_HH
