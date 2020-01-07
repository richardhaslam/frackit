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
 * \brief Contains functionality for computing the
 *        magnitude (length/area/volume) of geometries.
 */
#ifndef FRACKIT_MAGNITUDE_HH
#define FRACKIT_MAGNITUDE_HH

#include <type_traits>

#include "length.hh"
#include "area.hh"
#include "volume.hh"

namespace Frackit {

/*!
 * \brief Returns the magnitude of zero-dimensional geometries (points)
 */
template<class Geometry, std::enable_if_t<Geometry::myDimension() == 0, int> = 0>
typename Geometry::ctype computeMagnitude(const Geometry& geometry)
{ return 0.0; }

/*!
 * \brief Returns the length of one-dimensional geometries
 */
template<class Geometry, std::enable_if_t<Geometry::myDimension() == 1, int> = 0>
typename Geometry::ctype computeMagnitude(const Geometry& geometry)
{ return computeLength(geometry); }

/*!
 * \brief Returns the area of two-dimensional geometries
 */
template<class Geometry, std::enable_if_t<Geometry::myDimension() == 2, int> = 0>
typename Geometry::ctype computeMagnitude(const Geometry& geometry)
{ return computeArea(geometry); }

/*!
 * \brief Returns the volume of three-dimensional geometries
 */
template<class Geometry, std::enable_if_t<Geometry::myDimension() == 3, int> = 0>
typename Geometry::ctype computeMagnitude(const Geometry& geometry)
{ return computeVolume(geometry); }

} // end namespace Frackit

#endif // FRACKIT_MAGNITUDE_HH
