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
 *        area of two-dimensional geometries.
 */
#ifndef FRACKIT_MAGNITUDE_AREA_HH
#define FRACKIT_MAGNITUDE_AREA_HH

#include <cmath>

#include <frackit/geometry/disk.hh>
#include <frackit/geometry/cylindersurface.hh>

namespace Frackit {

//! \todo TODO doc me.
template<class ctype>
ctype computeArea(const Disk<ctype>& disk)
{ return M_PI*disk.majorAxisLength()*disk.minorAxisLength(); }

//! \todo TODO doc me.
template<class ctype>
ctype computeArea(const CylinderSurface<ctype>& cylSurface)
{ return 2.0*M_PI*cylSurface.radius()*cylSurface.height(); }

} // end namespace Frackit

#endif // FRACKIT_MAGNITUDE_AREA_HH
