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
 * \ingroup Magnitude
 * \brief Contains functionality for computing the
 *        magnitude (length/area/volume) of geometries.
 */
#ifndef FRACKIT_MAGNITUDE_HH
#define FRACKIT_MAGNITUDE_HH

#include <type_traits>

#include <gp_Pnt.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>

#include <frackit/precision/precision.hh>
#include "length.hh"
#include "area.hh"
#include "volume.hh"

namespace Frackit {

/*!
 * \ingroup Magnitude
 * \brief Returns the magnitude of zero-dimensional geometries (points)
 */
template<class Geometry, std::enable_if_t<Geometry::myDimension() == 0, int> = 0>
typename Geometry::ctype computeMagnitude(const Geometry& geometry)
{ return 0.0; }

/*!
 * \ingroup Magnitude
 * \brief Returns the length of one-dimensional geometries
 */
template<class Geometry, std::enable_if_t<Geometry::myDimension() == 1, int> = 0>
typename Geometry::ctype computeMagnitude(const Geometry& geometry)
{ return computeLength(geometry); }

/*!
 * \ingroup Magnitude
 * \brief Returns the area of two-dimensional geometries
 */
template<class Geometry, std::enable_if_t<Geometry::myDimension() == 2, int> = 0>
typename Geometry::ctype computeMagnitude(const Geometry& geometry)
{ return computeArea(geometry); }

/*!
 * \ingroup Magnitude
 * \brief Returns the volume of three-dimensional geometries
 */
template<class Geometry, std::enable_if_t<Geometry::myDimension() == 3, int> = 0>
typename Geometry::ctype computeMagnitude(const Geometry& geometry)
{ return computeVolume(geometry); }

/*!
 * \ingroup Magnitude
 * \brief Returns the length of a BRep edge
 */
template<class ctype = double>
ctype computeMagnitude(const TopoDS_Edge& edge)
{ return computeLength(edge); }

/*!
 * \ingroup Magnitude
 * \brief Returns the area of a BRep face
 * \param face The face
 * \param eps Tolerance value to be used
 * \param loc A location; defaults to the origin, however,
 *            higher precision is achieved if a point close
 *            to the actual face is chosen.
 */
template<class ctype = double>
ctype computeMagnitude(const TopoDS_Face& face,
                       ctype eps = Precision<ctype>::confusion(),
                       const gp_Pnt& loc = gp_Pnt())
{ return computeArea(face, eps, loc); }

/*!
 * \ingroup Magnitude
 * \brief Returns the volume of a BRep solid
 * \param solid The solid
 * \param eps Tolerance value to be used
 * \param loc A location; defaults to the origin, however,
 *            higher precision is achieved if a point close
 *            to the actual solid is chosen.
 */
template<class ctype = double>
ctype computeMagnitude(const TopoDS_Solid& solid,
                       ctype eps = Precision<ctype>::confusion(),
                       const gp_Pnt& loc = gp_Pnt())
{ return computeVolume(solid, eps, loc); }

} // end namespace Frackit

#endif // FRACKIT_MAGNITUDE_HH
