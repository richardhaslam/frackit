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
 * \ingroup Magnitude
 * \brief Contains functionality for computing the
 *        areas of two-dimensional geometries.
 */
#ifndef FRACKIT_MAGNITUDE_AREA_HH
#define FRACKIT_MAGNITUDE_AREA_HH

#include <cmath>

#include <gp_Pnt.hxx>
#include <TopoDS_Face.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>

#include <frackit/precision/precision.hh>

namespace Frackit {

/*!
 * \ingroup Magnitude
 * \brief Returns the area of an internal geometry
 *        that has an area() function available.
 */
template<class Geometry>
typename Geometry::ctype computeArea(const Geometry& geometry)
{ return geometry.area(); }

/*!
 * \ingroup Magnitude
 * \brief Returns the area of a TopoDS_Face.
 * \param face The face
 * \param eps Tolerance value to be used
 * \param loc A location; defaults to the origin, however,
 *            higher precision is achieved if a point close
 *            to the actual face is chosen.
 */
template<class ctype = double>
ctype computeArea(const TopoDS_Face& face,
                  ctype eps = Precision<ctype>::confusion(),
                  const gp_Pnt& loc = gp_Pnt())
{
    GProp_GProps gprops(loc);
    BRepGProp::SurfaceProperties(face, gprops, eps);
    return gprops.Mass();
}

} // end namespace Frackit

#endif // FRACKIT_MAGNITUDE_AREA_HH
