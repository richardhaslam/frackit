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

#include <gp_Pnt.hxx>
#include <TopoDS_Face.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>

#include <frackit/geometry/disk.hh>
#include <frackit/geometry/cylindersurface.hh>
#include <frackit/precision/precision.hh>

namespace Frackit {

//! \todo TODO doc me.
template<class Geometry>
typename Geometry::ctype computeArea(const Geometry& geometry)
{ return geometry.area(); }

//! \todo TODO doc me.
template<class ctype>
ctype computeArea(const CylinderSurface<ctype>& cylSurface)
{ return 2.0*M_PI*cylSurface.radius()*cylSurface.height(); }

//! \todo TODO doc me.
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