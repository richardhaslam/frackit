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
 * \ingroup GeometryUtilities
 * \brief Contains a free function to get the bounding box of a geometry.
 */
#ifndef FRACKIT_GEOMETRY_GET_BOUNDING_BOX_HH
#define FRACKIT_GEOMETRY_GET_BOUNDING_BOX_HH

#include <frackit/occ/breputilities.hh>

namespace Frackit {

/*!
 * \ingroup GeometryUtilities
 * \brief Free function to get the bounding box of a geometry.
 *        This is the default implementation that first casts the
 *        geometry into its BRep and then uses the OpenCascade functionalities
 *        to compute the bounding box of it.
 * \note More efficient overloads for specific geometries can be provided
 *       in the headers.
 */
template<class Geo>
auto getBoundingBox(const Geo& geo)
{ return OCCUtilities::getBoundingBox( OCCUtilities::getShape(geo) ); }

} // end namespace Frackit

#endif // FRACKIT_GEOMETRY_GET_BOUNDING_BOX_HH
