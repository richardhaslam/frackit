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
 * \brief Contains functionality for computing the distance
 *        of a geometry to the boundary of another geometry.
 */
#ifndef FRACKIT_DISTANCE_TO_BOUNDARY_HH
#define FRACKIT_DISTANCE_TO_BOUNDARY_HH

#include <Extrema_ExtAlgo.hxx>
#include <Extrema_ExtFlag.hxx>

#include <frackit/precision/precision.hh>
#include <frackit/geometry/disk.hh>

#include "distance.hh"

namespace Frackit {

/*!
 * \brief Compute the distance of a geometry
 *        to the boundary of another geometry.
 * \param geo1 The first geometry
 * \param geo2 The second geometry
 * \param deflection The epsilon used in the BrepExtrema command
 * \param extFlag The flag passed to the BrepExtrema command (MIN/MAX/MINMAX)
 * \param extAlgo The algorithm passed to the BrepExtrema command (TREE/GRAD)
 */
template<class Geo1, class Geo2>
Impl::PCT<Geo1, Geo2>
computeDistanceToBoundary(const Geo1& geo1,
                          const Geo2& geo2,
                          Impl::PCT<Geo1, Geo2> deflection = Precision<Impl::PCT<Geo1, Geo2>>::confusion(),
                          Extrema_ExtFlag extFlag = Extrema_ExtFlag_MINMAX,
                          Extrema_ExtAlgo extAlgo = Extrema_ExtAlgo_Grad)
{
    std::string msg = "Distance to boundary computation not implemented for ";
    msg += "\"" + Geo1::name() + "\"";
    msg += " and ";
    msg += "\"" + Geo2::name() + "\"";
}

/*!
 * \brief Compute the distance of a geometry
 *        to the bounding ellipse of a disk.
 * \param geo The geometry
 * \param disk The disk
 * \param deflection The epsilon used in the BrepExtrema command
 * \param extFlag The flag passed to the BrepExtrema command (MIN/MAX/MINMAX)
 * \param extAlgo The algorithm passed to the BrepExtrema command (TREE/GRAD)
 */
template<class Geo, class ctype>
Impl::PCT<Geo, Disk<ctype>>
computeDistanceToBoundary(const Geo& geo,
                          const Disk<ctype>& disk,
                          Impl::PCT<Geo, Disk<ctype>> deflection = Precision<Impl::PCT<Geo, Disk<ctype>>>::confusion(),
                          Extrema_ExtFlag extFlag = Extrema_ExtFlag_MINMAX,
                          Extrema_ExtAlgo extAlgo = Extrema_ExtAlgo_Grad)
{ return computeDistance(OCCUtilities::getShape(geo), OCCUtilities::getShape(disk.boundingEllipse())); }

} // end namespace Frackit

#endif // FRACKIT_DISTANCE_TO_BOUNDARY_HH
