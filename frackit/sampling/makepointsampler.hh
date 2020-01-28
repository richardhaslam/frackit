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
 * \ingroup Sampling
 * \brief Free functions to create point samplers from geometries.
 */
#ifndef FRACKIT_MAKE_POINT_SAMPLER_HH
#define FRACKIT_MAKE_POINT_SAMPLER_HH

#include "cylinderpointsampler.hh"
#include "boxpointsampler.hh"

namespace Frackit {

/*!
 * \ingroup Sampling
 * \brief Overload of the free function to create point
 *        samplers for cylinders.
 * \note Sampling along the radius coordinate occurs in the interval [0, r^2]
 *       and the square root is taken after sampling. Thus, the provided distribution
 *       must be for the squared radius rather than the radius itself.
 */
template< class Traits, class ctype >
CylinderPointSampler< ctype, Traits >
makePointSampler(const Cylinder<ctype>& cylinder)
{
    using Sampler = CylinderPointSampler<ctype, Traits>;
    return Sampler( cylinder,
                    typename Traits::DistributionBase1(0.0, cylinder.radius()*cylinder.radius()),
                    typename Traits::DistributionBase2(0.0, 2*M_PI),
                    typename Traits::DistributionBase3(0.0, cylinder.height()) );
}

/*!
 * \ingroup Sampling
 * \brief Overload of the free function to create point
 *        samplers for boxes.
 */
template< class Traits, class ctype >
BoxPointSampler< ctype, Traits >
makePointSampler(const Box<ctype>& box)
{
    using Sampler = BoxPointSampler<ctype, Traits>;
    return Sampler( typename Traits::DistributionBase1(box.xMin(), box.xMax()),
                    typename Traits::DistributionBase2(box.yMin(), box.yMax()),
                    typename Traits::DistributionBase3(box.zMin(), box.zMax()) );
}

} // end namespace Frackit

#endif // FRACKIT_MAKE_POINT_SAMPLER_HH
