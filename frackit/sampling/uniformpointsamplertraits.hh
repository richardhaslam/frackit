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
 * \ingroup Sampling
 * \brief Traits classes to be used for uniform point samplers.
 */
#ifndef FRACKIT_UNIFORM_POINT_SAMPLER_TRAITS_HH
#define FRACKIT_UNIFORM_POINT_SAMPLER_TRAITS_HH

#include <random>

namespace Frackit {

/*!
 * \ingroup Sampling
 * \brief Traits class to be used for sampling on
 *        uniform distributions in all coordinate directions.
 * \tparam ctype The type used for coordinates
 * \tparam worldDim The dimension of the coordinate space.
 */
template<class ctype, int worldDim>
struct UniformPointSamplerTraits;

/*!
 * \ingroup Sampling
 * \brief Traits class to be used for uniform sampling
 *         of points on 1-dimensional geometries.
 * \tparam ctype The type used for coordinates
 */
template<class ctype>
struct UniformPointSamplerTraits<ctype, 1>
{
    using DistributionBase1 = std::uniform_real_distribution<ctype>;
};

/*!
 * \ingroup Sampling
 * \brief Traits class to be used for uniform sampling
 *         of points on 2-dimensional geometries.
 * \tparam ctype The type used for coordinates
 */
template<class ctype>
struct UniformPointSamplerTraits<ctype, 2>
{
    using DistributionBase1 = std::uniform_real_distribution<ctype>;
    using DistributionBase2 = std::uniform_real_distribution<ctype>;
};

/*!
 * \ingroup Sampling
 * \brief Traits class to be used for uniform sampling
 *         of points on 3-dimensional geometries.
 * \tparam ctype The type used for coordinates
 */
template<class ctype>
struct UniformPointSamplerTraits<ctype, 3>
{
    using DistributionBase1 = std::uniform_real_distribution<ctype>;
    using DistributionBase2 = std::uniform_real_distribution<ctype>;
    using DistributionBase3 = std::uniform_real_distribution<ctype>;
};

} // end namespace Frackit

#endif // FRACKIT_UNIFORM_POINT_SAMPLER_TRAITS_HH
