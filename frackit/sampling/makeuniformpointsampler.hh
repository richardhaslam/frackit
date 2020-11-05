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
 * \brief Free functions to create point samplers from geometries with
 *        uniform distributions used for all coordinate directions.
 */
#ifndef FRACKIT_MAKE_UNIFORM_POINT_SAMPLER_HH
#define FRACKIT_MAKE_UNIFORM_POINT_SAMPLER_HH

#include <frackit/geometry/ctype.hh>
#include <frackit/common/extractdimension.hh>

#include "uniformpointsamplertraits.hh"
#include "makepointsampler.hh"

namespace Frackit {

/*!
 * \ingroup Sampling
 * \brief Free function to create point samplers,
 *        for the provided geometry, with uniform
 *        distributions for all coordinate directions.
 */
template< class Geometry, class... Args >
auto makeUniformPointSampler(const Geometry& geometry, Args&&... args)
{
    using ctype = typename CoordinateTypeTraits<Geometry>::type;
    static constexpr auto dim = DimensionalityTraits<Geometry>::geometryDimension();

    using Traits = UniformPointSamplerTraits<ctype, dim>;
    return makePointSampler<Traits>(geometry, std::forward<Args...>(args)...);
}

} // end namespace Frackit

#endif // FRACKIT_MAKE_UNIFORM_POINT_SAMPLER_HH
