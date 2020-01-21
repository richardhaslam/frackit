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
 * \brief Class to randomly generate geometries.
 */
#ifndef FRACKIT_GEOMETRY_SAMPLING_HH
#define FRACKIT_GEOMETRY_SAMPLING_HH

namespace Frackit {

//! Forward declaration of the default traits
template<class Geometry> struct DefaultSamplerTraits;

/*!
 * \brief Forward declaration of sampler classes used to
 *        randomly generate geometries. The sampler class
 *        is specialized for different geometries, and the
 *        required probability distribution functions are
 *        injected via corresponding traits classes.
 *        A default traits class should be provided by each
 *        sampler implementation.
 */
template<class Geometry, class T = DefaultSamplerTraits<Geometry>>
struct GeometrySampler;

} // end namespace Frackit

// implementations of the sampler classes
#include "disksampler.hh"

#endif // FRACKIT_GEOMETRY_SAMPLING_HH
