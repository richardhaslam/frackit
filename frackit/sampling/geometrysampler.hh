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
 * \brief Class that defines the interface for
 *        sampler classes of geometries.
 */
#ifndef FRACKIT_GEOMETRY_SAMPLER_HH
#define FRACKIT_GEOMETRY_SAMPLER_HH

namespace Frackit {

/*!
 * \ingroup Sampling
 * \brief Interface for geometry sampler classes.
 *        Sampler class implemetations must implement
 *        the () operator, with which a random geometry
 *        of the given type is created.
 * \tparam G The type of the geometry to be sampled.
 */
template<class G>
class GeometrySampler
{

public:
    //! Export the geometry type that is sampled
    using Geometry = G;

    //! every abstract base class has a virtual destructor
    virtual ~GeometrySampler () {}

    /*!
     * \brief Creates an object of Geometry.
     */
    virtual Geometry operator() () = 0;
};

} // end namespace Frackit

#endif // FRACKIT_GEOMETRY_SAMPLER_HH
