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
 * \brief Class to randomly generate points on boxes.
 */
#ifndef FRACKIT_BOX_POINT_SAMPLER_HH
#define FRACKIT_BOX_POINT_SAMPLER_HH

#include <random>
#include <frackit/geometry/box.hh>
#include "pointsampler.hh"

namespace Frackit {

/*!
 * \brief Class to randomly generate points within boxes
 *        following the provided distributions. The type of
 *        distribution to be used can be specified via traits,
 *        which are required to exporrt the distributions to be
 *        used for the x, y and z-coordinate, using the naming
 *        convention for three-dimensional geometries,
 *        i.e. they must export:
 *        using DistributionBase1 = ...; // Distribution for x-coordinate
 *        using DistributionBase2 = ...; // Distribution for y-coordinate
 *        using DistributionBase3 = ...; // Distribution for z-coordinate
 *
 * \tparam ctype The type used for coordinates
 * \tparam T The traits class containing the distributions
 */
template<class ctype, class T>
class BoxPointSampler
: public PointSampler< typename Box<ctype>::Point >
{
    using DistroX = typename T::DistributionBase1;
    using DistroY = typename T::DistributionBase2;
    using DistroZ = typename T::DistributionBase3;

public:
    //! export resulting point type
    using Point = typename Box<ctype>::Point;

    //! export traits class
    using Traits = T;

    /*!
     * \brief Constructor from distributions.
     */
    BoxPointSampler(const DistroX& px,
                    const DistroY& py,
                    const DistroZ& pz)
    : generator_(std::random_device{}())
    , p_x_(px)
    , p_y_(py)
    , p_z_(pz)
    {}

    /*!
     * \brief Sample a point from the distributions.
     */
    Point operator() () override
    {
        return Point( {p_x_(generator_),
                       p_y_(generator_),
                       p_z_(generator_)} );
    }

 private:
    std::default_random_engine generator_;
    DistroX p_x_;
    DistroY p_y_;
    DistroZ p_z_;
};

} // end namespace Frackit

#endif // FRACKIT_BOX_POINT_SAMPLER_HH
