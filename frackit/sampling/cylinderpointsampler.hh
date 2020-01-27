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
 * \brief Class to randomly generate points on cylinders.
 */
#ifndef FRACKIT_CYLINDER_POINT_SAMPLER_HH
#define FRACKIT_CYLINDER_POINT_SAMPLER_HH

#include <cmath>
#include <random>

#include <frackit/geometry/cylinder.hh>
#include <frackit/geometry/vector.hh>

#include "pointsampler.hh"

namespace Frackit {

/*!
 * \ingroup Sampling
 * \brief Class to randomly generate points within cylinders
 *        following the provided distributions. The type of
 *        distribution to be used can be specified via traits,
 *        which are required to exporrt the distributions to be
 *        used for the squared radius, the angle and the height,
 *        using the naming convention for three-dimensional geometries,
 *        i.e. they must export:
 *        using DistributionBase1 = ...; // Distribution for squared radius
 *        using DistributionBase2 = ...; // Distribution for the angle
 *        using DistributionBase3 = ...; // Distribution for the height
 *
 * \tparam ctype The type used for coordinates
 * \tparam T The traits class containing the distributions
 */
template<class ctype, class T>
class CylinderPointSampler
: public PointSampler< typename Cylinder<ctype>::Point >
{
    using DistroRadiusSquared = typename T::DistributionBase1;
    using DistroAngle = typename T::DistributionBase2;
    using DistroHeight = typename T::DistributionBase3;

public:
    //! export resulting point type
    using Point = typename Cylinder<ctype>::Point;

    //! export traits class
    using Traits = T;

    /*!
     * \brief Constructor from a cylinder and distributions.
     */
    CylinderPointSampler(const Cylinder<ctype>& cylinder,
                         const DistroRadiusSquared& pr,
                         const DistroAngle& pPhi,
                         const DistroHeight& ph)
    : bottom_(cylinder.bottomFace())
    , generator_(std::random_device{}())
    , p_radiusSquared_(pr)
    , p_angle_(pPhi)
    , p_height_(ph)
    {}

    /*!
     * \brief Sample a point from the distributions.
     * \note The distribution for the radius is expected
     *       to actually be a distribution for the squared
     *       radius, i.e. it is in the interval [0, r^2].
     *       After sampling, the square root is taken.
     */
    Point operator() () override
    {
        const auto phi = p_angle_(generator_);
        const auto h = p_height_(generator_);
        auto r = p_radiusSquared_(generator_);

        using std::sqrt;
        r = sqrt(r);

        using Vector = Frackit::Vector<ctype, 3>;
        auto a = Vector(bottom_.majorAxis());
        auto b = Vector(bottom_.minorAxis());
        auto n = Vector(bottom_.normal());

        using std::sin;
        using std::cos;
        a *= r*cos(phi);
        b *= r*sin(phi);
        n *= h;

        Point result = bottom_.center();
        result += a;
        result += b;
        result += n;

        return result;
    }

 private:
     typename Cylinder<ctype>::Disk bottom_;

     std::default_random_engine generator_;
     DistroRadiusSquared p_radiusSquared_;
     DistroAngle  p_angle_;
     DistroHeight p_height_;
};

} // end namespace Frackit

#endif // FRACKIT_CYLINDER_POINT_SAMPLER_HH
