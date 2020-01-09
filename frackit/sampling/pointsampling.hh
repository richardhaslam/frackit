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
 * \brief \todo TODO Doc me.
 */
#ifndef FRACKIT_POINT_SAMPLING_HH
#define FRACKIT_POINT_SAMPLING_HH

#include <cmath>
#include <random>

#include <frackit/geometry/cylinder.hh>
#include <frackit/geometry/vector.hh>

namespace Frackit {

/*!
 * \brief \todo TODO Doc me.
 */
template<class ctype, int worldDim>
struct PointSamplerTraits;

/*!
 * \brief \todo TODO Doc me.
 */
template<class ctype>
struct PointSamplerTraits<ctype, 3>
{
    static constexpr int worldDim = 3;
    using DistributionX = std::uniform_real_distribution<ctype>;
    using DistributionY = std::uniform_real_distribution<ctype>;
    using DistributionZ = std::uniform_real_distribution<ctype>;
};

/*!
 * \brief \todo TODO Doc me.
 */
template< class Geometry,
          class T = PointSamplerTraits<typename Geometry::ctype,
                                       Geometry::worldDimension()> >
class GeometryPointSampler;

//
// /*!
//  * \brief \todo TODO Doc me.
//  */
// template<class Scalar>
// class GeometryPointSampler< Box<Scalar> >
// {
//
// public:
//     //! export resulting point type
//     using Point = typename Box<Scalar>::Point;
//
//     /*!
//      * \brief \todo TODO Doc me.
//      */
//     GeometryPointSampler(const Box<Scalar>& box)
//     : generator_(std::random_device{}())
//     , p_x_(box.xMin(), box.xMax())
//     , p_y_(box.yMin(), box.yMax())
//     , p_z_(box.zMin(), box.zMax())
//     {}
//
//     /*!
//      * \brief \todo TODO Doc me.
//      */
//     Point operator() ()
//     {
//         return Point( {p_x_(generator_),
//                        p_y_(generator_),
//                        p_z_(generator_)} );
//     }
//
//  private:
//      std::default_random_engine generator_;
//      std::uniform_real_distribution<Scalar> p_x_;
//      std::uniform_real_distribution<Scalar> p_y_;
//      std::uniform_real_distribution<Scalar> p_z_;
// };

/*!
 * \brief \todo TODO Doc me.
 */
template<class ctype, class T>
class GeometryPointSampler< Cylinder<ctype>, T >
{
    using Disk = typename Cylinder<ctype>::Disk;
    using Vector = Frackit::Vector<ctype, 3>;

public:
    //! export resulting point type
    using Point = typename Cylinder<ctype>::Point;

    //! export traits class
    using Traits = T;

    /*!
     * \brief \todo TODO Doc me.
     */
    GeometryPointSampler(const Cylinder<ctype>& cylinder)
    : bottom_(cylinder.bottomFace())
    , generator_(std::random_device{}())
    , p_radius_(0.0, cylinder.radius())
    , p_angle_(0.0, 2*M_PI)
    , p_height_(0.0, cylinder.height())
    {}

    /*!
     * \brief \todo TODO Doc me.
     */
    Point operator() ()
    {
        const auto r = p_radius_(generator_);
        const auto phi = p_angle_(generator_);
        const auto h = p_height_(generator_);

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
     Disk bottom_;

     std::default_random_engine generator_;
     typename Traits::DistributionX p_radius_;
     typename Traits::DistributionY p_angle_;
     typename Traits::DistributionZ p_height_;
};

} // end namespace Frackit

#endif // FRACKIT_POINT_SAMPLING_HH
