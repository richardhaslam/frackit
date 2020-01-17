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
 * \brief Class to randomly generate points within geometries.
 *        The type of distribution functions to be used can be
 *        provided via a traits class, where a distribution for
 *        each coordinate of the geometry's local coordinate system
 *        is expected.
 */
#ifndef FRACKIT_POINT_SAMPLING_HH
#define FRACKIT_POINT_SAMPLING_HH

#include <cmath>
#include <random>

#include <frackit/geometry/cylinder.hh>
#include <frackit/geometry/box.hh>
#include <frackit/geometry/vector.hh>

namespace Frackit {

/*!
 * \brief Class to randomly generate points within geometries.
 * \tparam Geometry The geometry on which to sample points
 * \tparam T The traits class containing the type of distributions
 *           to be used.
 * \note The traits class is expected to provide the type of
 *       distribution functions to be used for each coordinate
 *       direction of the geometry's local basis. That is, for
 *       one-dimensional geometries the traits must define the
 *       alias:
 *       using DistroBase1 = ...;
 *       to export the distribution type used along the geometry's
 *       local coordinate. Similarly, for two and three-dimensional
 *       geometries the aliases
 *       using DistroBase2 = ...;
 *       using DistroBase3 = ...;
 *       for the second and/or third dimension of the geometry.
 */
template< class Geometry, class T >
class GeometryPointSampler;

/*!
 * \brief Traits class to be used for sampling on
 *        uniform distributions in all coordinate directions.
 * \tparam ctype The type used for coordinates
 * \tparam worldDim The dimension of the coordinate space.
 */
template<class ctype, int worldDim>
struct UniformPointSamplerTraits;

/*!
 * \brief Traits class to be used for uniform sampling
 *         of points on 1-dimensional geometries.
 * \tparam ctype The type used for coordinates
 */
template<class ctype>
struct UniformPointSamplerTraits<ctype, 1>
{
    using DistroBase1 = std::uniform_real_distribution<ctype>;
};

/*!
 * \brief Traits class to be used for uniform sampling
 *         of points on 2-dimensional geometries.
 * \tparam ctype The type used for coordinates
 */
template<class ctype>
struct UniformPointSamplerTraits<ctype, 2>
{
    using DistroBase1 = std::uniform_real_distribution<ctype>;
    using DistroBase2 = std::uniform_real_distribution<ctype>;
};

/*!
 * \brief Traits class to be used for uniform sampling
 *         of points on 3-dimensional geometries.
 * \tparam ctype The type used for coordinates
 */
template<class ctype>
struct UniformPointSamplerTraits<ctype, 3>
{
    using DistroBase1 = std::uniform_real_distribution<ctype>;
    using DistroBase2 = std::uniform_real_distribution<ctype>;
    using DistroBase3 = std::uniform_real_distribution<ctype>;
};

/*!
 * \brief Convenience function to create a point sampler that
 *        uniformly samples points on a box.
 */
template< class ctype >
GeometryPointSampler< Box<ctype>, UniformPointSamplerTraits<ctype, 3> >
makeUniformPointSampler(const Box<ctype>& box)
{
    using Traits = UniformPointSamplerTraits<ctype, 3>;
    using Sampler = GeometryPointSampler<Box<ctype>, Traits>;
    return Sampler( typename Traits::DistroBase1(box.xMin(), box.xMax()),
                    typename Traits::DistroBase2(box.yMin(), box.yMax()),
                    typename Traits::DistroBase3(box.zMin(), box.zMax()) );
}

/*!
 * \brief Convenience function to create a point sampler that
 *        uniformly samples points on a cylinder.
 */
template< class ctype >
GeometryPointSampler< Cylinder<ctype>, UniformPointSamplerTraits<ctype, 3> >
makeUniformPointSampler(const Cylinder<ctype>& cylinder)
{
    using Traits = UniformPointSamplerTraits<ctype, 3>;
    using Sampler = GeometryPointSampler<Cylinder<ctype>, Traits>;
    return Sampler( cylinder,
                    typename Traits::DistroBase1(0.0, cylinder.radius()),
                    typename Traits::DistroBase2(0.0, 2*M_PI),
                    typename Traits::DistroBase3(0.0, cylinder.height()) );
}

/*!
 * \brief Specialization of the sampler class for boxes.
 */
template<class ctype, class T>
class GeometryPointSampler< Box<ctype>, T >
{
    using DistroX = typename T::DistroBase1;
    using DistroY = typename T::DistroBase2;
    using DistroZ = typename T::DistroBase3;

public:
    //! export resulting point type
    using Point = typename Box<ctype>::Point;

    //! export traits class
    using Traits = T;

    /*!
     * \brief Constructor from distributions.
     */
    GeometryPointSampler(const DistroX& px,
                         const DistroY& py,
                         const DistroZ& pz)
    : generator_(std::random_device{}())
    , p_x_(px)
    , p_y_(py)
    , p_z_(pz)
    {}

    /*!
     * \brief \todo TODO Doc me.
     */
    Point operator() ()
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

/*!
 * \brief \todo TODO Doc me.
 */
template<class ctype, class T>
class GeometryPointSampler< Cylinder<ctype>, T >
{
    using DistroRadius = typename T::DistroBase1;
    using DistroAngle = typename T::DistroBase2;
    using DistroHeight = typename T::DistroBase3;

public:
    //! export resulting point type
    using Point = typename Cylinder<ctype>::Point;

    //! export traits class
    using Traits = T;

    /*!
     * \brief Constructor from a cylinder and distributions.
     */
    GeometryPointSampler(const Cylinder<ctype>& cylinder,
                         const DistroRadius& pr,
                         const DistroAngle& pPhi,
                         const DistroHeight& ph)
    : bottom_(cylinder.bottomFace())
    , generator_(std::random_device{}())
    , p_radius_(pr)
    , p_angle_(pPhi)
    , p_height_(ph)
    {}

    /*!
     * \brief Sample a point from the distributions.
     */
    Point operator() ()
    {
        const auto r = p_radius_(generator_);
        const auto phi = p_angle_(generator_);
        const auto h = p_height_(generator_);

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
     DistroRadius p_radius_;
     DistroAngle  p_angle_;
     DistroHeight p_height_;
};

} // end namespace Frackit

#endif // FRACKIT_POINT_SAMPLING_HH
