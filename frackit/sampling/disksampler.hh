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
#ifndef FRACKIT_DISK_SAMPLER_HH
#define FRACKIT_DISK_SAMPLER_HH

#include <random>

#include <frackit/common/math.hh>
#include <frackit/geometry/disk.hh>
#include <frackit/geometry/vector.hh>
#include <frackit/geometry/direction.hh>

#include "geometrysampling.hh"

namespace Frackit {

/*!
 * \brief Specialization of the default traits for the disk sampler.
 */
template<class ctype>
struct DefaultSamplerTraits<Disk<ctype>>
{
    using MajorAxisLengthDistribution = std::normal_distribution<ctype>;
    using MinorAxisLengthDistribution = std::normal_distribution<ctype>;

    using XAngleDistribution = std::normal_distribution<ctype>;
    using YAngleDistribution = std::normal_distribution<ctype>;
    using ZAngleDistribution = std::normal_distribution<ctype>;
};

/*!
 * \brief Geometry sampler for disks
 */
template< class ctype, class T >
class GeometrySampler< Disk<ctype>, T >
{
    using MajorAxisLengthDistribution = typename T::MajorAxisLengthDistribution;
    using MinorAxisLengthDistribution = typename T::MinorAxisLengthDistribution;
    using XAngleDistribution = typename T::XAngleDistribution;
    using YAngleDistribution = typename T::YAngleDistribution;
    using ZAngleDistribution = typename T::ZAngleDistribution;

    using Vector = Frackit::Vector<ctype, 3>;
    using Direction = Frackit::Direction<ctype, 3>;

public:
    //! export underlying geometry type
    using Disk = Frackit::Disk<ctype>;

    //! export traits class
    using Traits = T;

    /*!
     * \brief \todo TODO Doc me.
     */
    GeometrySampler(const MajorAxisLengthDistribution& majAxis,
                    const MinorAxisLengthDistribution& minAxis,
                    const XAngleDistribution& xAngle,
                    const YAngleDistribution& yAngle,
                    const ZAngleDistribution& zAngle)
    : generator_(std::random_device{}())
    , p_majorAxisLength_(majAxis)
    , p_minorAxisLength_(minAxis)
    , p_angle_x_(xAngle)
    , p_angle_y_(yAngle)
    , p_angle_z_(zAngle)
    {}

    /*!
     * \brief \todo TODO Doc me.
     */
    template<class PointSampler>
    Disk operator() (PointSampler& pointSampler)
    {
        auto a = p_majorAxisLength_(generator_);
        while (a <= 0.0) a = p_majorAxisLength_(generator_);

        auto b = p_minorAxisLength_(generator_);
        while (b <= 0.0) b = p_minorAxisLength_(generator_);

        if (b > a) b = a;

        const auto c = pointSampler();
        const auto alpha = p_angle_x_(generator_);
        const auto beta = p_angle_y_(generator_);
        const auto gamma = p_angle_z_(generator_);

        // find major/minor axis by rotations
        std::vector<Vector> axes({Vector(1.0, 0.0, 0.0),
                                  Vector(0.0, 1.0, 0.0)});

        auto e1 = Direction(axes[0]);
        auto e2 = Direction(axes[1]);
        auto e3 = Direction(Vector(0.0, 0.0, 1.0));

        rotate(axes[1], e1, alpha); // rotate minor axis around x
        rotate(axes, e2, beta);     // rotate both axes around y

        // minor axis is rotated around y in clockwise direction,
        // therefore fix sign. TODO: Check why this is the case.
        axes[1] = Vector(-1.0*axes[1].x(), axes[1].y(), axes[1].z());

        // rotate both axes around z
        rotate(axes, e3, gamma);

        return Disk(c, Direction(axes[0]), Direction(axes[1]), a, b);
    }

 private:
    std::default_random_engine generator_;

    MajorAxisLengthDistribution p_majorAxisLength_;
    MinorAxisLengthDistribution p_minorAxisLength_;
    XAngleDistribution p_angle_x_;
    YAngleDistribution p_angle_y_;
    ZAngleDistribution p_angle_z_;
};

} // end namespace Frackit

#endif // FRACKIT_DISK_SAMPLER_HH
