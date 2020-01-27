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
 * \brief Class to randomly generate disks in three-dimensional space.
 */
#ifndef FRACKIT_DISK_SAMPLER_HH
#define FRACKIT_DISK_SAMPLER_HH

#include <memory>
#include <random>
#include <type_traits>

#include <frackit/common/math.hh>
#include <frackit/geometry/point.hh>
#include <frackit/geometry/vector.hh>
#include <frackit/geometry/direction.hh>
#include <frackit/geometry/disk.hh>

#include "pointsampler.hh"
#include "geometrysampler.hh"

namespace Frackit {

/*!
 * \ingroup Sampling
 * \brief Default traits class for the disk sampler.
 *        Uses uniform distributions for all parameters.
 */
template<class ctype>
struct DefaultDiskSamplerTraits
{
    // Distribution used for the major axis length
    using MajorAxisLengthDistribution = std::normal_distribution<ctype>;
    // Distribution used for the minor axis length
    using MinorAxisLengthDistribution = std::normal_distribution<ctype>;

    // The following distribution determine the rotation angles
    // around the x-, y- and z- axis that have to be performed
    // subsequently (!) to rotate the standard basis of R^3 into the
    // disk local basis. The disk local basis is defined such that
    // the major axis is the first basis vector, the minor axis
    // the second basis vector and the normal vector the third.
    // We chose this approach as these three parameters can be
    // sampled independently, while the major/minor/normal axes
    // cannot be chosen independently.
    using XAngleDistribution = std::normal_distribution<ctype>;
    using YAngleDistribution = std::normal_distribution<ctype>;
    using ZAngleDistribution = std::normal_distribution<ctype>;
};

/*!
 * \ingroup Sampling
 * \brief Sampler for disk geometries.
 *        The disks are generated by sampling from provided
 *        distributions for:
 *            - major axis length
 *            - minor axis length
 *            - rotation angle of standard basis around x-axis
 *            - rotation angle of standard basis around y-axis
 *            - rotation angle of standard basis around z-axis
 *        The rotation angles are used to define the orientation of the
 *        local basis of the disk, composed of major axis, minor axis and normal.
 *        We start from the disk being defined such that its local basis is aligned
 *        with standard basis of R^3, i.e. the major axis is aligned with the x-axis,
 *        the minor axis is aligned with the y-axis and the normal of the disk is
 *        is aligned with the z-axis. Then, the angles state with which angles the
 *        standard basis of R^3 has to be rotated around the x-, y- and the z-axis,
 *        in this order (!), for it to describe the desired orientation of the
 *        disk-local basis. This approach was chosen since independent sampling of
 *        the axis vectors is not possible due to the requirement that they must form
 *        an orthonormal basis.
 */
template< class ctype = double, class T = DefaultDiskSamplerTraits<ctype> >
class DiskSampler : public GeometrySampler< Disk<ctype> >
{
    using MajorAxisLengthDistribution = typename T::MajorAxisLengthDistribution;
    using MinorAxisLengthDistribution = typename T::MinorAxisLengthDistribution;
    using XAngleDistribution = typename T::XAngleDistribution;
    using YAngleDistribution = typename T::YAngleDistribution;
    using ZAngleDistribution = typename T::ZAngleDistribution;

    using Point = Frackit::Point<ctype, 3>;
    using Vector = Frackit::Vector<ctype, 3>;
    using Direction = Frackit::Direction<ctype, 3>;

public:
    //! export underlying geometry type
    using Disk = Frackit::Disk<ctype>;
    using Geometry = Disk;

    //! export traits class
    using Traits = T;

    /*!
     * \brief Constructor.
     * \param pointSampler Samples the center points of the disks
     * \param majAxis Distribution used to sample major axis lengths
     * \param minAxis Distribution used to sample minor axis lengths
     * \param xAngle Distribution used to sample the angle of rotation around x-axis
     * \param yAngle Distribution used to sample the angle of rotation around y-axis
     * \param zAngle Distribution used to sample the angle of rotation around z-axis
     * \note For more info on the meaning of the rotation angles, see the description
     *       of this class.
     */
    template<class PointSamplerImpl>
    DiskSampler(const PointSamplerImpl& pointSampler,
                const MajorAxisLengthDistribution& majAxis,
                const MinorAxisLengthDistribution& minAxis,
                const XAngleDistribution& xAngle,
                const YAngleDistribution& yAngle,
                const ZAngleDistribution& zAngle)
    : pointSampler_(std::make_shared<PointSamplerImpl>(pointSampler))
    , generator_(std::random_device{}())
    , p_majorAxisLength_(majAxis)
    , p_minorAxisLength_(minAxis)
    , p_angle_x_(xAngle)
    , p_angle_y_(yAngle)
    , p_angle_z_(zAngle)
    {
        static_assert(std::is_base_of<PointSampler<Point>, PointSamplerImpl>::value,
                      "The provided point sampler does not inherit from the point sampler interface");
    }

    /*!
     * \brief Generate a random disk.
     */
    Disk operator() () override
    {
        auto a = p_majorAxisLength_(generator_);
        while (a <= 0.0) a = p_majorAxisLength_(generator_);

        auto b = p_minorAxisLength_(generator_);
        while (b <= 0.0) b = p_minorAxisLength_(generator_);

        if (b > a) b = a;

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

        // sample center point and make disk
        return Disk((*pointSampler_)(), Direction(axes[0]), Direction(axes[1]), a, b);
    }

 private:
    std::shared_ptr<PointSampler<Point>> pointSampler_; //!< pointer to the sampler for disk center points
    std::default_random_engine generator_;              //!< Random number generator

    MajorAxisLengthDistribution p_majorAxisLength_;     //!< Distribution used for major axis length
    MinorAxisLengthDistribution p_minorAxisLength_;     //!< Distribution used for minor axis length
    XAngleDistribution p_angle_x_;                      //!< Distribution used for x-axis rotation of local basis
    YAngleDistribution p_angle_y_;                      //!< Distribution used for y-axis rotation of local basis
    ZAngleDistribution p_angle_z_;                      //!< Distribution used for z-axis rotation of local basis
};

} // end namespace Frackit

#endif // FRACKIT_DISK_SAMPLER_HH
