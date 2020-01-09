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
 * \brief \todo TODO doc me.
 */
#ifndef FRACKIT_GEOMETRY_TRIANGLE_HH
#define FRACKIT_GEOMETRY_TRIANGLE_HH

#include <cassert>
#include <stdexcept>
#include <array>

#include <frackit/common/math.hh>
#include "point.hh"
#include "segment.hh"
#include "vector.hh"
#include "direction.hh"
#include "plane.hh"

namespace Frackit {

/*!
 * \brief \todo TODO doc me.
 */
template<class CT, int worldDim>
class Triangle;

/*!
 * \brief \todo TODO doc me.
 */
template<class CT>
class Triangle<CT, 3>
{
    using Vector = Frackit::Vector<CT, 3>;
    using Direction = Frackit::Direction<CT, 3>;

public:
    //! export dimensionality
    static constexpr int myDimension() { return 2; };
    static constexpr int worldDimension() { return 3; };

    //! export type used for coordinates
    using ctype = CT;

    //! export underlying geometry types
    using Point = Frackit::Point<ctype, 3>;
    using Segment = Frackit::Segment<ctype, 3>;
    using Plane = Frackit::Plane<ctype, 3>;

    /*!
     * \brief \todo TODO doc me.
     */
    Triangle(const Point& p1,
             const Point& p2,
             const Point& p3)
    : corners_({p1, p2, p3})
    , center_()
    {
        const auto x = 1.0/3.0*(p1.x() + p2.x() + p3.x());
        const auto y = 1.0/3.0*(p1.y() + p2.y() + p3.y());
        const auto z = 1.0/3.0*(p1.z() + p2.z() + p3.z());
        center_ = Point(x, y, z);

        const auto normalVec = crossProduct(Vector(p1, p2),
                                            Vector(p1, p3));
        normal_ = Direction( normalVec );
        area_ = 0.5*normalVec.length();
    }

    //! \todo TODO doc me.
    static std::string name() { return "Triangle"; }
    //! \todo TODO doc me.
    ctype area() const { return area_; }
    //! \todo TODO doc me.
    const Point& center() const { return center_; }
    //! \todo TODO doc me.
    const Direction& normal() const { return normal_; }

    //! \todo TODO doc me.
    static constexpr std::size_t numCorners() { return 3; }
    //! \todo TODO doc me.
    const Point& corner(unsigned int i) const
    {
        assert(i < numCorners());
        return corners_[i];
    }

    //! \todo TODO doc me.
    static constexpr std::size_t numEdges() { return 3; }
    //! \todo TODO doc me.
    Segment edge(unsigned int i) const
    {
        assert(i < numEdges());
        switch (i)
        {
            case 0: return Segment(corners_[0], corners_[1]);
            case 1: return Segment(corners_[1], corners_[2]);
            case 2: return Segment(corners_[2], corners_[0]);
            default: throw std::runtime_error(std::string("Invalid edge index"));
        }
    }

    //! \todo TODO doc me.
    Plane supportingPlane() const { return Plane(corner(0), normal()); }

private:
    std::array<Point, 3> corners_;
    Point center_;
    Direction normal_;
    ctype area_;
};

} // end namespace Frackit

#endif // FRACKIT_GEOMETRY_TRIANGLE_HH
