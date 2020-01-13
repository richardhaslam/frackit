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
#ifndef FRACKIT_GEOMETRY_QUADRILATERAL_HH
#define FRACKIT_GEOMETRY_QUADRILATERAL_HH

#include <cassert>
#include <array>
#include <stdexcept>

#include <frackit/common/math.hh>
#include <frackit/precision/precision.hh>

#include "point.hh"
#include "vector.hh"
#include "segment.hh"
#include "triangle.hh"
#include "plane.hh"

namespace Frackit {

/*!
 * \brief \todo TODO doc me.
 */
template<class CT, int worldDim>
class Quadrilateral;

/*!
 * \brief \todo TODO doc me.
 */
template<class CT>
class Quadrilateral<CT, 3>
{
    using Triangle = Frackit::Triangle<CT, 3>;
    using Vector = Frackit::Vector<CT, 3>;

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
     * \brief The constructor.
     * \note The node & edge ordering is as follows
     *          3
     *       2 --- 3
     *    2  |     | 1
     *       0 --- 1
     *          0
     */
    Quadrilateral(const Point& p0,
                  const Point& p1,
                  const Point& p2,
                  const Point& p3)
    : corners_({p0, p1, p2, p3})
    , supportPlane_(p0, p1, p2)
    {
        assert( supportPlane_.contains(p3, Precision<ctype>::confusion()
                                           *Vector(p0, p3).length()) );

        // compute center
        center_ = Point(0.0, 0.0, 0.0);
        Vector d1(center_, p0); d1 *= 0.25;
        Vector d2(center_, p1); d2 *= 0.25;
        Vector d3(center_, p2); d3 *= 0.25;
        Vector d4(center_, p3); d4 *= 0.25;

        center_ += d1;
        center_ += d2;
        center_ += d3;
        center_ += d4;

        // compute area
        area_ = Triangle(p0, p1, p3).area() + Triangle(p0, p3, p2).area();
    }

    //! \todo TODO doc me.
    static std::string name() { return "Quadrilateral"; }

    //! \todo TODO doc me.
    static constexpr std::size_t numCorners() { return 4; }
    //! \todo TODO doc me.
    const Point& corner(unsigned int i) const
    {
        assert(i < numCorners());
        return corners_[i];
    }

    //! \todo TODO doc me.
    static constexpr std::size_t numEdges() { return 4; }
    //! \todo TODO doc me.
    Segment edge(unsigned int i) const
    {
        assert(i < numEdges());
        switch (i)
        {
            case 0: return Segment(corners_[0], corners_[1]);
            case 1: return Segment(corners_[1], corners_[3]);
            case 2: return Segment(corners_[0], corners_[2]);
            case 3: return Segment(corners_[2], corners_[3]);
            default: throw std::runtime_error(std::string("Invalid edge index"));
        }
    }

    //! \todo TODO doc me.
    const Plane& supportingPlane() const { return supportPlane_; }
    //! \todo TODO doc me.
    ctype area() const { return area_; }

    //! Returns true if a point lies on the quadrilateral
    bool contains(const Point& p, ctype eps, bool checkIfOnPlane = true) const
    {
        if (checkIfOnPlane)
            if (!supportPlane_.contains(p, eps))
                return false;

        // sum up the areas of sub-triangles with the given point
        ctype otherArea(0.0);
        otherArea += Triangle(corners_[0], corners_[1], p).area();
        otherArea += Triangle(corners_[1], corners_[3], p).area();
        otherArea += Triangle(corners_[3], corners_[2], p).area();
        otherArea += Triangle(corners_[2], corners_[0], p).area();

        // if area is equal to the quadrilateral area, point is on it
        return abs(area() - otherArea) < eps*eps;
    }

    //! Returns true if a point lies on the quadrilateral
    bool contains(const Point& p,  bool checkIfOnPlane = true) const
    {
        auto eps = Precision<ctype>::confusion();
        eps *= Vector(corner(0), corner(3)).length();
        return contains(p, eps);
    }

private:
    std::array<Point, 4> corners_;
    Plane supportPlane_;
    Point center_;
    ctype area_;
};

} // end namespace Frackit

#endif // FRACKIT_GEOMETRY_QUADRILATERAL_HH