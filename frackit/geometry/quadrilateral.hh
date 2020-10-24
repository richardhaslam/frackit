// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup Geometry
 * \brief Classes that implement quadrilaterals in n-dimensional space.
 */
#ifndef FRACKIT_GEOMETRY_QUADRILATERAL_HH
#define FRACKIT_GEOMETRY_QUADRILATERAL_HH

#include <cmath>
#include <cassert>
#include <string>
#include <array>
#include <stdexcept>

#include <frackit/precision/precision.hh>

#include "geometry.hh"
#include "point.hh"
#include "vector.hh"
#include "segment.hh"
#include "triangle.hh"
#include "plane.hh"

namespace Frackit {

/*!
 * \ingroup Geometry
 * \brief Class that implements a quadrilateral
 *        in a coordinate space with the dimension worldDim.
 * \tparam CT The type used for coordinates
 * \tparam worldDim The dimension of the coordinate space
 */
template<class CT, int worldDim>
class Quadrilateral;

/*!
 * \ingroup Geometry
 * \brief Implementation of a quadrilateral in 3d space.
 * \tparam CT The type used for coordinates
 * \note This class describes planar quadrilaterals.
 *       Curved quadrilaterals are not supported.
 *       Upon construction it is checked if the provided
 *       corners lie on a common plane.
 */
template<class CT>
class Quadrilateral<CT, 3> : public Geometry
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

    //! Default constructor
    Quadrilateral() = default;

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

        const auto t1 = Triangle(p0, p1, p3);
        const auto t2 = Triangle(p0, p3, p2);

        // compute area
        area_ = t1.area() + t2.area();

        // compute center
        center_ = Point(0.0, 0.0, 0.0);
        auto d1 = Vector(center_, t1.center()); d1 *= t1.area();
        auto d2 = Vector(center_, t2.center()); d2 *= t2.area();
        center_ += d1; center_ += d2; center_ /= area_;
    }

    //! Return the name of this geometry
    std::string name() const override { return "Quadrilateral_3d"; }
    //! Return the number of corners
    static constexpr std::size_t numCorners() { return 4; }
    //! Return the i-th corner
    const Point& corner(unsigned int i) const
    {
        assert(i < numCorners());
        return corners_[i];
    }

    //! Return the number of edges
    static constexpr std::size_t numEdges() { return 4; }
    //! Return the i-th edge
    Segment edge(unsigned int i) const
    {
        assert(i < numEdges());
        switch (i)
        {
            case 0: return Segment(corners_[0], corners_[1]);
            case 1: return Segment(corners_[1], corners_[3]);
            case 2: return Segment(corners_[0], corners_[2]);
            case 3: return Segment(corners_[2], corners_[3]);
            default: throw std::runtime_error("Invalid edge index");
        }
    }

    //! Return the center of the quadrilateral
    const Point& center() const { return center_; }
    //! Return the plane this quadrilateral is embedded in
    const Plane& supportingPlane() const { return supportPlane_; }
    //! Return the area of the quadrilateral
    ctype area() const { return area_; }

    /*!
     * \brief Returns true if a point is on the quadrilateral
     * \param p The point to be checked
     * \param eps The tolerance to be used
     * \param checkIfOnPlane This can be set to false in case one
     *                       knows that the point lies on the supporting
     *                       plane in order to skip this check at this point.
     */
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
        using std::abs;
        return abs(area() - otherArea) < eps*eps;
    }

    /*!
     * \brief Returns true if a point is on the quadrilateral
     * \param p The point to be checked
     * \param checkIfOnPlane This can be set to false in case one
     *                       knows that the point lies on the supporting
     *                       plane in order to skip this check at this point.
     * \note This overload uses a default epsilon
     */
    bool contains(const Point& p,  bool checkIfOnPlane = true) const
    {
        auto eps = Precision<ctype>::confusion();
        eps *= Vector(corner(0), corner(3)).length();
        return contains(p, eps, checkIfOnPlane);
    }

private:
    std::array<Point, 4> corners_;
    Plane supportPlane_;
    Point center_;
    ctype area_;
};

} // end namespace Frackit

#endif // FRACKIT_GEOMETRY_QUADRILATERAL_HH
