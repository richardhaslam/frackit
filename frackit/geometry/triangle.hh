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
 * \brief Class that implements triangles in n-dimensional space.
 */
#ifndef FRACKIT_GEOMETRY_TRIANGLE_HH
#define FRACKIT_GEOMETRY_TRIANGLE_HH

#include <cassert>
#include <stdexcept>
#include <string>
#include <array>

#include <frackit/common/math.hh>
#include "geometry.hh"
#include "point.hh"
#include "segment.hh"
#include "vector.hh"
#include "direction.hh"
#include "plane.hh"

namespace Frackit {

/*!
 * \brief Class that implements triangles
 *        in a space of dimension worldDim
 * \tparam CT The type used for coordinates
 * \tparam worldDim The dimension of the coordinate space
 */
template<class CT, int worldDim>
class Triangle;

/*!
 * \brief Class that implements triangles in 3d space.
 * \tparam CT The type used for coordinates
 */
template<class CT>
class Triangle<CT, 3> : public Geometry
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
     * \brief Constructor from the corners.
     * \param p1 The first corner point
     * \param p2 The second corner point
     * \param p3 The third corner point
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

    //! Return the name of the geometry
    std::string name() const override { return "Triangle_3d"; }
    //! Return the triangle area
    ctype area() const { return area_; }
    //! Return center of the triangle
    const Point& center() const { return center_; }
    //! Return the unit normal vector
    const Direction& normal() const { return normal_; }
    //! Return the plane this quadrilateral is embedded in
    Plane supportingPlane() const { return Plane(corner(0), normal()); }

    //! Return the number of corners
    static constexpr std::size_t numCorners() { return 3; }
    //! Return the i-th corner
    const Point& corner(unsigned int i) const
    {
        assert(i < numCorners());
        return corners_[i];
    }

    //! Return the number of edges
    static constexpr std::size_t numEdges() { return 3; }
    //! Return the i-th edge
    Segment edge(unsigned int i) const
    {
        assert(i < numEdges());
        switch (i)
        {
            case 0: return Segment(corners_[0], corners_[1]);
            case 1: return Segment(corners_[1], corners_[2]);
            case 2: return Segment(corners_[2], corners_[0]);
            default: throw std::runtime_error("Invalid edge index");
        }
    }

    /*!
     * \brief Returns true if a point is on the triangle
     * \param p The point to be checked
     * \param eps The tolerance to be used
     * \param checkIfOnPlane This can be set to false in case one
     *                       knows that the point lies on the supporting
     *                       plane in order to skip this check at this point.
     */
    bool contains(const Point& p, ctype eps, bool checkIfOnPlane = true) const
    {
        if (checkIfOnPlane)
            if (!supportingPlane().contains(p, eps))
                return false;

        // sum up the areas of sub-triangles with the given point
        ctype otherArea(0.0);
        otherArea += Triangle<ctype, 3>(corners_[0], corners_[1], p).area();
        otherArea += Triangle<ctype, 3>(corners_[1], corners_[2], p).area();
        otherArea += Triangle<ctype, 3>(corners_[2], corners_[0], p).area();

        // if area is equal to the triangle area, point is on it
        return abs(area() - otherArea) < eps*eps;
    }

    /*!
     * \brief Returns true if a point is on the triangle
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
        return contains(p, eps);
    }

private:
    std::array<Point, 3> corners_;
    Point center_;
    Direction normal_;
    ctype area_;
};

} // end namespace Frackit

#endif // FRACKIT_GEOMETRY_TRIANGLE_HH
