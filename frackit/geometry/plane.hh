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
 * \brief Classes that implement planes in n-dimensional space.
 */
#ifndef FRACKIT_GEOMETRY_PLANE_HH
#define FRACKIT_GEOMETRY_PLANE_HH

#include <cassert>

#include <frackit/common/math.hh>
#include <frackit/precision/precision.hh>

#include "point.hh"
#include "vector.hh"
#include "segment.hh"
#include "line.hh"

namespace Frackit {

/*!
 * \brief Class that implements a plane
 *        in a coordinate space with dimension wd.
 * \tparam CT The type used for coordinates
 * \tparam wd The dimension of the coordinate space
 */
template<class CT, int wd>
class Plane;

/*!
 * \brief Implementation of a plane in 3d space.
 * \tparam CT The type used for coordinates
 */
template<class CT>
class Plane<CT, 3>
{
    using Vector = Frackit::Vector<CT, 3>;

public:
    //! export dimensionality
    static constexpr int myDimension() { return 2; };
    static constexpr int worldDimension() { return 3; };

    //! export type used for coordinates
    using ctype = CT;

    //! export underlying geometry types
    using Point = Frackit::Point<ctype, 3>;
    using Line = Frackit::Line<ctype, 3>;
    using Direction = Frackit::Direction<ctype, 3>;
    using Segment = Frackit::Segment<ctype, 3>;

    /*!
     * \brief Constructor
     * \param p A point on the plane
     * \param normal The normal of the plane
     */
    template<class C>
    Plane(const Frackit::Point<C, 3>& p,
          const Frackit::Direction<C, 3>& normal)
    {
        supportPoint_ = p;
        normal_ = normal;
        base1_ = makeOrthogonalVector(Vector(normal_));
        base2_ = crossProduct(Vector(normal_), Vector(base1_));
    }

    /*!
     * \brief Constructor from three points.
     * \param p1 The first point on the plane
     * \param p2 The second point on the plane
     * \param p3 The third point on the plane
     */
    Plane(const Point& p1, const Point& p2, const Point& p3)
    : supportPoint_(p1)
    {
        base1_ = Vector(p1, p2);
        normal_ = crossProduct(Vector(base1_), Vector(p1, p3));
        base2_ = crossProduct(Vector(normal_), Vector(base1_));
    }

    /*!
     * \brief Constructor specifying the
     *        basis vectors of the plane.
     * \param p A point on the plane
     * \param base1 First in-plane basis unit vector
     * \param base2 Second in-plane basis unit vector
     * \param normal The unit vector normal to the plane
     */
    Plane(const Point& p,
          const Direction& base1,
          const Direction& base2,
          const Direction& normal)
    : supportPoint_(p)
    , base1_(base1)
    , base2_(base2)
    , normal_(normal)
    {
        assert( std::abs(Vector(base1)*Vector(base2)) < 1e-6 );
        assert( std::abs(Vector(base1)*Vector(normal)) < 1e-6 );
    }

    //! Return the name of this geometry
    static std::string name() { return "Plane"; }
    //! Return the supporting point of the plane
    const Point& supportingPoint() const { return supportPoint_; }
    //! Return the unit normal vector
    const Direction& normal() const { return normal_; }
    //! Return the first in-plane basis vector
    const Direction& base1() const { return base1_; }
    //! Return the second in-plane basis vector
    const Direction& base2() const { return base2_; }

    //! Returns the projection of a point p onto the plane
    Point projection(const Point& p) const
    {
        auto d = Vector(normal());
        d *= d*Vector(supportPoint_, p);
        d *= -1.0;

        auto result = p;
        result += d;
        return result;
    }

    //! Returns the projection of a line l onto the plane
    Line projection(const Line& l) const
    {
        const auto newSupport = projection(l.supportPoint());
        auto p2 = l.supportingPoint();
        p2 += l.direction();

        return Line(newSupport, Vector(newSupport, projection(p2)));
    }

    //! Returns the projection of a segment s onto the plane
    Segment projection(const Segment& s) const
    {
        return Segment(projection(s.source()),
                       projection(s.target()));
    }

    /*!
     * \brief Returns true if a point is on the plane
     * \param p The point to be checked
     * \param eps The tolerance to be used
     */
    bool contains(const Point& p, ctype eps = Precision<ctype>::confusion()) const
    {
        auto d = Vector(supportPoint_, p);
        const auto length = d.length();
        if (length != 0.0)
            d /= length;

        using std::abs;
        return abs(d*Vector(normal_)) < Precision<ctype>::confusion();
    }

private:
    Point supportPoint_;
    Direction base1_;
    Direction base2_;
    Direction normal_;
};

} // end namespace Frackit

#endif // FRACKIT_GEOMETRY_PLANE_HH
