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
#ifndef FRACKIT_PLANE_HH
#define FRACKIT_PLANE_HH

#include <cassert>

#include <frackit/common/math.hh>
#include <frackit/geometry/precision.hh>

#include "point.hh"
#include "vector.hh"
#include "segment.hh"
#include "line.hh"

namespace Frackit {

/*!
 * \brief \todo TODO doc me.
 */
template<class CT, int wd>
class Plane;

/*!
 * \brief \todo TODO doc me.
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
     * \brief \todo TODO doc me.
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
     * \brief \todo TODO doc me.
     */
    Plane(const Point& p1, const Point& p2, const Point& p3)
    : supportPoint_(p1)
    {
        base1_ = Vector(p1, p2);
        normal_ = crossProduct(Vector(base1_), Vector(p1, p3));
        base2_ = crossProduct(Vector(normal_), Vector(base1_));
    }

    /*!
     * \brief \todo TODO doc me.
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

    //! \todo TODO doc me.
    static std::string name() { return "Plane"; }
    //! \todo TODO doc me.
    const Point& supportingPoint() const { return supportPoint_; }
    //! \todo TODO doc me.
    const Direction& normal() const { return normal_; }
    //! \todo TODO doc me.
    const Direction& base1() const { return base1_; }
    //! \todo TODO doc me.
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

    //! Returns true if a point lies on the plane
    //! \todo note about choice of eps
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

#endif // FRACKIT_PLANE_HH
