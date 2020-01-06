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
#ifndef FRACKIT_GEOMETRY_LINE_HH
#define FRACKIT_GEOMETRY_LINE_HH

#include <frackit/precision/precision.hh>

#include "point.hh"
#include "direction.hh"
#include "vector.hh"

namespace Frackit {

// Forward declarations
template<class CT, int wd> class Point;
template<class CT, int wd> class Direction;

/*!
 * \brief \todo TODO doc me.
 */
template<class CT, int wd>
class Line
{
    using Vector = Frackit::Vector<CT, wd>;

public:
    //! export dimensionality
    static constexpr int myDimension() { return 1; };
    static constexpr int worldDimension() { return wd; };

    //! export type used for coordinates
    using ctype = CT;

    //! export underlying geometry types
    using Point = typename Frackit::Point<ctype, wd>;
    using Direction = typename Frackit::Direction<ctype, wd>;

    /*!
     * \brief \todo TODO doc me.
     */
    Line(const Point& p, const Direction& dir)
    : supportPoint_(p)
    , direction_(dir)
    {}

    //! \todo TODO doc me.
    static std::string name() { return "Line"; }
    //! \todo TODO doc me.
    const Point& supportingPoint() const { return supportPoint_; }
    //! \todo TODO doc me.
    const Direction& direction() const { return direction_; }

    //! Returns the projection of p onto the line
    Point projection(const Point& p) const
    {
        auto d = Vector(direction());
        d *= d*Vector(supportPoint_, p);

        auto result = supportPoint_;
        result += d;

        return result;
    }

    //! Returns true if a point lies on the line
    bool contains(const Point& p, ctype eps) const
    { return Vector(p, projection(p)).squaredLength() < eps*eps; }

    //! Returns true if a point lies on the line (default eps)
    bool contains(const Point& p) const
    { return contains(p, Precision<ctype>::confusion()); }

private:
    Point supportPoint_;
    Direction direction_;
};

} // end namespace Fracit

#endif // FRACKIT_GEOMETRY_LINE_HH
