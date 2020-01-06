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
#ifndef FRACKIT_GEOMETRY_SEGMENT_HH
#define FRACKIT_GEOMETRY_SEGMENT_HH

#include <cmath>
#include <frackit/geometry/precision.hh>

#include "point.hh"
#include "vector.hh"
#include "direction.hh"
#include "line.hh"

namespace Frackit {

// Forward declarations
template<class CT, int wd> class Point;
template<class CT, int wd> class Vector;
template<class CT, int wd> class Direction;
template<class CT, int wd> class Line;

/*!
 * \brief \todo TODO doc me.
 */
template<class CT, int wd>
class Segment
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
    using Line = typename Frackit::Line<ctype, wd>;

    //! Default constructor
    Segment() = default;

    /*!
     * \brief \todo TODO doc me.
     */
    Segment(const Point& source, const Point& target)
    : source_(source)
    , target_(target)
    {}

    //! \todo TODO doc me.
    static std::string name() { return "Segment"; }
    //! \todo TODO doc me.
    const Point& source() const { return source_; }
    //! \todo TODO doc me.
    const Point& target() const { return target_; }

    //! Constructs the direction vector
    Direction direction() const { return Vector(source(), target()); }

    //! Constructs an object for the supporting line
    Line supportingLine() const
    { return Line(source(), direction()); }

    //! Returns true if a point lies on the segment
    //! \todo note about choice of eps
    bool contains(const Point& p, ctype eps, bool checkIfOnLine = true) const
    {
        const auto line = supportingLine();

        if (checkIfOnLine)
            if (!line.contains(p, eps))
                return false;

        const auto d1 = Vector(source(), p);
        const auto d2 = Vector(target(), p);

        const auto sp1 = Vector(line.direction())*d1;
        const auto sp2 = Vector(line.direction())*d2;

        if (std::signbit(sp1) != std::signbit(sp2))
            return true;

        // check if p is either source or target (for given eps)
        const auto eps2 = eps*eps;
        if (d1.squaredLength() < eps2) return true;
        if (d2.squaredLength() < eps2) return true;
        return false;
    }

    //! Returns true if a point lies on the segment
    //! \todo note about choice of eps
    bool contains(const Point& p, bool checkIfOnLine = true) const
    { return contains(p, Vector(source(), target()).length()*Precision<ctype>::confusion(), checkIfOnLine); }

    //! Returns the point on the segment for the given parameter
    //! \note It has to be 0.0 <= param <= 1.0, where 0.0
    //!       corresponds to the source and 1.0 to the target.
    Point getPoint(ctype param) const
    {
        assert(param >= 0.0 && param <= 1.0);

        auto d = Vector(source(), target());
        d *= param;

        return source() + d;
    }

private:
    Point source_;
    Point target_;
};

} // end namespace Frackit

#endif // FRACKIT_GEOMETRY_SEGMENT_HH
