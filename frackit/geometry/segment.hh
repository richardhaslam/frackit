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
 * \brief Class that describes a segment in n-dimensional space.
 */
#ifndef FRACKIT_GEOMETRY_SEGMENT_HH
#define FRACKIT_GEOMETRY_SEGMENT_HH

#include <string>
#include <frackit/precision/precision.hh>

#include "geometry.hh"
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
 * \ingroup Geometry
 * \brief Class that describes a segment in
 *        an n-dimensional space
 * \tparam CT The type used for coordinates.
 * \tparam wd The dimension of the coordinate space.
 */
template<class CT, int wd>
class Segment : public Geometry
{
    using ThisType = Segment<CT, wd>;
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
     * \brief Construct a segment from its corners.
     * \param source The first corner of the segment
     * \param target The second corner of the segment
     */
    Segment(const Point& source, const Point& target)
    : source_(source)
    , target_(target)
    {}

    //! Return the name of this geometry
    std::string name() const override { return "Segment_" + std::to_string(wd) + "d"; }
    //! Return the first corner of the segment
    const Point& source() const { return source_; }
    //! Return the second corner of the segment
    const Point& target() const { return target_; }
    //! Return the center of the segment
    Point center() const { return getPoint(0.5); }

    //! Constructs the direction vector
    Direction direction() const { return Vector(source(), target()); }
    //! Returns the length of the segment
    ctype length() const { return Vector(source(), target()).length(); }
    //! Constructs an object for the supporting line
    Line supportingLine() const
    { return Line(source(), direction()); }

    /*!
     * \brief Returns true if a point is on the quadrilateral
     * \param p The point to be checked
     * \param eps The tolerance to be used
     * \param checkIfOnLine This can be set to false in case one
     *                      knows that the point lies on the supporting
     *                      line in order to skip this check at this point.
     */
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

    /*!
     * \brief Returns true if a point is on the quadrilateral
     * \param p The point to be checked
     * \param checkIfOnLine This can be set to false in case one
     *                      knows that the point lies on the supporting
     *                      line in order to skip this check at this point.
     * \note This overload uses a default epsilon
     */
    bool contains(const Point& p, bool checkIfOnLine = true) const
    { return contains(p, Precision<ctype>::confusion()*length(), checkIfOnLine); }

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

    //! return the reversed segment
    Segment getReversed() const { return ThisType(target(), source()); }

private:
    Point source_;
    Point target_;
};

} // end namespace Frackit

#endif // FRACKIT_GEOMETRY_SEGMENT_HH
