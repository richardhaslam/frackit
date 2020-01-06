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
#ifndef FRACKIT_CYLINDER_SURFACE_HH
#define FRACKIT_CYLINDER_SURFACE_HH

#include <cmath>

#include "precision.hh"
#include "point.hh"
#include "segment.hh"
#include "direction.hh"
#include "vector.hh"
#include "circle.hh"
#include "cylinder.hh"

namespace Frackit {

//! Forward declaration
template<class CT> class Cylinder;

/*!
 * \brief \todo TODO doc me.
 */
template<class CT>
class CylinderSurface
{
    using Vector = typename Frackit::Direction<CT, 3>::Vector;

public:
    //! export dimensionality
    static constexpr int myDimension() { return 2; };
    static constexpr int worldDimension() { return 3; };

    //! export type used for coordinates
    using ctype = CT;

    //! export underlying geometry types
    using Point = Frackit::Point<ctype, 3>;
    using Segment = Frackit::Segment<ctype, 3>;
    using Direction = Frackit::Direction<ctype, 3>;
    using Circle = Frackit::Circle<ctype, 3>;
    using Cylinder = Frackit::Cylinder<ctype>;

    /*!
     * \brief \todo TODO doc me.
     */
    CylinderSurface(ctype radius, ctype height)
    : height_(height)
    , top_(Point({0.0, 0.0, height}), Direction(Vector(0.0, 0.0, 1.0)), radius)
    , bottom_(Point({0.0, 0.0, 0.0}), Direction(Vector(0.0, 0.0, 1.0)), radius)
    {}

    /*!
     * \brief \todo TODO doc me.
     */
    CylinderSurface(const Circle& bottom, ctype height)
    : height_(height)
    , top_(makeTopCircle_(bottom))
    , bottom_(bottom)
    {}

    //! \todo TODO doc me.
    static std::string name() { return "CylindricalSurface"; }

    //! \todo TODO doc me.
    const Direction& base1() const { return bottom_.base1(); }
    //! \todo TODO doc me.
    const Direction& base2() const { return bottom_.base2(); }
    //! \todo TODO doc me.
    const Direction& direction() const { return bottom_.normal(); }

    //! \todo TODO doc me.
    const Circle& upperBoundingCircle() const { return top_; }
    //! \todo TODO doc me.
    const Circle& lowerBoundingCircle() const { return bottom_; }
    //! \todo TODO doc me.
    Segment centerSegment() const { return Segment(bottom_.center(), top_.center()); }
    //! \todo TODO doc me.
    Cylinder cylinder() const { return Cylinder(bottom_, height_); }

    //! \todo TODO doc me.
    ctype height() const { return height_; }
    //! \todo TODO doc me.
    ctype radius() const { return bottom_.radius(); }

    //! Returns true if a point lies on the surface (given tolerance)
    //! \todo note about choice of eps
    bool contains(const Point& p, ctype eps) const
    {
        const auto segment = centerSegment();
        const auto proj = segment.supportingLine().projection(p);

        if (!centerSegment().contains(proj, eps))
            return false;

        using std::abs;
        return abs(Vector(proj, p).length() - radius()) < eps;
    }

    //! Returns true if a point lies on the surface (default tolerance)
    //! \todo note about choice of eps
    bool contains(const Point& p) const
    {
        using std::min;
        const auto lengthScale = min(radius(), height());
        return contains( p, Precision<ctype>::confusion()*lengthScale );
    }

private:
    //! \todo TODO doc me
    Circle makeTopCircle_(const Circle& bottom)
    {
        auto n = Vector(bottom.normal());
        n *= height_;

        auto topCenter = bottom.center();
        topCenter += n;

        return Circle(topCenter, bottom.normal(), bottom.radius());
    }

    // data members
    ctype height_;

    Circle top_;
    Circle bottom_;
};

} // end namespace Frackit

#endif // FRACKIT_CYLINDER_SURFACE_HH
