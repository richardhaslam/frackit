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
 * \brief Class that describes the lateral
 *        surface of a cylinder in 3d space.
 */
#ifndef FRACKIT_CYLINDER_SURFACE_HH
#define FRACKIT_CYLINDER_SURFACE_HH

#include <cmath>
#include <cassert>

#include <frackit/precision/precision.hh>
#include "point.hh"
#include "segment.hh"
#include "direction.hh"
#include "vector.hh"
#include "circle.hh"
#include "cylinder.hh"
#include "plane.hh"

namespace Frackit {

//! Forward declaration
template<class CT> class Cylinder;

/*!
 * \brief Class that describes the lateral
 *        surface of a cylinder in 3d space.
 * \tparam CT The type used for coordinates
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
    using Plane = Frackit::Plane<ctype, 3>;

    /*!
     * \brief Constructor.
     * \param radius The cylinder radius
     * \param height The cylinder height.
     */
    CylinderSurface(ctype radius, ctype height)
    : height_(height)
    , top_(Point({0.0, 0.0, height}), Direction(Vector(0.0, 0.0, 1.0)), radius)
    , bottom_(Point({0.0, 0.0, 0.0}), Direction(Vector(0.0, 0.0, 1.0)), radius)
    {}

    /*!
     * \brief Constructor.
     * \param bottom The circle describing the
     *               rim of the bottom boundary
     * \param height The cylinder height.
     */
    CylinderSurface(const Circle& bottom, ctype height)
    : height_(height)
    , top_(makeTopCircle_(bottom, height))
    , bottom_(bottom)
    {}

    //! Return the name of this geometry.
    static std::string name() { return "CylindricalSurface"; }

    //! Return the first basis vector (horizontal axis of the cylinder)
    const Direction& base1() const { return bottom_.base1(); }
    //! Return the second basis vector (horizontal axis of the cylinder)
    const Direction& base2() const { return bottom_.base2(); }
    //! Return the third basis vector (vertical axis of the cylinder)
    const Direction& base3() const { return direction(); }
    //! Return the vertical axis of the cylinder
    const Direction& direction() const { return bottom_.normal(); }

    //! Return the height of the cylinder
    ctype height() const { return height_; }
    //! Return the radius of the cylinder
    ctype radius() const { return bottom_.radius(); }

    //! Return the circle describing the top boundary of this surface
    const Circle& upperBoundingCircle() const { return top_; }
    //! Return the circle describing the bottom boundary of this surface
    const Circle& lowerBoundingCircle() const { return bottom_; }

    //! Return the segment describing the center of the cylinder
    Segment centerSegment() const
    { return Segment(bottom_.center(), top_.center()); }

    //! Return the cylinder that we are is the lateral surface of
    Cylinder cylinder() const
    { return Cylinder(bottom_, height_); }

    /*!
     * \brief Returns the plane that describes the
     *        tangent in a point p.
     * \param p The point at which to compute the tangent plane.
     * \note This function does not check if the point lies on
     *       the cylinder surface.
     * \note If the point lies on the upper or lower rim of the
     *       cylinder surface, the tangent plane is not uniquely
     *       defined. Here, we always return the tangent of the
     *       infinite cylinder surface.
     */
    Plane getTangentPlane(const Point& p) const
    {
        assert(contains(p, Precision<ctype>::confusion()*0.5*(radius() + height())));
        Direction n( Vector(centerSegment().supportingLine().projection(p), p) );
        return Plane(p, n);
    }

    /*!
     * \brief Returns true if a point lies on the surface.
     * \param p The point to be checked
     * \param eps The epsilon (tolerance) value to be used
     */
    bool contains(const Point& p, ctype eps) const
    {
        const auto segment = centerSegment();
        const auto proj = segment.supportingLine().projection(p);

        if (!centerSegment().contains(proj, eps))
            return false;

        using std::abs;
        return abs(Vector(proj, p).length() - radius()) < eps;
    }

    /*!
     * \brief Returns true if a point lies on the surface.
     * \param p The point to be checked
     * \note This overload uses a default epsilon.
     */
    bool contains(const Point& p) const
    { return contains( p, Precision<ctype>::confusion()*0.5*(radius() + height()) ); }

private:
    //! Creates the top circle for the given bottom and height
    Circle makeTopCircle_(const Circle& bottom, ctype height)
    {
        auto n = Vector(bottom.normal());
        n *= height;

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
