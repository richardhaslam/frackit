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
 * \brief Class that describes a cylinder in 3d space.
 */
#ifndef FRACKIT_GEOMETRY_CYLINDER_HH
#define FRACKIT_GEOMETRY_CYLINDER_HH

#include <cmath>

#include <frackit/precision/precision.hh>
#include "point.hh"
#include "segment.hh"
#include "direction.hh"
#include "vector.hh"
#include "circle.hh"
#include "disk.hh"
#include "cylindersurface.hh"

namespace Frackit {

// forward declaration
template<class CT> class CylinderSurface;

/*!
 * \brief Class that describes a cylinder in 3d space.
 * \tparam CT The type used for coordinates
 */
template<class CT>
class Cylinder
{
    using Vector = typename Frackit::Direction<CT, 3>::Vector;

public:
    //! export dimensionality
    static constexpr int myDimension() { return 3; };
    static constexpr int worldDimension() { return 3; };

    //! export type used for coordinates
    using ctype = CT;

    //! export underlying geometry types
    using Point = Frackit::Point<ctype, 3>;
    using Segment = Frackit::Segment<ctype, 3>;
    using Direction = Frackit::Direction<ctype, 3>;
    using Circle = Frackit::Circle<ctype, 3>;
    using Disk = Frackit::Disk<ctype>;
    using CylinderSurface = Frackit::CylinderSurface<ctype>;

    /*!
     * \brief Constructor.
     * \param radius The cylinder radius
     * \param height The cylinder height.
     */
    Cylinder(ctype radius, ctype height)
    : height_(height)
    , top_(Point({0.0, 0.0, height}),
           Direction(Vector(1.0, 0.0, 0.0)),
           Direction(Vector(0.0, 1.0, 0.0)),
           radius, radius)
    , bottom_(Point({0.0, 0.0, 0.0}),
              Direction(Vector(1.0, 0.0, 0.0)),
              Direction(Vector(0.0, 1.0, 0.0)),
              radius, radius)
    {}

    /*!
     * \brief Constructor.
     * \param bottom The circle describing the
     *               rim of the bottom boundary
     * \param height The cylinder height.
     */
    Cylinder(const Circle& bottom, ctype height)
    : Cylinder(Disk(bottom.center(),
                    bottom.base1(), bottom.base2(),
                    bottom.radius(), bottom.radius()), height)
    {}

    /*!
     * \brief Constructor.
     * \param bottom The disk describing the bottom boundary
     * \param height The cylinder height.
     */
    Cylinder(const Disk& bottom, ctype height)
    : height_(height)
    , top_(makeTopDisk_(bottom, height))
    , bottom_(bottom)
    {
        using std::abs;
        const auto eps = bottom.majorAxisLength()*Precision<ctype>::confusion();
        if (abs(bottom.majorAxisLength() - bottom.minorAxisLength()) > eps)
            throw std::runtime_error("Cylinder requires circular disks as base");
    }

    //! Return the name of the geometry
    static std::string name() { return "Cylinder"; }

    //! Return the first basis vector of the cylinder (horizontal axis)
    const Direction& base1() const { return bottom_.majorAxis(); }
    //! Return the second basis vector of the cylinder (horizontal axis)
    const Direction& base2() const { return bottom_.minorAxis(); }
    //! Return the third basis vector of the cylinder (vertical axis)
    const Direction& base3() const { return direction(); }
    //! Return the vertical axis of the cylinder
    const Direction& direction() const { return bottom_.normal(); }

    //! Return the top bounding face
    const Disk& topFace() const { return top_; }
    //! Return the bottom bounding face
    const Disk& bottomFace() const { return bottom_; }
    //! Return the lateral surface
    CylinderSurface lateralFace() const
    {
        return CylinderSurface(Circle(bottom_.center(),
                                      bottom_.normal(),
                                      bottom_.majorAxisLength()), height());
    }

    //! Return the segment describing the center of the cylinder
    Segment centerSegment() const
    { return Segment(bottom_.center(), top_.center()); }

    //! Return the height of the cylinder.
    ctype height() const { return height_; }
    //! Return the radius of the cylinder.
    ctype radius() const { return bottom_.majorAxisLength(); }
    //! Return the volume of the cylinder.
    ctype volume() const { return M_PI*this->radius()*this->radius()*this->height(); }

    /*!
     * \brief Returns true if a point lies within the cylinder.
     * \param p The point to be checked
     * \param eps The epsilon (tolerance) value to be used
     */
    bool contains(const Point& p, ctype eps) const
    {
        const auto segment = centerSegment();
        const auto proj = segment.supportingLine().projection(p);

        if (!centerSegment().contains(proj, eps))
            return false;

        return Vector(proj, p).length() <= radius();
    }

    /*!
     * \brief Returns true if a point lies within the cylinder.
     * \param p The point to be checked
     * \note This overload uses a default epsilon.
     */
    bool contains(const Point& p) const
    { return contains( p, Precision<ctype>::confusion()*0.5*(radius() + height()) ); }

private:
    //! Creates the top disk for the given bottom and height
    Disk makeTopDisk_(const Disk& bottom, ctype height)
    {
        auto n = Vector(bottom.normal());
        n *= height;

        auto topCenter = bottom.center();
        topCenter += n;

        return Disk(topCenter,
                    bottom.majorAxis(), bottom.minorAxis(),
                    bottom.majorAxisLength(), bottom.minorAxisLength());
    }

    // data members
    ctype height_;

    Disk top_;
    Disk bottom_;
};

} // end namespace Frackit

#endif // FRACKIT_GEOMETRY_CYLINDER_HH
