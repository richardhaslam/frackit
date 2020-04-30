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
 * \brief Class that describes a hollow cylinder in 3d space.
 */
#ifndef FRACKIT_GEOMETRY_HOLLOW_CYLINDER_HH
#define FRACKIT_GEOMETRY_HOLLOW_CYLINDER_HH

#include <cmath>
#include <string>
#include <memory>

#include <frackit/precision/precision.hh>

#include "geometry.hh"
#include "point.hh"
#include "segment.hh"
#include "direction.hh"
#include "circle.hh"
#include "vector.hh"
#include "cylinder.hh"
#include "cylindersurface.hh"

namespace Frackit {

// forward declaration
template<class CT> class CylinderSurface;

/*!
 * \ingroup Geometry
 * \brief Class that describes a hollow cylinder in 3d space.
 * \tparam CT The type used for coordinates
 */
template<class CT>
class HollowCylinder : public Geometry
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
    using Cylinder = Frackit::Cylinder<ctype>;
    using CylinderSurface = Frackit::CylinderSurface<ctype>;

    /*!
     * \brief Constructor.
     * \param innerRadius The inner radius
     * \param outerRadius The outer radius.
     * \param height The height of the hollow cylinder.
     */
    HollowCylinder(ctype innerRadius, ctype outerRadius, ctype height)
    : HollowCylinder(Point(0.0, 0.0, 0.0),
                     Direction(Vector(0.0, 0.0, 1.0)),
                     innerRadius, outerRadius, height)
    {}

    /*!
     * \brief Constructor.
     * \param bottomCenter Center point of the bottom.
     * \param axis The cylinder axis.
     * \param innerRadius The inner radius
     * \param outerRadius The outer radius
     */
    HollowCylinder(const Point& bottomCenter,
                   const Direction& axis,
                   ctype innerRadius,
                   ctype outerRadius,
                   ctype height)
    : bottomInnerCircle_(bottomCenter, axis, innerRadius)
    , topInnerCircle_(bottomCenter + axis*height, axis, innerRadius)
    , bottomOuterCircle_(nullptr)
    , topOuterCircle_(nullptr)
    , innerMantle_(nullptr)
    , outerMantle_(nullptr)
    , outerRadius_(outerRadius)
    , height_(height)
    {}

    //! Return the name of the geometry
    std::string name() const override { return "HollowCylinder"; }

    //! Return the first basis vector of the cylinder (horizontal axis)
    const Direction& base1() const { return bottomInnerCircle_.majorAxis(); }
    //! Return the second basis vector of the cylinder (horizontal axis)
    const Direction& base2() const { return bottomInnerCircle_.minorAxis(); }
    //! Return the third basis vector of the cylinder (vertical axis)
    const Direction& base3() const { return axis(); }
    //! Return the vertical axis of the cylinder
    const Direction& axis() const { return bottomInnerCircle_.normal(); }

    //! Return the inner circle of the bottom face
    const Circle& bottomInnerCircle() const { return bottomInnerCircle_; }
    //! Return the inner circle of the top face
    const Circle& topInnerCircle() const { return topInnerCircle_; }

    //! Return the outer circle of the bottom face
    const Circle& bottomOuterCircle() const
    {
        if (!bottomOuterCircle_)
            bottomOuterCircle_ = std::make_unique<Circle>(bottomInnerCircle_.center(), axis(), outerRadius());
        return *bottomOuterCircle_;
    }

    //! Return the outer circle of the top face
    const Circle& topOuterCircle() const
    {
        if (!topOuterCircle_)
            topOuterCircle_ = std::make_unique<Circle>(topInnerCircle_.center(), axis(), outerRadius());
        return *topOuterCircle_;
    }

    //! Return inner mantle surface
    const CylinderSurface& innerMantle() const
    {
        if (!innerMantle_)
            innerMantle_ = std::make_unique<CylinderSurface>(bottomInnerCircle(), height());
        return *innerMantle_;
    }

    //! Return outer mantle surface
    const CylinderSurface& outerMantle() const
    {
        if (!outerMantle_)
            outerMantle_ = std::make_unique<CylinderSurface>(bottomOuterCircle(), height());
        return *outerMantle_;
    }

    //! Return the segment describing the center of the hollow cylinder
    Segment centerSegment() const
    { return Segment(bottomInnerCircle_.center(), topInnerCircle_.center()); }

    //! Return the cylinder that this hollow cylinder is a part of
    Cylinder fullCylinder() const
    {
        Circle bottom(bottomInnerCircle_.center(),
                      bottomInnerCircle_.normal(),
                      outerRadius_);
        return Cylinder(bottom, height());
    }

    //! Return the height of the hollow cylinder.
    ctype height() const { return height_; }
    //! Return the inner radius of the hollow cylinder.
    ctype innerRadius() const { return bottomInnerCircle_.radius(); }
    //! Return the outer radius of the hollow cylinder.
    ctype outerRadius() const { return outerRadius_; }
    //! Return the volume of the cylinder.
    ctype volume() const
    {
        return M_PI*height()*(outerRadius()*outerRadius()
                            - innerRadius()*innerRadius());
    }

    /*!
     * \brief Returns true if a point lies within the hollow cylinder.
     * \param p The point to be checked
     * \param eps The epsilon (tolerance) value to be used
     */
    bool contains(const Point& p, ctype eps) const
    {
        const auto segment = centerSegment();
        const auto proj = segment.supportingLine().projection(p);

        if (!centerSegment().contains(proj, eps))
            return false;

        const auto squaredDist = Vector(proj, p).squaredLength();
        const auto squaredEps = eps*eps;
        return squaredDist <= outerRadius()*outerRadius() + squaredEps
               && squaredDist >= innerRadius()*innerRadius() - squaredEps;
    }

    /*!
     * \brief Returns true if a point lies within the hollow cylinder.
     * \param p The point to be checked
     * \note This overload uses a default epsilon.
     */
    bool contains(const Point& p) const
    { return contains( p, Precision<ctype>::confusion()*0.5*(outerRadius() + height()) ); }

private:
    Circle bottomInnerCircle_;
    Circle topInnerCircle_;

    // the outer circles are constructed on demand
    mutable std::unique_ptr<Circle> bottomOuterCircle_;
    mutable std::unique_ptr<Circle> topOuterCircle_;

    // the mantle surfaces are constructed on demand
    mutable std::unique_ptr<CylinderSurface> innerMantle_;
    mutable std::unique_ptr<CylinderSurface> outerMantle_;

    ctype outerRadius_;
    ctype height_;
};

} // end namespace Frackit

#endif // FRACKIT_GEOMETRY_CYLINDER_HH
