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
 * \brief Classes describing circles in n-dimensional space.
 * \note For the moment, only circles in three-dimensional space are implemented.
 */
#ifndef FRACKIT_CIRCLE_HH
#define FRACKIT_CIRCLE_HH

#include <cmath>
#include <string>

#include <frackit/common/math.hh>
#include <frackit/precision/precision.hh>

#include "geometry.hh"
#include "ellipticalgeometry.hh"
#include "vector.hh"

namespace Frackit {

/*!
 * \brief Class describing a circle
 * \tparam ctype The type used for coordinates.
 * \tparam wd the space dimension.
 */
template<class ctype, int wd>
class Circle;

/*!
 * \brief Class describing a circle in 3d space.
 * \tparam ctype The type used for coordinates
 */
template<class ctype>
class Circle<ctype, /*worldDim=*/3>
: public Geometry
, public EllipticalGeometry<ctype, /*worldDim=*/3>
{
    using ParentType = EllipticalGeometry<ctype, /*worldDim=*/3>;
    using Vector = Frackit::Vector<ctype, 3>;

public:
    //! export dimensionality
    static constexpr int myDimension() { return 1; }
    static constexpr int worldDimension() { return 3; }

    //! export underlying geometry types
    using typename ParentType::Point;
    using typename ParentType::Direction;

    /*!
     * \brief The constructor.
     * \param center The center point of the circle
     * \param normal The normal direction of the plane
     *               the circle is embedded in
     */
    Circle(const Point& center,
           const Direction& normal,
           ctype radius)
    : ParentType(makeBaseGeometry_(center, normal, radius))
    {}

    //! Return the name of this geometry
    std::string name() const override { return "Circle_3d"; }

    //! Return the basis vectors of the circle
    const Direction& base1() const { return this->majorAxis(); }
    const Direction& base2() const { return this->minorAxis(); }

    //! Return the radius of the circle
    ctype radius() const { return this->majorAxisLength(); }

    /*!
     * \brief Returns true if a point lies on the circle.
     * \param p The point to be checked
     * \param eps The epsilon (tolerance) value to be used
     */
    bool contains(const Point& p, ctype eps) const
    {
        if (!this->supportingPlane().contains(p, eps))
            return false;

        using std::abs;
        return abs(Vector(p, this->center()).length() - radius()) < eps;
    }

    /*!
     * \brief Returns true if a point lies on the circle.
     * \param p The point to be checked
     * \note This overload uses a default epsilon.
     */
    bool contains(const Point& p) const
    { return contains(p, Precision<ctype>::confusion()*radius()); }

private:
    //! constructs the base geometry
    ParentType makeBaseGeometry_(const Point& center,
                                 const Direction& normal,
                                 ctype radius)
    {
        const auto majAxis = Direction(makeOrthogonalVector(Vector(normal)));
        const auto minAxis = crossProduct(Vector(normal), Vector(majAxis));
        return ParentType(center, majAxis, minAxis, radius, radius);
    }
};

} // end namespace Frackit

#endif // FRACKIT_CIRCLE_HH
