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
#ifndef FRACKIT_CIRCLE_HH
#define FRACKIT_CIRCLE_HH

#include <cmath>
#include <frackit/common/math.hh>
#include <frackit/geometry/precision.hh>

#include "vector.hh"
#include "ellipticalgeometry.hh"

namespace Frackit {

/*!
 * \brief \todo TODO doc me.
 */
template<class ctype, int wd>
class Circle;

/*!
 * \brief \todo TODO doc me.
 */
template<class ctype>
class Circle<ctype, /*worldDim=*/3>
: public EllipticalGeometry<ctype, /*worldDim=*/3>
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
     * \brief \todo TODO doc me.
     */
    Circle(const Point& center,
           const Direction& normal,
           ctype radius)
    : ParentType(makeBaseGeometry_(center, normal, radius))
    {}

    //! \todo TODO doc me.
    static std::string name() { return "Circle"; }

    //! \todo TODO doc me.
    const Direction& base1() const { return this->majorAxis(); }
    const Direction& base2() const { return this->minorAxis(); }

    //! \todo TODO doc me.
    ctype radius() const { return this->majorAxisLength(); }

    //! \todo TODO doc me.
    //! \todo note about choice of eps
    bool contains(const Point& p, ctype eps) const
    {
        if (!this->supportingPlane().contains(p, eps))
            return false;

        using std::abs;
        return abs(Vector(p, this->center()).length() - radius()) < eps;
    }

    //! \todo TODO doc me.
    //! \todo note about choice of eps
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
