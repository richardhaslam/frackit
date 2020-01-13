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
#ifndef FRACKIT_ELLIPSE_HH
#define FRACKIT_ELLIPSE_HH

#include <cmath>

#include <frackit/geometry/precision.hh>
#include "ellipticalgeometry.hh"
#include "vector.hh"

namespace Frackit {
/*!
 * \brief \todo TODO doc me.
 */
template<class CT, int worldDim>
class Ellipse;

/*!
 * \brief \todo TODO doc me.
 */
template<class CT>
class Ellipse<CT, /*worldDim=*/3>
: public EllipticalGeometry<CT, /*worldDim=*/3>
{
    using ParentType = EllipticalGeometry<CT, /*worldDim=*/3>;
    using Vector = Frackit::Vector<CT, /*worldDim*/3>;

public:
    //! export dimensionality
    static constexpr int myDimension() { return 1; }
    static constexpr int worldDimension() { return 3; }

    //! export type used for coordinates
    using typename ParentType::ctype;

    //! export underlying geometry types
    using typename ParentType::Point;

    //! pull up base class' constructor
    using ParentType::ParentType;

    //! \todo TODO doc me.
    static std::string name() { return "Ellipse"; }

    //! \todo TODO doc me.
    ctype length() const
    { return M_PI*(this->majorAxisLength() + this->minorAxisLength()); }

    //! Returns true if a point is on the ellipse
    //! \todo note about choice of eps
    bool contains(const Point& p, ctype eps, bool checkIfOnPlane = true) const
    {
        if (checkIfOnPlane)
            if (!this->supportingPlane().contains(p, eps))
                return false;

        // check if point is on ellipse
        Vector deltaX(this->center(), p);
        const auto a = this->majorAxisLength();
        const auto b = this->minorAxisLength();
        const auto x = deltaX*Vector(this->majorAxis()); // x-coord in local basis
        const auto y = deltaX*Vector(this->minorAxis()); // y-coord in local basis

        using std::abs;
        return abs(x*x/(a*a) + y*y/(b*b) - 1.0) < eps;
    }

    //! Returns true if a point is on the ellipse
    //! \todo note about choice of eps
    bool contains(const Point& p, bool checkIfOnPlane = true) const
    {
        using std::min;
        auto eps = min(this->majorAxisLength(), this->minorAxisLength());
        eps *= Precision<ctype>::confusion();

        return contains(p, eps, checkIfOnPlane);
    }

    //! Returns the point on the arc for the given parameter
    //! \note It has to be 0.0 <= param <= 1.0, where 0.0
    //!       corresponds to the source and 1.0 to the target.
    Point getPoint(ctype param) const
    {
        assert(param >= 0.0 && param <= 1.0);
        return getPointFromAngle(2.0*M_PI*param);
    }

    //! Returns the point on the ellipse for a given angle
    //! \todo mention positive angle requirement
    Point getPointFromAngle(ctype angle) const
    {
        assert(!std::signbit(angle));

        using std::sin;
        using std::cos;
        const auto x = cos(angle)*this->majorAxisLength();
        const auto y = sin(angle)*this->minorAxisLength();

        auto xVec = Vector(this->majorAxis()); xVec *= x;
        auto yVec = Vector(this->minorAxis()); yVec *= y;

        auto result = this->center();
        result += xVec;
        result += yVec;
        return result;
    }

protected:
    //! Returns the angle corresponding to a point on the ellipse
    ctype getAngle(const Point& p, bool checkIfOnEllipse = true) const
    {
        if (checkIfOnEllipse)
            if (!contains(p))
                throw std::runtime_error(std::string("Point not on ellipse!"));

        const auto d = Vector(this->center(), p);
        const auto x = d*Vector(this->majorAxis());
        const auto y = d*Vector(this->minorAxis());

        using std::abs;
        if (abs(x) < Precision<ctype>::angular())
        {
            if (std::signbit(y)) return 3.0*M_PI/2.0;
            else return M_PI/2.0;
        }

        if (abs(y) < Precision<ctype>::angular())
        {
            if (std::signbit(x)) return M_PI;
            else return 0.0;
        }

        using std::acos;
        if (y > 0.0) return acos(x/d.length());
        else return 2.0*M_PI - acos(x/d.length());
    }
};

} // end namespace Frackit

#endif // FRACKIT_ELLIPSE_HH
