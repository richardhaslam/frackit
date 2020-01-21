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
 * \brief Classes that implement ellipses in n-dimensional space.
 */
#ifndef FRACKIT_GEOMETRY_ELLIPSE_HH
#define FRACKIT_GEOMETRY_ELLIPSE_HH

#include <cmath>

#include <frackit/precision/precision.hh>
#include "ellipticalgeometry.hh"
#include "vector.hh"

namespace Frackit {
/*!
 * \brief Class that implements an ellipse in a
 *        space with the dimension worldDim.
 * \tparam CT The type used for coordinates
 * \tparam wd The dimension of the space
 */
template<class CT, int worldDim>
class Ellipse;

/*!
 * \brief Class that implements an ellipse in 3d space.
 * \tparam CT The type used for coordinates
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

    //! Return the name of this geometry
    static std::string name() { return "Ellipse"; }

    /*!
     * \brief Returns true if a point is on the ellipse
     * \param p The point to be checked
     * \param eps The tolerance to be used
     * \param checkIfOnPlane Flag that can be set to false in case
     *                       it is known that the point is in-plane.
     * \note It is recommended to use epsilon values independent of
     *       the size of the disk here, as the coordinates are normalized
     *       by the minor & major axis lengths.
     */
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

    /*!
     * \brief Returns true if a point is on the ellipse
     * \param p The point to be checked
     * \param checkIfOnPlane Flag that can be set to false in case
     *                       it is known that the point is in-plane.
     * \note This overload uses a default epsilon.
     */
    bool contains(const Point& p, bool checkIfOnPlane = true) const
    { return contains(p, Precision<ctype>::confusion(), checkIfOnPlane); }

    /*!
     * \brief Returns the point on the ellipse for the given parameter
     * \param param The parameter (0.0 <= param <= 1.0)
     * \note param = 0.0 corresponds to the point on the ellipse
     *       for an angle of 0 w.r.t. the ellipse-local coordinate
     *       system consisting of the major and minor axes.
     *       Similarly, param = 1.0 corresponds to 2*Pi.
     */
    Point getPoint(ctype param) const
    {
        assert(param >= 0.0 && param <= 1.0);
        return getPointFromAngle(2.0*M_PI*param);
    }

    /*!
     * \brief Returns the point on the ellipse for a given angle
     * \param angle The angle in radians
     * \note The angle is to be understood with respect to the
     *       ellipse local coordinate system consisting of the
     *       major and minor axes.
     */
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

    /*!
     * \brief Returns the parameter associated with the given point
     *        on the ellipse. This parameter describes the fraction
     *        of the ellipse curve that lies between the major axis
     *        and the point.
     * \param p The point on the ellipse
     * \param checkIfOnEllipse If set to false, it is not checked if the
     *                         provided point actually is on the ellipse.
     * \returns Returns the parameter 0 <= param <= 1.0 along the ellipse
     */
    ctype getParam(const Point& p, bool checkIfOnEllipse = true) const
    { return getAngle(p, checkIfOnEllipse)/(2.0*M_PI); }

    /*!
     * \brief Returns the angle for which the given point on the ellipse
     *        can be computed on the basis of its parametrization.
     * \param p The point on the ellipse
     * \param checkIfOnEllipse If set to false, it is not checked if the
     *                         provided point actually is on the ellipse.
     * \returns Returns the ellipse angle 0 <= angle <= 2.0*Pi
     */
    ctype getAngle(const Point& p, bool checkIfOnEllipse = true) const
    {
        if (checkIfOnEllipse)
            if (!contains(p))
                throw std::runtime_error("Point not on ellipse!");

        const auto d = Vector(this->center(), p);
        const auto x = d*Vector(this->majorAxis());
        const auto y = d*Vector(this->minorAxis());

        using std::abs;
        if (abs(x) < this->majorAxisLength()*Precision<ctype>::confusion())
            return std::signbit(y) ? 3.0*M_PI/2.0 : M_PI/2.0;
        if (abs(y) < this->minorAxisLength()*Precision<ctype>::confusion())
            return std::signbit(x) ? M_PI : 0.0;

        const auto ratio = this->majorAxisLength()/this->minorAxisLength();
        const auto isNegX = std::signbit(x);
        const auto isNegY = std::signbit(y);

        using std::atan;
        if (!isNegY && !isNegX) return atan(y*ratio/x);
        else if (isNegY && !isNegX) return 2.0*M_PI + atan(y*ratio/x);
        else return M_PI + atan(y*ratio/x);
    }
};

} // end namespace Frackit

#endif // FRACKIT_GEOMETRY_ELLIPSE_HH
