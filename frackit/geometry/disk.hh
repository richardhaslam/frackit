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
 * \brief Class that implements elliptical disks in 3d space.
 */
#ifndef FRACKIT_DISK_HH
#define FRACKIT_DISK_HH

#include <cmath>

#include <frackit/precision/precision.hh>
#include "ellipticalgeometry.hh"
#include "ellipse.hh"
#include "vector.hh"

namespace Frackit {

/*!
 * \brief Class that implements elliptical disks in 3d space.
 * \tparam CT The type used for coordinates
 */
template<class CT>
class Disk
: public EllipticalGeometry<CT, /*worldDim=*/3>
{
    using ParentType = EllipticalGeometry<CT, /*worldDim=*/3>;
    using Vector = Frackit::Vector<CT, 3>;

public:
    //! export dimensionality
    static constexpr int myDimension() { return 2; };
    static constexpr int worldDimension() { return 3; };

    //! export coordinate type
    using typename ParentType::ctype;

    //! export underlying geometry types
    using typename ParentType::Point;
    using typename ParentType::Direction;
    using Ellipse = Frackit::Ellipse<CT, 3>;

    //! pull up base class' constructor
    using ParentType::ParentType;

    /*!
     * \brief Constructor from an ellipse
     */
    Disk(const Ellipse& ellipse)
    : ParentType(ellipse.center(),
                 ellipse.majorAxis(),
                 ellipse.minorAxis(),
                 ellipse.majorAxisLength(),
                 ellipse.minorAxisLength())
    {}

    //! Return the name of this geometry
    static std::string name() { return "Disk"; }
    //! Return the disk area
    ctype area() const
    {
        return M_PI*this->majorAxisLength()
                   *this->minorAxisLength();
    }

    /*!
     * \brief Returns true if a point is on the disk
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

        // check if point is inside ellipse
        Vector deltaX(this->center(), p);
        const auto a = this->majorAxisLength();
        const auto b = this->minorAxisLength();
        const auto x = deltaX*Vector(this->majorAxis()); // x-coord in local basis
        const auto y = deltaX*Vector(this->minorAxis()); // y-coord in local basis

        return x*x/(a*a) + y*y/(b*b) <= 1.0 + eps*eps;
    }

    /*!
     * \brief Returns true if a point is on the disk
     * \param p The point to be checked
     * \param checkIfOnPlane Flag that can be set to false in case
     *                       it is known that the point is in-plane.
     * \note This overload uses a default epsilon.
     */
    bool contains(const Point& p, bool checkIfOnPlane = true) const
    { return contains(p, Precision<ctype>::confusion(), checkIfOnPlane); }

    /*!
     * \brief Returns the point on the disk for the given parameters
     * \param param1 The first parameter (angular fraction -> 0 <= param1 <= 1)
     * \param param2 The second parameter (radial fraction -> 0 <= param1 <= 1.0)
     */
    Point getPoint(ctype param1, ctype param2) const
    {
        assert(param1 >= 0.0 && param1 <= 1.0);
        assert(param2 >= 0.0 && param2 <= 1.0);

        const auto p = boundingEllipse().getPoint(param1);
        Vector d(this->center(), p);
        d *= param2;
        return p + d;
    }

    //! Returns the ellipse describing the disk's boundary
    Ellipse boundingEllipse() const
    {
        return Ellipse(this->center(),
                       this->majorAxis(),
                       this->minorAxis(),
                       this->majorAxisLength(),
                       this->minorAxisLength());
    }
};

} // end namespace OpenFrack

#endif // FRACKIT_DISK_HH
