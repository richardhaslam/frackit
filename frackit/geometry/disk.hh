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
#ifndef FRACKIT_DISK_HH
#define FRACKIT_DISK_HH

#include <cmath>

#include <frackit/precision/precision.hh>
#include "ellipticalgeometry.hh"
#include "ellipse.hh"
#include "vector.hh"

namespace Frackit {

/*!
 * \brief \todo TODO doc me.
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
     * \brief \todo TODO doc me.
     * \todo check if points are on ellipse
     */
    Disk(const Ellipse& ellipse)
    : ParentType(ellipse.center(),
                 ellipse.majorAxisLength(),
                 ellipse.minorAxisLength(),
                 ellipse.majorAxis(),
                 ellipse.minorAxis())
    {}

    //! \todo TODO doc me.
    static std::string name() { return "Disk"; }

    /*!
     * \brief Returns true if a point is on the disk
     * \param p The point to be checked
     * \param eps The tolerance to be used
     * \param checkIfOnPlane Flag that can be set to false in case
     *                       it is known that the point is in-plane.
     * \todo note about choice of eps
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
     * \todo note about choice of eps
     */
    bool contains(const Point& p, bool checkIfOnPlane = true) const
    { return contains(p, Precision<ctype>::confusion(), checkIfOnPlane); }

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
