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
 * \brief Class that implements spheres in 3d space.
 */
#ifndef FRACKIT_GEOMETRY_SPHERE_HH
#define FRACKIT_GEOMETRY_SPHERE_HH

#include <cassert>
#include <cmath>
#include <string>

#include <frackit/precision/precision.hh>
#include "point.hh"
#include "vector.hh"

namespace Frackit {

/*!
 * \ingroup Geometry
 * \brief Class that implements spheres in 3d space.
 * \tparam CT The type used for coordinates
 */
template<class CT>
class Sphere : public Geometry
{
public:
    //! export dimensionality
    static constexpr int myDimension() { return 3; };
    static constexpr int worldDimension() { return 3; };

    //! export type used for coordinates
    using ctype = CT;

    //! export underlying geometry types
    using Point = Frackit::Point<ctype, 3>;

    /*!
     * \brief Construct a sphere from a center an radius.
     */
    Sphere(const Point& c, ctype radius)
    : center_(c)
    , radius_(radius)
    {
        assert(radius > 0.0);
    }

    //! Return the name of this geometry.
    std::string name() const override { return "Sphere"; }

    //! Return the center point.
    const Point& center() const { return center_; }

    //! Return the radius.
    ctype radius() const { return radius_; }

    //! Return the volume of the aphere
    ctype volume() const { return 4.0/3.0*M_PI*radius()*radius()*radius(); }

    /*!
     * \brief Returns true if a point lies within the sphere.
     * \param p The point to be checked
     * \param eps The epsilon (tolerance) value to be used
     */
    bool contains(const Point& p, ctype eps) const
    {
        const auto v = Vector<ctype, 3>(center(), p);
        return v.squaredLength() - radius()*radius() < eps*eps;
    }

    /*!
     * \brief Returns true if a point lies within the sphere.
     * \param p The point to be checked
     * \note This overload uses a default epsilon.
     */
    bool contains(const Point& p) const
    { return contains(p, radius_*Precision<ctype>::confusion()); }

private:
    Point center_;
    ctype radius_;
};

} // end namespace Frackit

#endif // FRACKIT_GEOMETRY_SPHERE_HH
