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
 * \brief Base class for elliptical geometries.
 */
#ifndef FRACKIT_GEOMETRY_ELLIPTICAL_GEOMETRY_HH
#define FRACKIT_GEOMETRY_ELLIPTICAL_GEOMETRY_HH

#include <cassert>

#include <frackit/common/math.hh>

#include "point.hh"
#include "direction.hh"
#include "vector.hh"
#include "plane.hh"

namespace Frackit {

//! Forward declarations
template<class CT, int wd> class Point;
template<class CT, int wd> class Direction;
template<class CT, int wd> class Vector;
template<class CT, int wd> class Plane;

/*!
 * \ingroup Geometry
 * \brief Base class for elliptical geometries
 *        in n-dimensional space.
 * \tparam CT The type used for coordinates
 * \tparam wd The dimension of the space
 */
template<class CT, int worldDim>
class EllipticalGeometry;

/*!
 * \ingroup Geometry
 * \brief Base class for ellipticcal geometries
 *        in three-dimensional space.
 * \tparam CT The type used for coordinates
 */
template<class CT>
class EllipticalGeometry<CT, /*worldDim=*/3>
{

public:
    //! export type used for coordinates
    using ctype = CT;

    //! export underlying geometry types
    using Point = Frackit::Point<CT, 3>;
    using Direction = Frackit::Direction<CT, 3>;
    using Plane = Frackit::Plane<CT, 3>;

    //! default constructor
    EllipticalGeometry() = default;

    /*!
     * \brief Constructor.
     * \param center The center of the ellipse
     * \param majAxis The major axis
     * \param minAxis The minor axis
     * \param majAxisLength The major axis length
     * \param minAxisLength The minor axis length
     */
    EllipticalGeometry(const Point& center,
                       const Direction& majAxis,
                       const Direction& minAxis,
                       ctype majAxisLength,
                       ctype minAxisLength)
    : center_(center)
    , majorAxis_(majAxis)
    , minorAxis_(minAxis)
    , majorAxisLength_(majAxisLength)
    , minorAxisLength_(minAxisLength)
    {
        using Vector = Frackit::Vector<ctype, 3>;
        assert(std::abs(Vector(minAxis)*Vector(majAxis)) < 1e-7);
        normal_ = Direction( crossProduct(Vector(majAxis), Vector(minAxis)) );
    }

    //! Return the center of the ellipse
    const Point& center() const { return center_; }

    //! Return the major axis
    const Direction& majorAxis() const { return majorAxis_; }
    //! Return the minor axis
    const Direction& minorAxis() const { return minorAxis_; }
    //! Return the normal vector of the supporting plane
    const Direction& normal() const { return normal_; }

    //! Return the major axis length
    ctype majorAxisLength() const { return majorAxisLength_; }
    //! Return the minor axis length
    ctype minorAxisLength() const { return minorAxisLength_; }

    //! Return the plane the ellipse is embedded in
    Plane supportingPlane() const
    {
        return Plane(center(),
                     majorAxis(),
                     minorAxis(),
                     normal());
    }

private:
    Point center_;

    Direction majorAxis_;
    Direction minorAxis_;
    Direction normal_;

    ctype majorAxisLength_;
    ctype minorAxisLength_;
};

} // end namespace Frackit

#endif // FRACKIT_GEOMETRY_ELLIPTICAL_GEOMETRY_HH
