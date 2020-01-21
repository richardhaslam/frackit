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
 * \brief Class that defines the interface for
 *        sampler classes of geometries.
 */
#ifndef FRACKIT_GEOMETRY_SAMPLER_HH
#define FRACKIT_GEOMETRY_SAMPLER_HH

#include <frackit/common/extractctype.hh>
#include <frackit/common/extractdimension.hh>
#include <frackit/geometry/point.hh>

namespace Frackit {

/*!
 * \brief Interface for geometry sampler classes.
 *        The interface is defined by a () operator that
 *        receives a point, which is the point around which
 *        the geometry is to be created, and it returns an
 *        object of the geometry.
 */
template<class Geometry>
class GeometrySampler
{
    using ctype = typename CoordinateTypeTraits<Geometry>::type;
    static constexpr int worldDim = DimensionalityTraits<Geometry>::worldDimension();

public:
    /*!
     * \brief Creates an object of Geometry around the given point.
     * \param point The point around which the geometry is to be created
     */
    virtual Geometry operator() (const Point<ctype, worldDim>& point) = 0;
};

} // end namespace Frackit

#endif // FRACKIT_GEOMETRY_SAMPLER_HH
