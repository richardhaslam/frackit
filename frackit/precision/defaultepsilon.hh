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
 * \brief Defines epsilons to be used for floating point
 *        arithmetic on geometries, e.g. to determine if
 *        a point lies on a geometry or for computing the
 *        intersection of two geometries. The values for
 *        the epsilons are computed based on the magnitude
 *        of the provided geometry.
 */
#ifndef FRACKIT_DEFAULT_EPSILON_HH
#define FRACKIT_DEFAULT_EPSILON_HH

#include <frackit/geometry/segment.hh>
#include <frackit/geometry/circle.hh>
#include <frackit/geometry/ellipse.hh>
#include <frackit/geometry/disk.hh>
#include <frackit/geometry/cylindersurface.hh>

#include "precision.hh"

namespace Frackit {

/*!
 * \brief Default epsilon for geometries.
 *        We use the base epsilon here and provide
 *        overloads for geometries having a magnitude.
 */
template<class Geometry>
typename Geometry::ctype defaultEpsilon(const Geometry& geom)
{ return Precision<typename Geometry::ctype>::confusion(); }

/*!
 * \brief Default epsilon for operations on segments.
 */
template<class ctype, int worldDim>
ctype defaultEpsilon(const Segment<ctype, worldDim>& seg)
{ return Precision<ctype>::confusion()*seg.length(); }

/*!
 * \brief Default epsilon for operations on circles.
 */
template<class ctype, int worldDim>
ctype defaultEpsilon(const Circle<ctype, worldDim>& circle)
{ return Precision<ctype>::confusion()*circle.radius(); }

/*!
 * \brief Default epsilon for operations on ellipses.
 */
template<class ctype, int worldDim>
ctype defaultEpsilon(const Ellipse<ctype, worldDim>& ellipse)
{
    return Precision<ctype>::confusion()
           *0.5*(ellipse.majorAxisLength() + ellipse.minorAxisLength());
}

/*!
 * \brief Default epsilon for operations on disks.
 */
template<class ctype>
ctype defaultEpsilon(const Disk<ctype>& disk)
{ return defaultEpsilon(disk.boundingEllipse()); }

/*!
 * \brief Default epsilon for operations on cylinder surfaces.
 */
template<class ctype>
ctype defaultEpsilon(const CylinderSurface<ctype>& cylSurface)
{
    return Precision<ctype>::confusion()
           *0.5*(cylSurface.radius() + cylSurface.height());
}

} // end namespace Frackit

#endif // FRACKIT_DEFAULT_EPSILON_HH
