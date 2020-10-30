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
 * \ingroup GeometryUtilities
 * \brief Helper structs to check whether a geometry is planar. This is
 *        true for non-curved two-dimensional geometries embedded in 3d space.
 */
#ifndef FRACKIT_GEOMETRY_UTILITIES_ISPLANAR_HH
#define FRACKIT_GEOMETRY_UTILITIES_ISPLANAR_HH

#include <TopoDS_Face.hxx>
#include <TopoDS_Shell.hxx>

#include <frackit/geometry/cylindersurface.hh>

namespace Frackit {

/*!
 * \ingroup GeometryUtilities
 * \brief Helper structs to check whether a geometry is planar.
 *        Per default, two-dimensional geometries embedded in 3d space
 *        are considered to be planar.
 */
template<class G>
struct IsPlanarGeometry
{
    static constexpr bool value = DimensionalityTraits<G>::geometryDimension() == 2
                                  && DimensionalityTraits<G>::worldDimension() == 3;
};

//! Cylinder surfaces are curved
template<class ct>
struct IsPlanarGeometry<CylinderSurface<ct>> { static constexpr bool value = false; };

//! TopoDS_Face is generic, so it can represent a planar geometry,
//! but does not necessarily do so. Thus, we export false here.
template<>
struct IsPlanarGeometry<TopoDS_Face>
{ static constexpr bool value = false; };

//! TopoDS_Shell is generic, so it can represent a planar geometry,
//! but does not necessarily do so. Thus, we export false here.
template<>
struct IsPlanarGeometry<TopoDS_Shell>
{ static constexpr bool value = false; };

/*!
 * \ingroup GeometryUtilities
 * \brief Returns true if a geometry is planar.
 */
template<class G>
constexpr bool isPlanarGeometry(const G& g) { return IsPlanarGeometry<G>::value; }

} // end namespace Frackit

#endif // FRACKIT_GEOMETRY_UTILITIES_ISPLANAR_HH
