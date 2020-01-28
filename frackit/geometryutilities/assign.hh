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
 * \brief Utility functionality to cast pointers on objects of
 *        the geometry interface into the corresponding geometry type.
 */
#ifndef FRACKIT_GEOMETRY_UTILITY_ASSIGN_HH
#define FRACKIT_GEOMETRY_UTILITY_ASSIGN_HH

#include <memory>
#include <type_traits>

#include <frackit/geometry/geometry.hh>
#include <frackit/geometryutilities/name.hh>

namespace Frackit {

/*!
 * \ingroup GeometryUtilities
 * \brief Try to cast a pointer on a geometrical object
 *        into the provided instance of a geometry.
 * \param geoPtr Pointer to an object of the geometry interface
 * \param geo The geometry object in which it should be casted
 */
template<class GeoType>
bool assign(Geometry* geoPtr, GeoType& geo)
{
    static_assert(std::is_base_of<Geometry, GeoType>::value,
                  "Provided geometry does not inherit from the geometry interface");

    if (geometryName(geo) == geoPtr->name())
    {
        GeoType* actualGeoPtr = dynamic_cast<GeoType*>(geoPtr);
        geo = GeoType(*actualGeoPtr);
        return true;
    }

    return false;
}

/*!
 * \brief Overload for std::shared_ptr.
 */
template<class GeoType>
bool assign(std::shared_ptr<Geometry> geoPtr, GeoType& geo)
{ return assign(geoPtr.get(), geo); }

/*!
 * \brief Overload for std::unique_ptr.
 */
template<class GeoType>
bool assign(std::unique_ptr<Geometry> geoPtr, GeoType& geo)
{ return assign(geoPtr.get(), geo); }

} // end namespace Frackit

#endif // FRACKIT_GEOMETRY_UTILITY_ASSIGN_HH
