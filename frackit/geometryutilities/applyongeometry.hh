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
 * \brief Utility functionality to apply functions on
 *        geometrical objects of which only the pointer
 *        to the virtual base class is available. In a
 *        first step, the geometry is cast into ist actual
 *        implementation and then the provided lambda function
 *        is evaluated on the actual geometry.
 */
#ifndef FRACKIT_GEOMETRY_UTILITY_APPLY_ON_GEOMETRY_HH
#define FRACKIT_GEOMETRY_UTILITY_APPLY_ON_GEOMETRY_HH

#include <memory>
#include <string>
#include <stdexcept>
#include <type_traits>

#include <frackit/geometry/geometry.hh>
#include <frackit/geometryutilities/name.hh>

// The supported geometry types.
// If a new geometry is added, it has to be included here and the
// implementation in the applyOnGeometry() function has to be extended.
#include <frackit/geometry/disk.hh>

#include "assign.hh"

namespace Frackit {

/*!
 * \brief Apply a function to the geometry referenced
 *        to by the given pointer on the geometry interface.
 * \param geoPtr Pointer to an object of the geometry interface
 * \param applyFunc The function to be executed on the geometry
 * \note All possible geometry types must be hardcoded here
 */
template<class ApplyFunc, class ctype = double>
auto applyOnGeometry(ApplyFunc&& applyFunc, Geometry* geoPtr)
{
    if (geoPtr->name() == "Disk")
    {
        Disk<ctype> disk;
        if (!assign(geoPtr, disk))
            throw std::runtime_error("applyOnGeometry(): could not assign disk");
        return applyFunc(disk);
    }

    std::string msg = "applyOnGeometry() function not implemented ";
    msg += "for geometry with name \"";
    msg += geoPtr->name();
    msg += "\"\n";
    throw std::runtime_error(msg);
}

/*!
 * \brief Overload for std::shared_ptr.
 */
template<class ApplyFunc, class ctype = double>
auto applyOnGeometry(ApplyFunc&& applyFunc, std::shared_ptr<Geometry> geoPtr)
{ return applyOnGeometry<ApplyFunc, ctype>(applyFunc, geoPtr.get()); }

/*!
 * \brief Overload for std::unique_ptr.
 */
template<class ApplyFunc, class ctype = double>
auto applyOnGeometry(std::unique_ptr<Geometry> geoPtr, ApplyFunc& applyFunc)
{ return applyOnGeometry<ApplyFunc, ctype>(applyFunc, geoPtr.get()); }

} // end namespace Frackit

#endif // FRACKIT_GEOMETRY_UTILITY_APPLY_ON_GEOMETRY_HH
