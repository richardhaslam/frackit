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
 * \brief Implementation details on intersection dimension constraint enforcements.
 */
#ifndef FRACKIT_ENTITYNETWORK_CONSTRAINT_IMPL_ADMISSIBLE_DIMENSION_HH
#define FRACKIT_ENTITYNETWORK_CONSTRAINT_IMPL_ADMISSIBLE_DIMENSION_HH

#include <variant>
#include <vector>

#include <TopoDS_Vertex.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Solid.hxx>

#include <frackit/common/extractdimension.hh>

namespace Frackit {
namespace ConstraintImpl {

    /*!
     * \brief Evaluates if the dimension of a geometry less or equal as maxDim.
     */
    template<class Geo>
    bool isAdmissibleDimension(const Geo& isGeometry, int maxDim)
    { return getDimension(isGeometry) <= maxDim; }

    /*!
    * \brief Overload for an std::variant
    */
    template<class... T>
    bool isAdmissibleDimension(const std::variant<T...>& intersection, int maxDim)
    { return std::visit([&] (auto&& is) { return isAdmissibleDimension(is, maxDim); }, intersection); }

    /*!
     * \brief Overload for a vector of intersections
     */
    template<class T>
    bool isAdmissibleDimension(const std::vector<T>& intersection, int maxDim)
    {
        return std::all_of(intersection.begin(),
                           intersection.end(),
                           [&] (const auto& is) { return isAdmissibleDimension(is, maxDim); });
    }

} // end namespace ConstraintImpl
} // end namespace Frackit

#endif // FRACKIT_ENTITYNETWORK_CONSTRAINT_IMPL_ADMISSIBLE_DIMENSION_HH
