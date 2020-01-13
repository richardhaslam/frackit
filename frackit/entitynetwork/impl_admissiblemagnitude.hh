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
 * \brief Implementation details on magnitude constraint enforcements.
 */
#ifndef FRACKIT_ENTITYNETWORK_CONSTRAINT_IMPL_ADMISSIBLE_MAGNITUDE_HH
#define FRACKIT_ENTITYNETWORK_CONSTRAINT_IMPL_ADMISSIBLE_MAGNITUDE_HH

#include <variant>
#include <vector>

#include <TopoDS_Vertex.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Solid.hxx>

#include <frackit/magnitude/magnitude.hh>

namespace Frackit {
namespace ConstraintImpl {

    /*!
     * \brief Evaluates if the magnitude of a geometry is above a given threshold.
     */
    template<class Geo, class ctype>
    bool isAdmissibleMagnitude(const Geo& isGeometry, ctype threshold)
    {
        // zero dimensional (point) intersections fulfill any magnitude constraint
        if (Geo::myDimension() == 0) return true;
        // empty intersections fulfill any magnitude constraint
        if (IsEmptyIntersection<Geo>::value) return true;
        // any other intersection geometry
        return computeMagnitude(isGeometry) >= threshold;
    }

    /*!
     * \brief Evaluates if the magnitude of a vertex shape is above a given threshold.
     */
    template<class ctype>
    bool isAdmissibleMagnitude(const TopoDS_Vertex& isVertex, ctype threshold)
    { return true; }

    /*!
     * \brief Evaluates if the magnitude of an edge is above a given threshold.
     */
    template<class ctype>
    bool isAdmissibleMagnitude(const TopoDS_Edge& isEdge, ctype threshold)
    { return computeMagnitude(isEdge) >= threshold; }

    /*!
     * \brief Evaluates if the magnitude of a face is above a given threshold.
     */
    template<class ctype>
    bool isAdmissibleMagnitude(const TopoDS_Face& isFace, ctype threshold)
    { return computeMagnitude(isFace) >= threshold; }

    /*!
     * \brief Evaluates if the magnitude of a solid is above a given threshold.
     */
    template<class ctype>
    bool isAdmissibleMagnitude(const TopoDS_Solid& isSolid, ctype threshold)
    { return computeMagnitude(isSolid) >= threshold; }

    /*!
    * \brief Overload for an std::variant
    */
    template<class... T, class ctype>
    bool isAdmissibleMagnitude(const std::variant<T...>& intersection, ctype threshold)
    { return std::visit([&] (auto&& is) { return isAdmissibleMagnitude(is, threshold); }, intersection); }

    /*!
     * \brief Overload for a vector of intersections
     */
    template<class T, class ctype>
    bool isAdmissibleMagnitude(const std::vector<T>& intersection, ctype threshold)
    {
        return std::all_of(intersection.begin(),
                           intersection.end(),
                           [&] (const auto& is) { return isAdmissibleMagnitude(is, threshold); });
    }

} // end namespace ConstraintImpl
} // end namespace Frackit

#endif // FRACKIT_ENTITYNETWORK_CONSTRAINT_IMPL_ADMISSIBLE_MAGNITUDE_HH
