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
 * \brief Implementation details on distance to boundary constraints.
 */
#ifndef FRACKIT_ENTITYNETWORK_CONSTRAINT_IMPL_ADMISSIBLE_BOUNDARY_DISTANCE_HH
#define FRACKIT_ENTITYNETWORK_CONSTRAINT_IMPL_ADMISSIBLE_BOUNDARY_DISTANCE_HH

#include <stdexcept>
#include <variant>
#include <vector>

#include <frackit/distance/distance.hh>
#include <frackit/distance/distancetoboundary.hh>
#include <frackit/distance/pointonboundary.hh>

#include <frackit/magnitude/magnitude.hh>
#include <frackit/intersection/intersectiontraits.hh>

#include <frackit/geometry/point.hh>
#include <frackit/geometry/segment.hh>

namespace Frackit {
namespace ConstraintImpl {

    /*!
     * \brief Evaluates if the distance of an intersection geometry to
     *        the boundary of an entity is above a given threshold.
     */
    template<class IsGeom, class Geo, class ctype>
    bool isAdmissibleDistanceToBoundary(const IsGeom& is, const Geo& entity, ctype threshold)
    {
        std::string msg = "Distance to boundary not implemented for intersection geometry ";
        msg += "\"" + IsGeom::name() + "\"";
        msg += " and entity geometry ";
        msg += "\"" + Geo::name() + "\"";
        throw std::runtime_error(std::string(msg));
    }

    /*!
     * \brief Overload for empty intersections. These fulfill the distance constraint.
     */
    template<int wd, class Geo, class ctype>
    bool isAdmissibleDistanceToBoundary(const EmptyIntersection<wd>& is,
                                        const Geo& entity,
                                        ctype threshold)
    { return true; }

    /*!
     * \brief Overload for point intersections.
     *        Points that lie on the boundary of the entity are considered
     *        to fulfill the distance constraint.
     */
    template<class ctype, int wd, class Geo, class ctype2>
    bool isAdmissibleDistanceToBoundary(const Point<ctype, wd>& is,
                                        const Geo& entity,
                                        ctype2 threshold)
    {
        if (pointOnGeometryBoundary(is, entity))
            return true;
        return computeDistanceToBoundary(is, entity) >= threshold;
    }

    /*!
     * \brief Overload for point intersections.
     *        Points that lie on the boundary of the entity are considered
     *        to fulfill the distance constraint.
     */
    template<class ctype, int wd, class Geo, class ctype2>
    bool isAdmissibleDistanceToBoundary(const Segment<ctype, wd>& seg,
                                        const Geo& entity,
                                        ctype2 threshold)
    {
        return isAdmissibleDistanceToBoundary(seg.source(), entity, threshold)
               && isAdmissibleDistanceToBoundary(seg.target(), entity, threshold);
    }

    /*!
     * \brief Overload for std::variant.
     */
    template<class... T, class Geo, class ctype>
    bool isAdmissibleDistanceToBoundary(const std::variant<T...>& intersection,
                                        const Geo& entity,
                                        ctype threshold)
    {
        return std::visit([&] (auto&& is)
                          { return isAdmissibleDistanceToBoundary_(is, entity, threshold); },
                          intersection);
    }

    /*!
     * \brief Overload for a vector of intersections.
     */
    template<class T, class Geo, class ctype>
    bool isAdmissibleDistanceToBoundary(const std::vector<T>& intersection,
                                        const Geo& entity,
                                        ctype threshold)
    {
        return std::all_of(intersection.begin(),
                           intersection.end(),
                           [&] (const auto& is)
                           { return isAdmissibleDistanceToBoundary(is, entity, threshold); });
    }

} // end namespace ConstraintImpl
} // end namespace Frackit

#endif // FRACKIT_ENTITYNETWORK_CONSTRAINT_IMPL_ADMISSIBLE_BOUNDARY_DISTANCE_HH
