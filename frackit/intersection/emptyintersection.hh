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
 * \ingroup Intersection
 * \brief Class used to define empty intersection.
 */
#ifndef FRACKIT_EMPTY_INTERSECTION_HH
#define FRACKIT_EMPTY_INTERSECTION_HH

#include <string>
#include <variant>

namespace Frackit {

/*!
 * \ingroup Intersection
 * \brief Class used to define empty intersection.
 * \tparam wd The dimension of the coordinate space.
 * \tparam CT The type used for coordinates
 */
template<int wd = 3, class CT = double>
struct EmptyIntersection
{
    //! Technically, this does not have a defined dimension
    static constexpr int myDimension() { return 0; }

    //! Return the dimension of the coordinate space
    static constexpr int worldDimension() { return wd; }

    //! Return the name of this geometry
    static std::string name() { return "EmptyIntersection"; }

    //! Export coordinate type for compatibility
    using ctype = CT;
};

/*!
 * \ingroup Intersection
 * \brief Traits class to identify empty intersections
 * \tparam IsGeometry The geometry type of an intersection
 */
template<class IsGeometry>
struct IsEmptyIntersection
{ static constexpr bool value = false; };

//! Specialization for the empty intersection class
template<int wd>
struct IsEmptyIntersection<EmptyIntersection<wd>>
{ static constexpr bool value = true; };

/*!
 * \ingroup Intersection
 * \brief Returns true if a geometry describes an empty intersection.
 * \tparam IsGeometry the geometry of an intersection
 */
template<class IsGeometry>
constexpr bool isEmptyIntersection(const IsGeometry& is)
{ return IsEmptyIntersection<IsGeometry>::value; }

/*!
 * \ingroup Intersection
 * \brief Overload for intersection variant
 */
template<class... T>
bool isEmptyIntersection(const std::variant<T...>& intersection)
{ return std::visit([&] (auto&& is) { return isEmptyIntersection(is); }, intersection); }

/*!
 * \ingroup Intersection
 * \brief Overload for general intersections possibly
 *        containing possibly various types.
 */
template<class Geo>
bool isEmptyIntersection(const std::vector<Geo>& intersections)
{
    return std::all_of(intersections.begin(),
                       intersections.end(),
                       [&] (const auto& is) { return isEmptyIntersection(is); } );
}

} // end namespace Frackit

#endif // FRACKIT_EMPTY_INTERSECTION_HH
