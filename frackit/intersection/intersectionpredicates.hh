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
 * \brief \todo TODO doc me.
 */
#ifndef FRACKIT_INTERSECTION_PREDICATES_HH
#define FRACKIT_INTERSECTION_PREDICATES_HH

#include <cmath>
#include <cassert>
#include <variant>
#include <stdexcept>

#include <frackit/common/promotedtype.hh>
#include <frackit/geometry/plane.hh>
#include <frackit/geometry/disk.hh>

#include "emptyintersection.hh"
#include "intersect.hh"

namespace Frackit {
namespace IntersectionPredicates {

/*!
 * \brief Returns true if a geometry
 *        describes an empty intersection.
 * \tparam IsGeometry the geometry of an intersection
 */
template<class IsGeometry>
constexpr bool isEmpty(const IsGeometry& is)
{ return false; }

//! Overload for empty intersections
template<int wd>
constexpr bool isEmpty(const EmptyIntersection<wd>& empty)
{ return true; }

//! Overload for intersection variant
template<class... T>
bool isEmpty(const std::variant<T...>& intersection)
{ return std::visit([&] (auto&& is) { return isEmpty(is); }, intersection); }

//! Overload for general intersections containing possibly various types
template<class Geo>
bool isEmpty(const std::vector<Geo>& intersections)
{
    return std::all_of(intersections.begin(),
                       intersections.end(),
                       [&] (const auto& is) { return isEmpty(is); } );
}

/*!
 * \brief Returns the angle in which two geometries intersect
 */
template<class Geo1, class Geo2>
double angle(const Geo1& geo1, const Geo2& geo2)
{
    std::string msg = "IntersectionPredicates::angle() not implemented for ";
    msg += "\"" + Geo1::name() + "\"";
    msg += " and ";
    msg += "\"" + Geo2::name() + "\"";
    throw std::runtime_error( msg );
}

/*!
 * \brief Returns the angle in which two planes intersect
 */
template<class ctype1, int wd, class ctype2>
PromotedType<ctype1, ctype2>
angle(const Plane<ctype1, wd>& plane1, const Plane<ctype2, wd>& plane2)
{
    assert( !isEmpty(intersect(plane1, plane2)) );

    using std::abs;
    using std::acos;
    using Vector = Vector<PromotedType<ctype1, ctype2>, wd>;
    return acos( abs(Vector(plane1.normal())*Vector(plane2.normal())) );
}

/*!
 * \brief Returns the angle in which two disks intersect
 */
template<class ctype1, class ctype2>
PromotedType<ctype1, ctype2>
angle(const Disk<ctype1>& disk1, const Disk<ctype2>& disk2)
{
    assert( !isEmpty(intersect(disk1, disk2)) );
    return angle(disk1.supportingPlane(), disk2.supportingPlane());
}

} // end namespace IntersectionPredicates
} // end namespace Frackit

#endif // FRACKIT_INTERSECTION_PREDICATES_HH
