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
#include <algorithm>
#include <cassert>
#include <variant>
#include <vector>
#include <stdexcept>
#include <type_traits>
#include <limits>

#include <frackit/common/promotedtype.hh>
#include <frackit/geometry/point.hh>
#include <frackit/geometry/segment.hh>
#include <frackit/geometry/ellipse.hh>
#include <frackit/geometry/ellipsearc.hh>
#include <frackit/geometry/plane.hh>
#include <frackit/geometry/disk.hh>
#include <frackit/geometry/line.hh>
#include <frackit/geometry/cylindersurface.hh>

#include "emptyintersection.hh"
#include "intersectiontraits.hh"
#include "intersect.hh"

namespace Frackit {
namespace IntersectionPredicates {

/*!
 * \brief Returns true if a geometry
 *        describes an empty intersection.
 * \tparam IsGeometry the geometry of an intersection
 */
template<class IsGeometry,
         std::enable_if_t<!IsEmptyIntersection<IsGeometry>::value, int> = 0>
constexpr bool isEmpty(const IsGeometry& is)
{ return false; }

//! Overload for empty intersections
template<class IsGeometry,
         std::enable_if_t<IsEmptyIntersection<IsGeometry>::value, int> = 0>
constexpr bool isEmpty(const IsGeometry& is)
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
 * \note This overload is active when no specialization is available
 */
template<class Geo1, class Geo2, class IsGeometry>
double angle(const Geo1& geo1, const Geo2& geo2, const IsGeometry& isGeom)
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
template<class ctype1, int wd, class ctype2, class ctype3>
PromotedType<ctype1, ctype2>
angle(const Plane<ctype1, wd>& plane1,
      const Plane<ctype2, wd>& plane2,
      const Line<ctype3, wd>& isLine)
{
    using std::abs;
    using std::acos;
    using Vector = Vector<PromotedType<ctype1, ctype2>, wd>;
    return acos( abs(Vector(plane1.normal())*Vector(plane2.normal())) );
}

/*!
 * \brief Returns the angle in which two planes intersect
 * \note This overload is used when one knows that the planes
 *       intersect but does not want to compute the intersection.
 */
template<class ctype1, int wd, class ctype2>
PromotedType<ctype1, ctype2>
angle(const Plane<ctype1, wd>& plane1,
      const Plane<ctype2, wd>& plane2)
{
    assert( !isEmpty(intersect(plane1, plane2)) );

    using std::abs;
    using std::acos;
    using Vector = Vector<PromotedType<ctype1, ctype2>, wd>;
    return acos( abs(Vector(plane1.normal())*Vector(plane2.normal())) );
}

/*!
 * \brief Returns the angle in which two disks intersect in a point
 */
template<class ctype1, class ctype2, class ctype3>
PromotedType<ctype1, ctype2>
angle(const Disk<ctype1>& disk1,
      const Disk<ctype2>& disk2,
      const Point<ctype3, 3>& isPoint)
{ return angle(disk1.supportingPlane(), disk2.supportingPlane()); }

/*!
 * \brief Returns the angle in which two disks intersect in a segment
 */
template<class ctype1, class ctype2, class ctype3>
PromotedType<ctype1, ctype2>
angle(const Disk<ctype1>& disk1,
      const Disk<ctype2>& disk2,
      const Segment<ctype3, 3>& isPoint)
{ return angle(disk1.supportingPlane(), disk2.supportingPlane()); }

/*!
 * \brief Returns the angle in which a disk and a cylinder surface intersect
 * \note This overload is for the two intersecting in a point
 */
template<class ctype1, class ctype2, class ctype3, int wd>
PromotedType<ctype1, ctype2>
angle(const Disk<ctype1>& disk,
      const CylinderSurface<ctype2>& cylSurface,
      const Point<ctype3, wd>& isPoint)
{ return angle(disk.supportingPlane(), cylSurface.getTangentPlane(isPoint)); }

/*!
 * \brief Overload with swapped arguments
 */
template<class ctype1, class ctype2, class ctype3, int wd>
PromotedType<ctype1, ctype2>
angle(const CylinderSurface<ctype1>& cylSurface,
      const Disk<ctype2>& disk,
      const Point<ctype3, wd>& isPoint)
{ return angle(disk, cylSurface, isPoint); }

/*!
 * \brief Returns the angle in which a disk and a cylinder surface intersect
 * \note This overload is for the two intersecting in a segment
 */
template<class ctype1, class ctype2, class ctype3, int wd>
PromotedType<ctype1, ctype2>
angle(const Disk<ctype1>& disk,
      const CylinderSurface<ctype2>& cylSurface,
      const Segment<ctype3, wd>& isSeg)
{ return angle(disk.supportingPlane(), cylSurface.getTangentPlane(isSeg.source())); }

/*!
 * \brief Overload with swapped arguments
 */
template<class ctype1, class ctype2, class ctype3, int wd>
PromotedType<ctype1, ctype2>
angle(const CylinderSurface<ctype1>& cylSurface,
      const Disk<ctype2>& disk,
      const Segment<ctype3, wd>& isSeg)
{ return angle(disk, cylSurface, isSeg); }

/*!
 * \brief Returns the angle in which a disk and a cylinder surface intersect
 * \note This overload is for the two intersecting in an ellipse arc
 */
template<class ctype1, class ctype2, class ctype3, int wd>
PromotedType<ctype1, ctype2>
angle(const Disk<ctype1>& disk,
      const CylinderSurface<ctype2>& cylSurface,
      const EllipseArc<ctype3, wd>& isArc)
{
    // use the minimum angle between the disk plane and the
    // tangent plane on the surface at four sample points
    using std::min;
    using ctype = PromotedType<ctype1, ctype2>;
    std::vector<ctype> params({0.0, 0.25, 0.5, 0.75, 1.0});

    const auto& diskPlane = disk.supportingPlane();
    ctype resultAngle = std::numeric_limits<ctype>::max();
    for (auto param : params)
    {
        const auto p = isArc.getPoint(param);
        const auto tangentPlane = cylSurface.getTangentPlane(p);
        resultAngle = min(resultAngle, angle(diskPlane, tangentPlane));
    }

    return resultAngle;
}

/*!
 * \brief Overload with swapped arguments
 */
template<class ctype1, class ctype2, class ctype3, int wd>
PromotedType<ctype1, ctype2>
angle(const CylinderSurface<ctype1>& cylSurface,
      const Disk<ctype2>& disk,
      const EllipseArc<ctype3, wd>& isArc)
{ return angle(disk, cylSurface, isArc); }

/*!
 * \brief Returns the angle in which a disk and a cylinder surface intersect
 * \note This overload is for the two intersecting in an ellipse
 */
template<class ctype1, class ctype2, class ctype3, int wd>
PromotedType<ctype1, ctype2>
angle(const Disk<ctype1>& disk,
      const CylinderSurface<ctype2>& cylSurface,
      const Ellipse<ctype3, wd>& isEllipse)
{
    // use the minimum angle between the disk plane and the
    // tangent plane on the surface at eight sample points
    using std::min;
    using ctype = PromotedType<ctype1, ctype2>;
    std::vector<ctype> params({0.0, 0.125, 0.25, 0.375, 0.5, 0.625,  0.75, 0.875, 1.0});

    const auto& diskPlane = disk.supportingPlane();
    ctype resultAngle = std::numeric_limits<ctype>::max();
    for (auto param : params)
    {
        const auto p = isEllipse.getPoint(param);
        const auto tangentPlane = cylSurface.getTangentPlane(p);
        resultAngle = min(resultAngle, angle(diskPlane, tangentPlane));
    }

    return resultAngle;
}

/*!
 * \brief Overload with swapped arguments
 */
template<class ctype1, class ctype2, class ctype3, int wd>
PromotedType<ctype1, ctype2>
angle(const CylinderSurface<ctype1>& cylSurface,
      const Disk<ctype2>& disk,
      const Ellipse<ctype3, wd>& isEllipse)
{ return angle(disk, cylSurface, isEllipse); }

//! Overload for intersection variant
template<class Geo1, class Geo2, class... T>
PromotedType<typename Geo1::ctype, typename Geo2::ctype>
angle(const Geo1& geo1,
      const Geo2& geo2,
      const std::variant<T...>& intersection)
{ return std::visit([&] (auto&& is) { return angle(geo1, geo2, is); }, intersection); }

//! Overload for general intersections containing possibly various types
template<class Geo1, class Geo2, class T>
PromotedType<typename Geo1::ctype, typename Geo2::ctype>
angle(const Geo1& geo1,
      const Geo2& geo2,
      const std::vector<T>& intersections)
{
    using std::min;
    using ctype = PromotedType<typename Geo1::ctype, typename Geo2::ctype>;

    ctype result = std::numeric_limits<ctype>::max();
    for (const auto& is : intersections)
        result = min(result, angle(geo1, geo2, is));

    return result;
}

} // end namespace IntersectionPredicates
} // end namespace Frackit

#endif // FRACKIT_INTERSECTION_PREDICATES_HH
