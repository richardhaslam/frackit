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
#ifndef FRACKIT_INTERSECTION_TRAITS_HH
#define FRACKIT_INTERSECTION_TRAITS_HH

#include <type_traits>
#include <variant>
#include <vector>

#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shell.hxx>

#include <frackit/geometry/point.hh>
#include <frackit/geometry/segment.hh>
#include <frackit/geometry/plane.hh>
#include <frackit/geometry/disk.hh>
#include <frackit/geometry/ellipse.hh>
#include <frackit/geometry/ellipsearc.hh>
#include <frackit/geometry/cylindersurface.hh>

#include "emptyintersection.hh"

namespace Frackit {

/*!
 * \brief Traits class to identify empty intersections
 * \tparam IsGeometry The geometry type of an intersection
 */
template<class IsGeometry>
struct IsEmptyIntersection
{ static constexpr bool value = false; };

//! Specialization for empty intersections
template<int wd>
struct IsEmptyIntersection<EmptyIntersection<wd>>
{ static constexpr bool value = true; };

//! \todo TODO Doc me.
template<class Geometry1, class Geometry2>
struct IntersectionTraits;

//! Convenience alias to obtain intersection type
template<class Geometry1, class Geometry2>
using Intersection = typename IntersectionTraits<Geometry1, Geometry2>::type;

//! \todo TODO Doc me.
template<class ctype, int wd>
struct IntersectionTraits< Segment<ctype, wd>, Segment<ctype, wd> >
{
    using type = std::variant< Segment<ctype, wd>,
                               Point<ctype, wd>,
                               EmptyIntersection<wd> >;
};

//! \todo TODO Doc me.
template<class ctype>
struct IntersectionTraits< Plane<ctype, 3>, Plane<ctype, 3> >
{
    using type = std::variant< Line<ctype, 3>,
                               Plane<ctype, 3>,
                               EmptyIntersection<3> >;
};

//! \todo TODO Doc me.
template<class ctype>
struct IntersectionTraits< Plane<ctype, 3>, Line<ctype, 3> >
{
    using type = std::variant< Point<ctype, 3>,
                               Line<ctype, 3>,
                               EmptyIntersection<3> >;
};

//! \todo TODO Doc me.
template<class ctype>
struct IntersectionTraits< Line<ctype, 3>, Plane<ctype, 3> >
: public IntersectionTraits< Plane<ctype, 3>, Line<ctype, 3> >
{};

//! \todo TODO Doc me.
template<class ctype>
struct IntersectionTraits< Disk<ctype>, Line<ctype, 3> >
{
    static constexpr int wd = Disk<ctype>::worldDimension();
    static_assert(wd == 3, "World dimension of 3 expected");

    using type = std::variant< Point<ctype, wd>,
                               Segment<ctype, wd>,
                               EmptyIntersection<wd> >;
};

//! \todo TODO Doc me.
template<class ctype>
struct IntersectionTraits< Line<ctype, 3>, Disk<ctype> >
: public IntersectionTraits< Disk<ctype>, Line<ctype, 3> >
{};

//! \todo TODO Doc me.
template<class ctype>
struct IntersectionTraits< Disk<ctype>, Disk<ctype> >
{
    static constexpr int wd = Disk<ctype>::worldDimension();
    static_assert(wd == 3, "World dimension of 3 expected");

    using type = std::variant< Point<ctype, wd>,
                               Segment<ctype, wd>,
                               Disk<ctype>,
                               EmptyIntersection<wd> >;
};

//! \todo TODO Doc me.
template<class ctype>
struct IntersectionTraits< CylinderSurface<ctype>, Disk<ctype> >
{
    static constexpr int wd = Disk<ctype>::worldDimension();
    static constexpr int wd2 = CylinderSurface<ctype>::worldDimension();
    static_assert(wd == 3, "World dimension of 3 expected");
    static_assert(wd == wd2, "World dimension of the geometries must match");

    using BaseType = std::variant< Point<ctype, wd>,
                                   Segment<ctype, wd>,
                                   EllipseArc<ctype, wd>,
                                   Ellipse<ctype, wd>,
                                   EmptyIntersection<wd> >;
    using type = std::vector<BaseType>;
};

//! \todo TODO Doc me.
template<class ctype>
struct IntersectionTraits< Disk<ctype>, CylinderSurface<ctype> >
: public IntersectionTraits< CylinderSurface<ctype>, Disk<ctype> >
{};

//! \todo TODO Doc me.
template<class ctype>
struct IntersectionTraits< TopoDS_Shell, Disk<ctype> >
{
    using BaseType = std::variant< Point<ctype, 3>,
                                   TopoDS_Edge,
                                   TopoDS_Face,
                                   EmptyIntersection<3> >;
    using type = std::vector<BaseType>;
};

//! \todo TODO Doc me.
template<class ctype>
struct IntersectionTraits< Disk<ctype>, TopoDS_Shell >
: public IntersectionTraits< TopoDS_Shell, Disk<ctype> >
{};

//! \todo TODO Doc me.
template<class ctype>
struct IntersectionTraits< TopoDS_Face, Disk<ctype> >
{
    using BaseType = std::variant< Point<ctype, 3>,
                                   TopoDS_Edge,
                                   TopoDS_Face,
                                   EmptyIntersection<3> >;
    using type = std::vector<BaseType>;
};

//! \todo TODO Doc me.
template<class ctype>
struct IntersectionTraits< Disk<ctype>, TopoDS_Face >
: public IntersectionTraits< TopoDS_Face, Disk<ctype> >
{};

} // end namespace Frackit

#endif // FRACKIT_INTERSECTION_TRAITS_HH
