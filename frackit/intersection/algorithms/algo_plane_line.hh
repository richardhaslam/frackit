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
 * \brief Contains the intersection algorithm
 *        between a plane and a line.
 */
#ifndef FRACKIT_PLANE_LINE_INTERSECTION_HH
#define FRACKIT_PLANE_LINE_INTERSECTION_HH

#include <stdexcept>

#include <Geom_Line.hxx>
#include <Geom_Surface.hxx>
#include <GeomAPI_IntCS.hxx>
#include <Standard_Handle.hxx>

#include <frackit/geometry/plane.hh>
#include <frackit/geometry/line.hh>

#include <frackit/precision/precision.hh>
#include <frackit/occ/breputilities.hh>
#include <frackit/occ/gputilities.hh>

#include <frackit/intersection/intersectiontraits.hh>
#include <frackit/intersection/emptyintersection.hh>

namespace Frackit {
namespace IntersectionAlgorithms {

//! Intersect a plane and a line
//! The result can be:
//! - a line
//! - a point
//! - no intersection
template<class ctype>
Intersection< Plane<ctype, 3>, Line<ctype, 3> >
intersect_plane_line(const Plane<ctype, 3>& plane,
                     const Line<ctype, 3>& line,
                     ctype eps)
{
    using ResultType = Intersection< Plane<ctype, 3>, Line<ctype, 3> >;

    // if the line is not parallel to the plane, the result can only be a point
    const auto lineDir = OCCUtilities::direction(line.direction());
    const auto planeNormal = OCCUtilities::direction(plane.normal());
    if (!lineDir.IsNormal(planeNormal, Precision<ctype>::angular()))
    {
        const auto linePoint = OCCUtilities::point(line.supportingPoint());
        const auto planePoint = OCCUtilities::point(plane.supportingPoint());

        // let the the geometry package compute the intersection
        Handle(Geom_Surface) planeHandle = new Geom_Plane(planePoint, planeNormal);
        Handle(Geom_Line) lineHandle = new Geom_Line(linePoint, lineDir);
        GeomAPI_IntCS interSection(lineHandle, planeHandle);
        if (!interSection.IsDone())
            throw std::runtime_error("Could not perform disk-line intersection");

        assert(interSection.NbSegments() == 0);
        assert(interSection.NbPoints() == 1);

        // check if point is on the disk
        const auto p = OCCUtilities::point(interSection.Point(1));
        if (plane.contains(p, eps))
            return ResultType( p );
        else
            return ResultType( EmptyIntersection<3>() );
    }

    // The line is parallel. If the distance is > eps, there is no intersection
    const auto d = Vector<ctype, 3>(line.supportingPoint(),
                                    plane.projection(line.supportingPoint()));
    if (d.length() > eps)
        return ResultType( EmptyIntersection<3>() );

    // the intersection is the line
    return ResultType( line );
}

} // end namespace IntersectionAlgorithms
} // end namespace Frackit

#endif // FRACKIT_PLANE_LINE_INTERSECTION_HH
