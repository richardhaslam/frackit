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
 * \brief Contains the intersection algorithm
 *        between two planes.
 */
#ifndef FRACKIT_PLANE_PLANE_INTERSECTION_HH
#define FRACKIT_PLANE_PLANE_INTERSECTION_HH

#include <stdexcept>

#include <gp_Vec.hxx>
#include <Geom_Plane.hxx>
#include <GeomAPI_IntSS.hxx>
#include <Standard_Handle.hxx>

#include <frackit/geometry/plane.hh>
#include <frackit/precision/precision.hh>
#include <frackit/occ/gputilities.hh>

#include "intersectiontraits.hh"
#include "emptyintersection.hh"

namespace Frackit {
namespace IntersectionAlgorithms {

//! Intersect two planes
//! The result can be:
//! - a plane
//! - a line
//! - no intersection
template<class ctype>
Intersection< Plane<ctype, 3>, Plane<ctype, 3> >
intersect_plane_plane(const Plane<ctype, 3>& plane1,
                      const Plane<ctype, 3>& plane2,
                      ctype eps)
{
    using ResultType = Intersection< Plane<ctype, 3>, Plane<ctype, 3> >;

    const auto d1 = OCCUtilities::direction(plane1.normal());
    const auto d2 = OCCUtilities::direction(plane2.normal());

    // if planes are identical, the result is a plane or empty
    if (gp_Vec(d1).IsParallel(gp_Vec(d2), Precision<ctype>::angular()))
    {
        if (plane1.contains(plane2.supportingPoint(), eps))
            return ResultType( plane1 );
        else
            return ResultType( EmptyIntersection<3>() );
    }

    // find the intersection line between the two planes
    const auto c1 = OCCUtilities::point(plane1.supportingPoint());
    const auto c2 = OCCUtilities::point(plane2.supportingPoint());

    Handle(Geom_Plane) gpPlane1 = new Geom_Plane(c1, d1);
    Handle(Geom_Plane) gpPlane2 = new Geom_Plane(c2, d2);
    GeomAPI_IntSS interSection(gpPlane1, gpPlane2, eps);
    if (!interSection.IsDone())
        throw std::runtime_error("Could not perform plane-plane intersection");

    assert(interSection.NbLines() == 1);
    assert(!interSection.Line(1)->IsClosed());
    assert(interSection.Line(1)->IsCN(1));

    // The result is an infinite line.
    // The parameter range is -inf -> +inf
    // We evaluate a support point and direction for param = 0.
    const auto& sp = interSection.Line(1)->Value(0.0);
    const auto& dir = interSection.Line(1)->DN(0.0, /*derivative order*/1);

    // create Line object from support point and direction
    using Line = Frackit::Line<ctype, 3>;
    return ResultType( Line(OCCUtilities::point(sp),
                            OCCUtilities::direction(dir)) );
}

} // end namespace IntersectionAlgorithms
} // end namespace Frackit

#endif // FRACKIT_PLANE_PLANE_INTERSECTION_HH
