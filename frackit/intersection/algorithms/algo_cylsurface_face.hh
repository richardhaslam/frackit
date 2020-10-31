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
 *        between a lateral cylinder surface and
 *        a TopoDS_Face in 3d space.
 */
#ifndef FRACKIT_CYLINDERSURFACE_FACE_INTERSECTION_HH
#define FRACKIT_CYLINDERSURFACE_FACE_INTERSECTION_HH

#include <TopoDS_Face.hxx>
#include <frackit/geometry/cylindersurface.hh>
#include <frackit/occ/breputilities.hh>

#include <frackit/intersection/intersectiontraits.hh>
#include "algo_face_face_3d.hh"

namespace Frackit {
namespace IntersectionAlgorithms {

//! Intersect a lateral cylinder surface and a TopoDS_Face
//! The result can be:
//!     - points
//!     - edge shapes
//!     - face shapes
//!     Multiple of the above are possible since the faces shape might be curved.
template<class ctype>
Intersection< CylinderSurface<ctype>, TopoDS_Face >
intersect_cylinderSurface_face(const CylinderSurface<ctype>& cylSurface,
                               const TopoDS_Face& face,
                               ctype eps)
{ return intersect_face_face_3d(OCCUtilities::getShape(cylSurface), face, eps); }

} // end namespace IntersectionAlgorithms
} // end namespace Frackit

#endif // FRACKIT_CYLINDERSURFACE_FACE_INTERSECTION_HH
