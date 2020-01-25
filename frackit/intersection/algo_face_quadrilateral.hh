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
 *        between a face shape and a quadrilateral in 3d space.
 */
#ifndef FRACKIT_FACE_QUADRILATERAL_INTERSECTION_HH
#define FRACKIT_FACE_QUADRILATERAL_INTERSECTION_HH

#include "algo_face_planargeom.hh"
#include "intersectiontraits.hh"

namespace Frackit {
namespace IntersectionAlgorithms {

//! Intersect a face shape and a quadrilateral in 3d space.
//! The result can be composed of:
//! - points
//! - edges
//! - faces
//! Multiple of the above are possible
//! since the face shape might be curved.
template<class ctype>
Intersection< TopoDS_Face, Quadrilateral<ctype, 3> >
intersect_face_quadrilateral(const TopoDS_Face& face,
                             const Quadrilateral<ctype, 3>& quad,
                             ctype eps)
{ return intersect_face_planarGeometry(face, quad, eps); }

} // end namespace IntersectionAlgorithms
} // end namespace Frackit

#endif // FRACKIT_FACE_QUADRILATERAL_INTERSECTION_HH
