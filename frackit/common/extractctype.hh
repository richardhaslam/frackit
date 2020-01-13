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
 * \brief Type traits to extract the coordinate type of geometric objects.
 */
#ifndef FRACKIT_COMMON_COORDINATE_TYPE_HH
#define FRACKIT_COMMON_COORDINATE_TYPE_HH

#include <utility>

#include <TopoDS_Shape.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Solid.hxx>
#include <Standard.hxx>

namespace Frackit {

/*!
* \brief Traits class to extract the coordinate type from a geometry.
*/
template<class Geom>
struct CoordinateTypeTraits
{ using type = typename Geom::ctype; };

/*!
* \brief Specialization for Brep shapes.
*/
template<>
struct CoordinateTypeTraits<TopoDS_Shape>
{ using type = Standard_Real; };

/*!
 * \brief Specialization for Brep vertices.
 */
template<>
struct CoordinateTypeTraits<TopoDS_Vertex>
: public CoordinateTypeTraits<TopoDS_Shape>
{};

/*!
 * \brief Specialization for Brep edges.
 */
template<>
struct CoordinateTypeTraits<TopoDS_Edge>
: public CoordinateTypeTraits<TopoDS_Shape>
{};

/*!
 * \brief Specialization for Brep faces.
 */
template<>
struct CoordinateTypeTraits<TopoDS_Face>
: public CoordinateTypeTraits<TopoDS_Shape>
{};

/*!
 * \brief Specialization for Brep solids.
 */
template<>
struct CoordinateTypeTraits<TopoDS_Solid>
: public CoordinateTypeTraits<TopoDS_Shape>
{};

} // end namespace Frackit

#endif // FRACKIT_COMMON_COORDINATE_TYPE_HH
