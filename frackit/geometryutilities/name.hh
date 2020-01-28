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
 * \ingroup GeometryUtilities
 * \brief Utility functionality to get the name of geometries.
 *        This provides an interface that is compatible also
 *        with classes from the BRep package.
 */
#ifndef FRACKIT_GEOMETRY_NAME_HH
#define FRACKIT_GEOMETRY_NAME_HH

#include <stdexcept>
#include <string>

#include <TopoDS_Solid.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Shape.hxx>

namespace Frackit {

/*!
 * \ingroup GeometryUtilities
 * \brief Return the name of an internal geometry class.
 * \note Here, we simply forward to the name() function of the geometries
 */
template<class Geo>
std::string geometryName(const Geo& geo)
{ return geo.name(); }

/*!
 * \ingroup GeometryUtilities
 * \brief Return the name for a TopoDS_Solid.
 */
std::string geometryName(const TopoDS_Solid& s)
{ return "TopoDS_Solid"; }

/*!
 * \ingroup GeometryUtilities
 * \brief Return the name for a TopoDS_Shell.
 */
std::string geometryName(const TopoDS_Shell& s)
{ return "TopoDS_Shell"; }

/*!
 * \ingroup GeometryUtilities
 * \brief Return the name for a TopoDS_Face.
 */
std::string geometryName(const TopoDS_Face& f)
{ return "TopoDS_Face"; }

/*!
 * \ingroup GeometryUtilities
 * \brief Return the name for a TopoDS_Wire.
 */
std::string geometryName(const TopoDS_Wire& w)
{ return "TopoDS_Wire"; }

/*!
 * \ingroup GeometryUtilities
 * \brief Return the name for a TopoDS_Edge.
 */
std::string geometryName(const TopoDS_Edge& e)
{ return "TopoDS_Edge"; }

/*!
 * \ingroup GeometryUtilities
 * \brief Return the name for a TopoDS_Vertex.
 */
std::string geometryName(const TopoDS_Vertex& v)
{ return "TopoDS_Vertex"; }

/*!
 * \ingroup GeometryUtilities
 * \brief Return the name for a TopoDS_Shape.
 */
std::string geometryName(const TopoDS_Shape& v)
{ return "TopoDS_Shape"; }

} // end namespace Frackit

#endif // FRACKIT_GEOMETRY_NAME_HH
