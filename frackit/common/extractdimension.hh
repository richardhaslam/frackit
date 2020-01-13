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
 * \brief Type traits to extract the dimensionality of geometric objects.
 */
#ifndef FRACKIT_COMMON_EXTRACT_DIMENSION_HH
#define FRACKIT_COMMON_EXTRACT_DIMENSION_HH

#include <utility>

#include <TopoDS_Shape.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>
#include <Standard.hxx>

namespace Frackit {

/*!
 * \brief Traits class to extract the dimensionalty of a geometry.
 * \note This is the specialization used for internal geometry types
 *       which carry the information of the dimensions themselves.
 */
template<class Geom>
struct DimensionalityTraits
{
    static constexpr int geomDim = Geom::myDimension();
    static constexpr int worldDim = Geom::worldDimension();
};

/*!
 * \brief Specialization for Brep shapes.
 */
template<>
struct DimensionalityTraits<TopoDS_Shape>
{
    static constexpr int worldDim = 3;
};

/*!
 * \brief Specialization for Brep vertices.
 */
template<>
struct DimensionalityTraits<TopoDS_Vertex>
: public DimensionalityTraits<TopoDS_Shape>
{
    static constexpr int geomDim = 0;
};

/*!
 * \brief Specialization for Brep edges.
 */
template<>
struct DimensionalityTraits<TopoDS_Edge>
: public DimensionalityTraits<TopoDS_Shape>
{
    static constexpr int geomDim = 1;
};

/*!
 * \brief Specialization for Brep wires.
 */
template<>
struct DimensionalityTraits<TopoDS_Wire>
: public DimensionalityTraits<TopoDS_Shape>
{
    static constexpr int geomDim = 1;
};

/*!
 * \brief Specialization for Brep faces.
 */
template<>
struct DimensionalityTraits<TopoDS_Face>
: public DimensionalityTraits<TopoDS_Shape>
{
    static constexpr int geomDim = 2;
};

/*!
 * \brief Specialization for Brep shells.
 */
template<>
struct DimensionalityTraits<TopoDS_Shell>
: public DimensionalityTraits<TopoDS_Shape>
{
    static constexpr int geomDim = 2;
};

/*!
 * \brief Specialization for Brep solids.
 */
template<>
struct DimensionalityTraits<TopoDS_Solid>
: public DimensionalityTraits<TopoDS_Shape>
{
    static constexpr int geomDim = 3;
};

/*!
 * \brief Free function to return the dimension of a geometry.
 */
template<class Geom>
constexpr int getDimension(const Geom& g)
{ return DimensionalityTraits<Geom>::geomDim; }

/*!
 * \brief Free function to return the dimension of the space
 *        in which a geometry is described in.
 */
template<class Geom>
constexpr int getWorldDimension(const Geom& g)
{ return DimensionalityTraits<Geom>::worldDim; }

} // end namespace Frackit

#endif // FRACKIT_COMMON_EXTRACT_DIMENSION_HH
