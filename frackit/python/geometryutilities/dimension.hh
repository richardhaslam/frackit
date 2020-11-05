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
#ifndef FRACKIT_PYTHON_GEOMETRY_DIMENSION_HH
#define FRACKIT_PYTHON_GEOMETRY_DIMENSION_HH

#include <type_traits>

#include <TopoDS_Shape.hxx>
#include <TopoDS_Compound.hxx>

#include <frackit/common/extractdimension.hh>
#include <frackit/python/geometry/brepwrapper.hh>

namespace Frackit::Python {

//! import traits from Frackit namespace
template<class Geometry>
struct DimensionalityTraits
: public Frackit::DimensionalityTraits<Geometry>
{};

//! dimensionality traits for brep shape wrappers
template<class Shape>
struct DimensionalityTraits<BRepWrapper<Shape>>
: public Frackit::DimensionalityTraits<Shape>
{};

//! Helper struct to check if a geometry has dimensionality known at compile-time
template<class Geometry>
struct HasFixedDimensionality : public std::true_type {};
template<>
struct HasFixedDimensionality<TopoDS_Shape> : public std::false_type {};
template<>
struct HasFixedDimensionality<TopoDS_Compound> : public std::false_type {};

template<class Shape>
struct HasFixedDimensionality<BRepWrapper<Shape>>
: public HasFixedDimensionality<Shape> {};

} // end namespace Frackit::Python

#endif
