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
#ifndef FRACKIT_PYTHON_GEOMETRY_BOX_HH
#define FRACKIT_PYTHON_GEOMETRY_BOX_HH

#include <pybind11/pybind11.h>

#include <frackit/geometry/geometry.hh>
#include <frackit/geometry/box.hh>
#include "registerdimensionproperties.hh"

namespace Frackit::Python {

namespace py = pybind11;

template<class ctype>
void registerBox(py::module& module)
{
    using namespace py::literals;
    using Box = Box<ctype>;
    using Point = typename Box::Point;

    // define constructors
    py::class_<Box, Geometry> cls(module, "Box");
    cls.def(py::init<ctype, ctype, ctype, ctype, ctype, ctype>(),
            "xmin"_a, "ymin"_a, "zmin"_a, "xmax"_a, "ymax"_a, "zmax"_a);
    cls.def(py::init<const Point&, Point>(), "firstCorner"_a, "secondCorner"_a);

    // dimensionality properties
    registerDimensionProperties(cls);

    // getter functions
    cls.def("name", &Box::name, "name of the geometry class");
    cls.def("volume", &Box::volume, "returns the volume of the box");
    cls.def("xMin", &Box::xMin, "returns the minimum x-coordinate of the box");
    cls.def("yMin", &Box::yMin, "returns the minimum y-coordinate of the box");
    cls.def("zMin", &Box::zMin, "returns the minimum z-coordinate of the box");
    cls.def("xMax", &Box::xMax, "returns the maximum x-coordinate of the box");
    cls.def("yMax", &Box::yMax, "returns the maximum y-coordinate of the box");
    cls.def("zMax", &Box::zMax, "returns the maximum z-coordinate of the box");

    cls.def_static("numCorners", &Box::numCorners, "returns the number of corners of the box");
    cls.def("corner", &Box::corner, "cornerIdx"_a, "returns the corner with the provided index");
    cls.def("center", &Box::center, "returns the center of the box");

    cls.def_static("numEdges", &Box::numEdges, "returns the number of edges of the box");
    cls.def("edge", &Box::edge, "edgeIdx"_a, "returns the edge with the provided index");

    cls.def_static("numFaces", &Box::numFaces, "returns the number of faces of the box");
    cls.def("face", &Box::face, "faceIdx"_a, "returns the face with the provided index");

    // contains queries
    cls.def("contains", py::overload_cast<const Point&>(&Box::contains, py::const_),
            "point"_a, "returns true if the given point is within the box (default tolerance)");
    cls.def("contains", py::overload_cast<const Point&, ctype>(&Box::contains, py::const_),
            "point"_a, "eps"_a, "returns true if the given point is within the box (given tolerance)");

    cls.def("__repr__", [&] (const Box& b) { return "Frackit::Box"; });
}

} // end namespace Frackit::Python

#endif
