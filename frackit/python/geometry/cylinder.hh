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
#ifndef FRACKIT_PYTHON_GEOMETRY_CYLINDER_HH
#define FRACKIT_PYTHON_GEOMETRY_CYLINDER_HH

#include <string>

#include <pybind11/pybind11.h>
#include <frackit/geometry/geometry.hh>
#include <frackit/geometry/cylinder.hh>

namespace Frackit::Python {

namespace py = pybind11;

template<class ctype>
void registerCylinder(py::module& module)
{
    using namespace py::literals;
    using Cylinder = Cylinder<ctype>;
    using Circle = typename Cylinder::Circle;
    using Disk = typename Cylinder::Disk;
    using Point = typename Cylinder::Point;

    // register class and define constructorss
    py::class_<Cylinder, Geometry> cls(module, "Cylinder");
    cls.def(py::init<ctype, ctype>(), "radius"_a, "height"_a);
    cls.def(py::init<const Circle&, ctype>(), "bottomCircle"_a, "height"_a);
    cls.def(py::init<const Disk&, ctype>(), "bottom"_a, "height"_a);

    // getter functions
    cls.def("name", &Cylinder::name, "name of the geometry class");
    cls.def("base1", &Cylinder::base1, "first basis vector orthogonal to center line");
    cls.def("base2", &Cylinder::base2, "second basis vector orthogonal to center line");
    cls.def("base3", &Cylinder::base3, "basis vector parallel to center line");
    cls.def("direction", &Cylinder::direction, "direction of the cylinder axis");

    cls.def("topFace", &Cylinder::topFace, "returns the upper face");
    cls.def("bottomFace", &Cylinder::bottomFace, "returns the lower face");
    cls.def("lateralFace", &Cylinder::lateralFace, "returns the lateral cylinder surface");
    cls.def("centerSegment", &Cylinder::centerSegment, "returns the segment describing the cylinder axis");

    cls.def("height", &Cylinder::height, "returns the height of the cylinder");
    cls.def("radius", &Cylinder::radius, "returns the radius of the cylinder");
    cls.def("volume", &Cylinder::volume, "returns the volume of the cylinder");

    // contains queries
    cls.def("contains", py::overload_cast<const Point&>(&Cylinder::contains, py::const_),
            "point"_a, "returns true if the given point is within the cylinder (default tolerance)");
    cls.def("contains", py::overload_cast<const Point&, ctype>(&Cylinder::contains, py::const_),
            "point"_a, "eps"_a, "returns true if the given point is within the cylinder (given tolerance)");
}

} // end namespace Frackit::Python

#endif
