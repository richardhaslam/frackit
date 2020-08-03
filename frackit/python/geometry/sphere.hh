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
#ifndef FRACKIT_PYTHON_GEOMETRY_SPHERE_HH
#define FRACKIT_PYTHON_GEOMETRY_SPHERE_HH

#include <pybind11/pybind11.h>

#include <frackit/geometry/geometry.hh>
#include <frackit/geometry/sphere.hh>
#include "registerdimensionproperties.hh"

namespace Frackit::Python {

namespace py = pybind11;

template<class ctype>
void registerSphere(py::module& module)
{
    using namespace py::literals;
    using Sphere = Sphere<ctype>;
    using Point = typename Sphere::Point;

    // define constructors
    py::class_<Sphere, Geometry> cls(module, "Sphere");
    cls.def(py::init<const Point&, ctype>(), "center"_a, "radius"_a);

    // dimensionality properties
    registerDimensionProperties(cls);

    // getter functions
    cls.def("name", &Sphere::name, "name of the geometry class");
    cls.def("center", &Sphere::center, "returns the center point of the sphere");
    cls.def("radius", &Sphere::radius, "returns the radius of the sphere");
    cls.def("volume", &Sphere::volume, "returns the volume of the sphere");

    // contains queries
    cls.def("contains", py::overload_cast<const Point&>(&Sphere::contains, py::const_),
            "point"_a, "returns true if the given point is within the sphere (default tolerance)");
    cls.def("contains", py::overload_cast<const Point&, ctype>(&Sphere::contains, py::const_),
            "point"_a, "eps"_a, "returns true if the given point is within the sphere (given tolerance)");

    cls.def("__repr__", [&] (const Sphere& s) { return "Frackit::Sphere"; });
}

} // end namespace Frackit::Python

#endif
