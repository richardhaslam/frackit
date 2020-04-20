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
#ifndef FRACKIT_PYTHON_GEOMETRY_CIRCLE_HH
#define FRACKIT_PYTHON_GEOMETRY_CIRCLE_HH

#include <string>

#include <pybind11/pybind11.h>

#include <frackit/geometry/geometry.hh>
#include <frackit/geometry/ellipticalgeometry.hh>
#include <frackit/geometry/circle.hh>
#include "registerdimensionproperties.hh"

namespace Frackit::Python {

namespace py = pybind11;

namespace Detail {

    template<class ctype, int worldDim>
    void registerCircle(py::module& module)
    {
        using Circle = Circle<ctype, worldDim>;
        using EllipticalGeometry = EllipticalGeometry<ctype, worldDim>;
        std::string className("Circle_" + std::to_string(worldDim));
        py::class_<Circle, Geometry, EllipticalGeometry> cls(module, className.c_str());

        // dimensionality properties
        registerDimensionProperties(cls);

        // define constructors
        using namespace py::literals;
        using Point = typename Circle::Point;
        using Direction = typename Circle::Direction;

        cls.def(py::init<const Point&, const Direction&, ctype>(),
                "center"_a, "normal"_a, "radius"_a);
        cls.def("name", &Circle::name, "name of the geometry class");
        cls.def("radius", &Circle::radius, "returns the radius of the circle");
        cls.def("base1", &Circle::base1, "returns the first unit basis vector");
        cls.def("base2", &Circle::base2, "returns the second unit basis vector");

        // contains queries
        cls.def("contains", py::overload_cast<const Point&>(&Circle::contains, py::const_),
                "point"_a, "returns true if the given point is on the circle (default tolerance)");
        cls.def("contains", py::overload_cast<const Point&, ctype>(&Circle::contains, py::const_),
                "point"_a, "eps"_a, "returns true if the given point is on the circle (given tolerance)");

        // get a point from local coordinate or angle
        cls.def("getPoint", &Circle::getPoint, "localCoordinate"_a,
                "return the point on the circle for the given local coordinate (0.0 <= localCoordinate <= 1)");
        cls.def("getPointFromAngle", &Circle::getPointFromAngle, "angle"_a,
                "return the point on the circle for the given angle w.r.t. the circle-local coordinate system (angle >= 0.0)");

        cls.def("__repr__", [&] (const Circle& c) { return "Frackit::Circle<" + std::to_string(worldDim) + ">"; });
    }

} // end namespace detail

template<class ctype>
void registerCircle(py::module& module)
{
    Detail::registerCircle<ctype, 3>(module);
}

} // end namespace Frackit::Python

#endif
