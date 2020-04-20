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
#ifndef FRACKIT_PYTHON_GEOMETRY_PLANE_HH
#define FRACKIT_PYTHON_GEOMETRY_PLANE_HH

#include <string>

#include <pybind11/pybind11.h>

#include <frackit/geometry/geometry.hh>
#include <frackit/geometry/plane.hh>
#include "registerdimensionproperties.hh"

namespace Frackit::Python {

namespace py = pybind11;

namespace Detail {

    template<class ctype, int worldDim>
    void registerPlane(py::module& module)
    {
        using namespace py::literals;
        using Plane = Plane<ctype, worldDim>;
        using Point = typename Plane::Point;
        using Line = typename Plane::Line;
        using Segment = typename Plane::Segment;
        using Direction = typename Plane::Direction;

        // register class and define constructors
        std::string className("Plane_" + std::to_string(worldDim));
        py::class_<Plane, Geometry> cls(module, className.c_str());
        cls.def(py::init<>());
        cls.def(py::init<const Point&, const Direction&>(), "supportPoint"_a, "normal"_a);
        cls.def(py::init<const Point&, const Point&,  const Point&>(), "p1"_a, "p2"_a, "p3"_a);
        cls.def(py::init<const Point&, const Direction&,  const Direction&, const Direction&>(),
                "supportPoint"_a, "base1"_a, "base2"_a, "normal"_a);

        // dimensionality properties
        registerDimensionProperties(cls);

        // getter functions
        cls.def("name", &Plane::name, "name of the geometry");
        cls.def("normal", &Plane::normal, "normal direction");
        cls.def("base1", &Plane::base1, "first in-plane basis vector");
        cls.def("base2", &Plane::base2, "second in-plane basis vector");
        cls.def("supportingPoint", &Plane::supportingPoint, "supporting point on the plane");

        // projections
        cls.def("projection", py::overload_cast<const Point&>(&Plane::projection, py::const_), "point"_a, "project a point on the plane");
        cls.def("projection", py::overload_cast<const Segment&>(&Plane::projection, py::const_), "segment"_a, "project a segment on the plane");
        cls.def("projection", py::overload_cast<const Line&>(&Plane::projection, py::const_), "line"_a, "project a line on the plane");

        // contains query
        cls.def("contains", &Plane::contains, "point"_a, "eps"_a = Precision<ctype>::confusion(),
                "returns true if the given point is on the plane");

        cls.def("__repr__", [&] (const Plane& p) { return "Frackit::Plane<" + std::to_string(worldDim) + ">"; });
    }

} // end namespace detail

template<class ctype>
void registerPlane(py::module& module)
{
    Detail::registerPlane<ctype, 3>(module);
}

} // end namespace Frackit::Python

#endif
