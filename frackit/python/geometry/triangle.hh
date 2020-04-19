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
#ifndef FRACKIT_PYTHON_GEOMETRY_TRIANGLE_HH
#define FRACKIT_PYTHON_GEOMETRY_TRIANGLE_HH

#include <string>

#include <pybind11/pybind11.h>

#include <frackit/geometry/geometry.hh>
#include <frackit/geometry/triangle.hh>

namespace Frackit::Python {

namespace py = pybind11;

namespace Detail {

    template<class ctype, int worldDim>
    void registerTriangle(py::module& module)
    {
        using namespace py::literals;
        using Triangle = Triangle<ctype, worldDim>;
        using Point = typename Triangle::Point;

        // register class and define constructor
        std::string className("Triangle_" + std::to_string(worldDim));
        py::class_<Triangle, Geometry> cls(module, className.c_str());
        cls.def(py::init<const Point&, const Point&, const Point&>(), "p1"_a, "p2"_a, "p3"_a);

        // member functions
        cls.def("name", &Triangle::name, "name of the geometry");
        cls.def("area", &Triangle::area, "area of the triangle");
        cls.def("center", &Triangle::center, "center point of the triangle");
        cls.def("normal", &Triangle::normal, "direction normal to the triangle plane");
        cls.def("supportingPlane", &Triangle::supportingPlane, "returns the supporting plane");

        cls.def_static("numCorners", &Triangle::numCorners, "returns the number of corners of the triangle");
        cls.def("corner", &Triangle::corner, "cornerIdx"_a, "returns the corner with the provided index");

        cls.def_static("numEdges", &Triangle::numEdges, "returns the number of edges of the triangle");
        cls.def("edge", &Triangle::edge, "edgeIdx"_a, "returns the edge with the provided index");

        // contains queries
        cls.def("contains", py::overload_cast<const Point&, bool>(&Triangle::contains, py::const_),
                "point"_a, "checkIfOnPlane"_a = true,
                "returns true if the given point is on the triangle (given tolerance)");
        cls.def("contains", py::overload_cast<const Point&, ctype, bool>(&Triangle::contains, py::const_),
                "point"_a, "eps"_a, "checkIfOnPlane"_a = true,
                "returns true if the given point is on the triangle (given tolerance)");
    }

} // end namespace detail

template<class ctype>
void registerTriangle(py::module& module)
{
    Detail::registerTriangle<ctype, 3>(module);
}

} // end namespace Frackit::Python

#endif
