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
#ifndef FRACKIT_PYTHON_GEOMETRY_QUADRILATERAL_HH
#define FRACKIT_PYTHON_GEOMETRY_QUADRILATERAL_HH

#include <string>

#include <pybind11/pybind11.h>

#include <frackit/geometry/geometry.hh>
#include <frackit/geometry/quadrilateral.hh>
#include "registerdimensionproperties.hh"

namespace Frackit::Python {

namespace py = pybind11;

namespace Detail {

    template<class ctype, int worldDim>
    void registerQuadrilateral(py::module& module)
    {
        using namespace py::literals;
        using Quad = Quadrilateral<ctype, worldDim>;
        using Point = typename Quad::Point;

        // register class and define constructors
        std::string className("Quadrilateral_" + std::to_string(worldDim));
        py::class_<Quad, Geometry> cls(module, className.c_str());
        cls.def(py::init<const Point&, const Point&, const Point&, const Point&>(),
                "p1"_a, "p2"_a, "p3"_a, "p4"_a);

        // dimensionality properties
        registerDimensionProperties(cls);

        // member functions
        cls.def("name", &Quad::name, "name of the geometry");
        cls.def("area", &Quad::area, "area of the quadrilateral");
        cls.def("center", &Quad::center, "center point of the quadrilateral");
        cls.def("supportingPlane", &Quad::supportingPlane, "returns the supporting plane");
        cls.def("isConvex", &Quad::isConvex, "returns true if the quadrilateral is convex");

        cls.def_static("numCorners", &Quad::numCorners, "returns the number of corners of the quadrilateral");
        cls.def("corner", &Quad::corner, "cornerIdx"_a, "returns the corner with the provided index");

        cls.def_static("numEdges", &Quad::numEdges, "returns the number of edges of the quadrilateral");
        cls.def("edge", &Quad::edge, "edgeIdx"_a, "returns the edge with the provided index");

        // contains queries
        cls.def("contains", py::overload_cast<const Point&, bool>(&Quad::contains, py::const_),
                "point"_a, "checkIfOnPlane"_a = true,
                "returns true if the given point is on the quadrilateral (default tolerance)");
        cls.def("contains", py::overload_cast<const Point&, ctype, bool>(&Quad::contains, py::const_),
                "point"_a, "eps"_a, "checkIfOnPlane"_a = true,
                "returns true if the given point is on the quadrilateral (given tolerance)");

        cls.def("__repr__", [&] (const Quad& q) { return "Frackit::Quadrilateral<" + std::to_string(worldDim) + ">"; });
    }

} // end namespace detail

template<class ctype>
void registerQuadrilateral(py::module& module)
{
    Detail::registerQuadrilateral<ctype, 3>(module);
}

} // end namespace Frackit::Python

#endif
