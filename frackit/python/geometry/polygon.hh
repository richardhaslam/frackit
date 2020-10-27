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
#ifndef FRACKIT_PYTHON_GEOMETRY_POLYGON_HH
#define FRACKIT_PYTHON_GEOMETRY_POLYGON_HH

#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <frackit/geometry/geometry.hh>
#include <frackit/geometry/polygon.hh>
#include "registerdimensionproperties.hh"

namespace Frackit::Python {

namespace py = pybind11;

namespace Detail {

    template<class ctype, int worldDim>
    void registerPolygon(py::module& module)
    {
        using namespace py::literals;
        using Polygon = Frackit::Polygon<ctype, worldDim>;
        using Point = typename Polygon::Point;

        // register class and define constructors
        std::string className("Polygon_" + std::to_string(worldDim));
        py::class_<Polygon, Geometry> cls(module, className.c_str());
        cls.def(py::init<const std::vector<Point>&>(), "corners"_a);

        // dimensionality properties
        registerDimensionProperties(cls);

        // member functions
        cls.def("name", &Polygon::name, "name of the geometry");
        cls.def("area", &Polygon::area, "area of the polygon");
        cls.def("center", &Polygon::center, "centroid of the polygon");
        cls.def("supportingPlane", &Polygon::supportingPlane, "returns the supporting plane");
        cls.def("isConvex", &Polygon::isConvex, "returns true if the polygon is convex");

        cls.def("numCorners", &Polygon::numCorners, "returns the number of corners of the polygon");
        cls.def("corner", &Polygon::corner, "cornerIdx"_a, "returns the corner with the provided index");

        cls.def("numEdges", &Polygon::numEdges, "returns the number of edges of the polygon");
        cls.def("edge", &Polygon::edge, "edgeIdx"_a, "returns the edge with the provided index");

        // contains queries
        cls.def("contains", py::overload_cast<const Point&, bool>(&Polygon::contains, py::const_),
                "point"_a, "checkIfOnPlane"_a = true,
                "returns true if the given point is on the polygon (default tolerance)");
        cls.def("contains", py::overload_cast<const Point&, ctype, bool>(&Polygon::contains, py::const_),
                "point"_a, "eps"_a, "checkIfOnPlane"_a = true,
                "returns true if the given point is on the polygon (given tolerance)");

        cls.def("__repr__", [&] (const Polygon& p) { return "Frackit::Polygon<" + std::to_string(worldDim) + ">"; });
    }

} // end namespace detail

template<class ctype>
void registerPolygon(py::module& module)
{
    Detail::registerPolygon<ctype, 3>(module);
}

} // end namespace Frackit::Python

#endif
