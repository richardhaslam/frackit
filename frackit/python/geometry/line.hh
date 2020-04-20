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
#ifndef FRACKIT_PYTHON_GEOMETRY_LINE_HH
#define FRACKIT_PYTHON_GEOMETRY_LINE_HH

#include <string>

#include <pybind11/pybind11.h>

#include <frackit/geometry/geometry.hh>
#include <frackit/geometry/line.hh>
#include "registerdimensionproperties.hh"

namespace Frackit::Python {

namespace py = pybind11;

namespace Detail {

    template<class ctype, int worldDim>
    void registerLine(py::module& module)
    {
        using namespace py::literals;
        using Line = Line<ctype, worldDim>;
        using Point = typename Line::Point;
        using Direction = typename Line::Direction;

        // register class and define constructor
        std::string className("Line_" + std::to_string(worldDim));
        py::class_<Line, Geometry> cls(module, className.c_str());
        cls.def(py::init<const Point&, const Direction&>(), "supportPoint"_a, "direction"_a);

        // dimensionality properties
        registerDimensionProperties(cls);

        // member functions
        cls.def("name", &Line::name, "name of the geometry");
        cls.def("direction", &Line::direction, "direction of the line");
        cls.def("supportingPoint", &Line::supportingPoint, "supporting point on the line");
        cls.def("projection", &Line::projection, "p"_a, "project a point on the line");
        cls.def("contains", py::overload_cast<const Point&, ctype>(&Line::contains, py::const_),
                "point"_a, "eps"_a, "returns true if the given point is on the line (given tolerance)");
        cls.def("contains", py::overload_cast<const Point&>(&Line::contains, py::const_),
                "point"_a, "returns true if the given point is on the line (default tolerance)");
    }

} // end namespace detail

template<class ctype>
void registerLine(py::module& module)
{
    Detail::registerLine<ctype, 1>(module);
    Detail::registerLine<ctype, 2>(module);
    Detail::registerLine<ctype, 3>(module);
}

} // end namespace Frackit::Python

#endif
