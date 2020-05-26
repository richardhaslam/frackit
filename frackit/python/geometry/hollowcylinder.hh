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
#ifndef FRACKIT_PYTHON_GEOMETRY_HOLLOW_CYLINDER_HH
#define FRACKIT_PYTHON_GEOMETRY_HOLLOW_CYLINDER_HH

#include <string>

#include <pybind11/pybind11.h>
#include <frackit/geometry/geometry.hh>
#include <frackit/geometry/hollowcylinder.hh>
#include "registerdimensionproperties.hh"

namespace Frackit::Python {

namespace py = pybind11;

template<class ctype>
void registerHollowCylinder(py::module& module)
{
    using namespace py::literals;
    using HollowCylinder = HollowCylinder<ctype>;
    using Point = typename HollowCylinder::Point;
    using Direction = typename HollowCylinder::Direction;

    // register class and define constructorss
    py::class_<HollowCylinder, Geometry> cls(module, "HollowCylinder");
    cls.def(py::init<ctype, ctype, ctype>(), "innerRadius"_a, "outerRadius"_a, "height"_a);
    cls.def(py::init<const Point&, const Direction&, ctype, ctype, ctype>(),
            "bottomCenter"_a, "axis"_a, "innerRadius"_a, "outerRadius"_a, "height"_a);

    // dimensionality properties
    registerDimensionProperties(cls);

    // getter functions
    cls.def("name", &HollowCylinder::name, "name of the geometry class");

    cls.def("base1", &HollowCylinder::base1, "first basis vector orthogonal to center line");
    cls.def("base2", &HollowCylinder::base2, "second basis vector orthogonal to center line");
    cls.def("base3", &HollowCylinder::base3, "basis vector parallel to center line");
    cls.def("axis", &HollowCylinder::axis, "direction of the cylinder axis");

    cls.def("bottomInnerCircle", &HollowCylinder::bottomInnerCircle, "returns the inner circle of the bottom face");
    cls.def("topInnerCircle", &HollowCylinder::topInnerCircle, "returns the inner circle of the top face");
    cls.def("bottomOuterCircle", &HollowCylinder::bottomOuterCircle, "returns the outer circle of the bottom face");
    cls.def("topOuterCircle", &HollowCylinder::topOuterCircle, "returns the outer circle of the top face");

    cls.def("innerLateralFace", &HollowCylinder::innerLateralFace, "returns the inner lateral surface (mantle) of the hollow cylinder");
    cls.def("outerLateralFace", &HollowCylinder::outerLateralFace, "returns the outer lateral surface (mantle) of the hollow cylinder");

    cls.def("centerSegment", &HollowCylinder::centerSegment, "returns the segment describing the cylinder axis");
    cls.def("fullCylinder", &HollowCylinder::fullCylinder, "returns the full cylinder this hollow cylinder is a part of");

    cls.def("height", &HollowCylinder::height, "returns the height of the hollow cylinder");
    cls.def("innerRadius", &HollowCylinder::innerRadius, "returns the inner radius of the hollow cylinder");
    cls.def("outerRadius", &HollowCylinder::outerRadius, "returns the outer radius of the hollow cylinder");
    cls.def("volume", &HollowCylinder::volume, "returns the volume of the hollow cylinder");

    // contains queries
    cls.def("contains", py::overload_cast<const Point&>(&HollowCylinder::contains, py::const_),
            "point"_a, "returns true if the given point is within the hollow cylinder (default tolerance)");
    cls.def("contains", py::overload_cast<const Point&, ctype>(&HollowCylinder::contains, py::const_),
            "point"_a, "eps"_a, "returns true if the given point is within the hollow cylinder (given tolerance)");

    cls.def("__repr__", [&] (const HollowCylinder& c) { return "Frackit::HollowCylinder"; });
}

} // end namespace Frackit::Python

#endif
