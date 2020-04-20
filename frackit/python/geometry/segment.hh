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
#ifndef FRACKIT_PYTHON_GEOMETRY_SEGMENT_HH
#define FRACKIT_PYTHON_GEOMETRY_SEGMENT_HH

#include <string>

#include <pybind11/pybind11.h>

#include <frackit/geometry/geometry.hh>
#include <frackit/geometry/segment.hh>
#include "registerdimensionproperties.hh"

namespace Frackit::Python {

namespace py = pybind11;

namespace Detail {

    template<class ctype, int worldDim>
    void registerSegment(py::module& module)
    {
        using namespace py::literals;
        using Segment = Segment<ctype, worldDim>;
        using Point = typename Segment::Point;

        // register class and define constructor
        std::string className("Segment_" + std::to_string(worldDim));
        py::class_<Segment, Geometry> cls(module, className.c_str());
        cls.def(py::init<const Point&, const Point&>(), "source"_a, "target"_a);

        // dimensionality properties
        registerDimensionProperties(cls);

        // member functions
        cls.def("name", &Segment::name, "name of the geometry");
        cls.def("source", &Segment::source, "source point of the segment");
        cls.def("target", &Segment::target, "target point of the segment");
        cls.def("direction", &Segment::direction, "direction of the segment");
        cls.def("length", &Segment::length, "length of the segment");
        cls.def("supportingLine", &Segment::supportingLine, "returns the supporting line");

        // obtain a point on the segment between 0 <= localCoordinate <= 1
        cls.def("getPoint", &Segment::getPoint, "localCoordinate"_a,
                "returns the point on the line for the given local coordinate (0 <= localCoordinate <= 1)");

        // contains queries
        cls.def("contains", py::overload_cast<const Point&, bool>(&Segment::contains, py::const_),
                "point"_a, "checkIfOnLine"_a = true,
                "returns true if the given point is on the segment (default tolerance)");
        cls.def("contains", py::overload_cast<const Point&, ctype, bool>(&Segment::contains, py::const_),
                "point"_a, "eps"_a, "checkIfOnLine"_a = true,
                "returns true if the given point is on the segment (given tolerance)");
    }

} // end namespace detail

template<class ctype>
void registerSegment(py::module& module)
{
    Detail::registerSegment<ctype, 1>(module);
    Detail::registerSegment<ctype, 2>(module);
    Detail::registerSegment<ctype, 3>(module);
}

} // end namespace Frackit::Python

#endif
