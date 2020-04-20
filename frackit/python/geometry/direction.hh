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
#ifndef FRACKIT_PYTHON_GEOMETRY_DIRECTION_HH
#define FRACKIT_PYTHON_GEOMETRY_DIRECTION_HH

#include <string>

#include <pybind11/pybind11.h>
#include <frackit/geometry/geometry.hh>
#include <frackit/geometry/direction.hh>
#include "registerdimensionproperties.hh"

namespace Frackit::Python {

namespace py = pybind11;

namespace Detail {

    template<class ctype, int worldDim>
    void registerDirection(py::module& module)
    {
        using namespace py::literals;
        using Direction = Direction<ctype, worldDim>;
        using Vector = typename Direction::Vector;

        // register class and define constructors
        std::string className("Direction_" + std::to_string(worldDim));
        py::class_<Direction, Geometry> cls(module, className.c_str());
        cls.def(py::init<>());
        cls.def(py::init<const Vector&>(), "vector"_a);

        // dimensionality properties
        registerDimensionProperties(cls);

        // member functions
        cls.def("name", &Direction::name, "name of the geometry class");
        cls.def("invert", &Direction::invert, "inverts the direction");

        // returns a vector in this direction with provided length
        cls.def(py::self * ctype());

        // define retrieval functions for the coordinates
        cls.def("x", &Direction::x, "x-coordinate of the direction");
        if constexpr (worldDim > 1)
            cls.def("y", &Direction::y, "y-coordinate of the direction");
        if constexpr (worldDim > 2)
            cls.def("z", &Direction::z, "z-coordinate of the direction");
    }

} // end namespace detail

template<class ctype>
void registerDirection(py::module& module)
{
    Detail::registerDirection<ctype, 1>(module);
    Detail::registerDirection<ctype, 2>(module);
    Detail::registerDirection<ctype, 3>(module);
}

} // end namespace Frackit::Python

#endif
