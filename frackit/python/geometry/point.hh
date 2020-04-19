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
#ifndef FRACKIT_PYTHON_GEOMETRY_POINT_HH
#define FRACKIT_PYTHON_GEOMETRY_POINT_HH

#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <frackit/geometry/geometry.hh>
#include <frackit/geometry/point.hh>

namespace Frackit::Python {

namespace py = pybind11;

namespace Detail {

    template<class Base, class Impl>
    void registerPointBase(py::module& module)
    {
        using ctype = typename Base::ctype;
        using Vector = Frackit::Vector<ctype, Impl::worldDimension()>;

        // register class and define constructor
        std::string className("PointBase_" + std::to_string(Base::worldDimension()));
        py::class_<Base, Geometry> cls(module, className.c_str());
        cls.def(py::init<>());

        // member function
        cls.def("name", &Base::name, "name of the geometry class");

        // addition/subtraction overloads
        cls.def(py::self += Vector());
        cls.def(py::self + Vector());
        cls.def(py::self -= Vector());
        cls.def(py::self - Vector());

        // equality checks
        using namespace py::literals;
        cls.def("isEqual",
                py::overload_cast<const Impl&>(&Base::isEqual, py::const_),
                "p"_a, "equality check with default tolerance");
        cls.def("isEqual",
                py::overload_cast<const Impl&, ctype>(&Base::isEqual, py::const_),
                "p"_a, "b"_a, "equality check with default tolerance");
    }

    template<class ctype, int worldDim>
    void registerPoint(py::module& module)
    {
        using namespace py::literals;
        using Point = Point<ctype, worldDim>;
        using Base = PointBase<Point, ctype, worldDim>;
        registerPointBase<Base, Point>(module);

        // register class and define constructors
        std::string className("Point_" + std::to_string(worldDim));
        py::class_<Point, Base> cls(module, className.c_str());
        cls.def(py::init<>());
        if constexpr (worldDim == 1)
            cls.def(py::init<ctype>(), "x"_a);
        else if constexpr (worldDim == 2)
            cls.def(py::init<ctype, ctype>(), "x"_a, "y"_a);
        else if constexpr (worldDim == 3)
            cls.def(py::init<ctype, ctype, ctype>(), "x"_a, "y"_a, "z"_a);

        // define retrieval functions for the coordinates
        cls.def("x", &Point::x, "x-coordinate of the point");
        if constexpr (worldDim > 1)
            cls.def("y", &Point::y, "y-coordinate of the point");
        if constexpr (worldDim > 2)
            cls.def("z", &Point::z, "z-coordinate of the point");
    }

} // end namespace detail

template<class ctype>
void registerPoint(py::module& module)
{
    Detail::registerPoint<ctype, 1>(module);
    Detail::registerPoint<ctype, 2>(module);
    Detail::registerPoint<ctype, 3>(module);
}

} // end namespace Frackit::Python

#endif
