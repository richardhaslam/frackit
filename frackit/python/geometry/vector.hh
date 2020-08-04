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
#ifndef FRACKIT_PYTHON_GEOMETRY_VECTOR_HH
#define FRACKIT_PYTHON_GEOMETRY_VECTOR_HH

#include <array>
#include <string>
#include <algorithm>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <frackit/precision/precision.hh>
#include <frackit/geometry/geometry.hh>
#include <frackit/geometry/vector.hh>
#include "registerdimensionproperties.hh"

namespace Frackit::Python {

namespace py = pybind11;

namespace Detail {

    template<class Base, class Impl>
    void registerVectorBase(py::module& module)
    {
        using ctype = typename Base::ctype;

        std::string className("VectorBase_" + std::to_string(Base::worldDimension()));
        py::class_<Base, Geometry> cls(module, className.c_str());

        // dimensionality properties
        registerDimensionProperties(cls);

        cls.def("name", &Base::name, "name of the geometry class");
        cls.def("squaredLength", &Base::squaredLength, "squared length of the vector");
        cls.def("length", &Base::length, "length of the vector");

        using namespace py::literals;
        cls.def("isParallel", &Base::isParallel, "v"_a, "eps"_a = Precision<ctype>::angular(),
                "returns true if the given vector is parallel to this one");

        cls.def(py::self * Impl());   // scalar product
        cls.def(py::self + Impl());   // vector addition
        cls.def(py::self - Impl());   // vector subtraction
        cls.def(py::self *= ctype()); // multiplication by scalar
        cls.def(py::self /= ctype()); // division by scalar
    }

    template<class ctype, int worldDim>
    void registerVector(py::module& module)
    {
        using Vector = Vector<ctype, worldDim>;
        using Base = VectorBase<Vector, ctype, worldDim>;
        registerVectorBase<Base, Vector>(module);

        std::string className("Vector_" + std::to_string(worldDim));
        py::class_<Vector, Base> cls(module, className.c_str());

        // define constructors
        using namespace py::literals;
        using Point = typename Vector::Point;
        using Direction = typename Vector::Direction;

        cls.def(py::init<>());
        cls.def(py::init<const Point&, const Point&>(), "source"_a, "target"_a);
        cls.def(py::init<const Direction&>(), "direction"_a);
        cls.def(py::init( [] (py::list x) -> Vector
                          {
                            // create array of matching size and fill with values
                            // from the list. Use zeros at the end if list is shorter.
                            std::array<ctype, worldDim> vals;
                            std::fill(vals.begin(), vals.end(), 0.0);
                            unsigned int minDim = std::min(worldDim, static_cast<int>(x.size()));
                            for (unsigned int i = 0; i < minDim; ++i)
                                vals[i] = x[i].template cast<ctype>();

                            if constexpr (worldDim == 1)
                                return {vals[0]};
                            else if constexpr (worldDim == 2)
                                return {vals[0], vals[1]};
                            else
                                return {vals[0], vals[1], vals[2]};
                          } ), "coordinates"_a);

        if constexpr (worldDim == 1)
            cls.def(py::init<ctype>(), "x"_a);
        else if constexpr (worldDim == 2)
            cls.def(py::init<ctype, ctype>(), "x"_a, "y"_a);
        else if constexpr (worldDim == 3)
            cls.def(py::init<ctype, ctype, ctype>(), "x"_a, "y"_a, "z"_a);

        // define retrieval functions for the coordinates
        cls.def("x", &Vector::x, "x-coordinate of the vector");
        if constexpr (worldDim > 1) cls.def("y", &Vector::y, "y-coordinate of the vector");
        if constexpr (worldDim > 2) cls.def("z", &Vector::z, "z-coordinate of the vector");

        using std::to_string;
        if constexpr (worldDim == 1)
            cls.def("__repr__", [&] (const Vector& v)
                                    { return "Frackit::Vector<1> (" + to_string(v.x()) + ")"; });
        else if constexpr (worldDim == 2)
            cls.def("__repr__", [&] (const Vector& v)
                                    { return "Frackit::Vector<2> (" + to_string(v.x()) + ", "
                                                                    + to_string(v.y()) + ")"; });
        else
            cls.def("__repr__", [&] (const Vector& v)
                                    { return "Frackit::Vector<3> (" + to_string(v.x()) + ", "
                                                                    + to_string(v.y()) + ", "
                                                                    + to_string(v.z()) + ")"; });
    }

} // end namespace detail

template<class ctype>
void registerVector(py::module& module)
{
    Detail::registerVector<ctype, 1>(module);
    Detail::registerVector<ctype, 2>(module);
    Detail::registerVector<ctype, 3>(module);
}

} // end namespace Frackit::Python

#endif
