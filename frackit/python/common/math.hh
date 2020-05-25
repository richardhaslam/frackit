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
#ifndef FRACKIT_PYTHON_COMMON_MATH_HH
#define FRACKIT_PYTHON_COMMON_MATH_HH

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <frackit/geometry/vector.hh>
#include <frackit/geometry/direction.hh>
#include <frackit/common/math.hh>

namespace Frackit::Python {

namespace py = pybind11;

template<class ctype>
void registerMath(py::module& module)
{
    module.def("toDegrees", &Frackit::toDegrees<ctype>, "Converts radians into degrees");
    module.def("toRadians", &Frackit::toRadians<ctype>, "Converts degrees into radians");

    // Register rotation overload for single vector
    using namespace py::literals;
    using Vector_3 = Vector<ctype, 3>;
    using Direction_3 = Direction<ctype, 3>;
    module.def("rotate",
               py::overload_cast<Vector_3&, const Direction_3&, ctype>(&rotate<ctype>),
               "Rotates a vector around the given axis & angle",
               "vector"_a, "axis"_a, "angle"_a);
}

} // end namespace Frackit::Python

#endif
