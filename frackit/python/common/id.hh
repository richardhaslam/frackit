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
#ifndef FRACKIT_PYTHON_COMMON_ID_HH
#define FRACKIT_PYTHON_COMMON_ID_HH

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <frackit/common/id.hh>

namespace Frackit::Python {

namespace py = pybind11;

void registerId(py::module& module)
{
    py::class_<Id> cls(module, "Id");
    cls.def(py::init<std::size_t>());
    cls.def("get", &Id::get, "return the id");
    cls.def(py::self == Id());
    cls.def(hash(py::self));
}

} // end namespace Frackit::Python

#endif
