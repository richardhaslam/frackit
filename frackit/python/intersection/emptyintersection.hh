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
#ifndef FRACKIT_PYTHON_EMPTY_INTERSECTION_HH
#define FRACKIT_PYTHON_EMPTY_INTERSECTION_HH

#include <string>
#include <pybind11/pybind11.h>

#include <frackit/intersection/emptyintersection.hh>
#include <frackit/python/geometry/registerdimensionproperties.hh>

namespace Frackit::Python {

namespace py = pybind11;

namespace Detail {

template<class ct, int wd>
void registerEmptyIntersection(py::module& module)
{
    using EI = Frackit::EmptyIntersection<wd, ct>;
    const std::string className = "EmptyIntersection_" + std::to_string(wd);

    py::class_<EI> cls(module, className.c_str());
    cls.def("name", &EI::name, "return the name of this (empty) geometry");
    registerDimensionProperties(cls);
}

} // end namespace Detail

template<class ctype>
void registerEmptyIntersection(py::module& module)
{
    Detail::registerEmptyIntersection<ctype, 1>(module);
    Detail::registerEmptyIntersection<ctype, 2>(module);
    Detail::registerEmptyIntersection<ctype, 3>(module);
}

} // end namespace Frackit::Python

#endif
