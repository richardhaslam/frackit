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
#ifndef FRACKIT_PYTHON_PRECISION_HH
#define FRACKIT_PYTHON_PRECISION_HH

#include <pybind11/pybind11.h>
#include <frackit/precision/precision.hh>

namespace Frackit::Python {

template<class ctype>
void registerPrecision(pybind11::module& module)
{
    pybind11::class_<Precision<ctype>> cls(module, "Precision");
    cls.def_property_readonly_static("confusion", [] (pybind11::object /*self*/) { return Precision<ctype>::confusion(); }, "tolerance value for distance comparisons");
    cls.def_property_readonly_static("angular", [] (pybind11::object /*self*/) { return Precision<ctype>::angular(); }, "tolerance value for angle comparisons");
}

} // end namespace Frackit::Python

#endif // FRACKIT_PYTHON_PRECISION_HH
