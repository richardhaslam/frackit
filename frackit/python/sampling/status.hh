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
#ifndef FRACKIT_PYTHON_SAMPLING_STATUS_HH
#define FRACKIT_PYTHON_SAMPLING_STATUS_HH

#include <pybind11/pybind11.h>

#include <frackit/common/id.hh>
#include <frackit/sampling/status.hh>

namespace Frackit::Python {

namespace py = pybind11;

void registerSamplingStatus(py::module& module)
{
    py::class_<SamplingStatus> cls(module, "SamplingStatus");
    cls.def(py::init());
    cls.def("setTargetCount",
            &SamplingStatus::setTargetCount,
            "Define the target number of samples for the given id");
    cls.def("finished",
            py::overload_cast<>(&SamplingStatus::finished),
            "Returns true if overall target sample count is reached");
    cls.def("finished",
            py::overload_cast<const Id&>(&SamplingStatus::finished),
            "Returns true if the sample count for the given id is reached");
    cls.def("increaseCounter",
            &SamplingStatus::increaseCounter,
            "Register that a new sample for the given id was taken");
    cls.def("resetCounters",
            &SamplingStatus::resetCounters,
            "Reset all counters");
    cls.def("resetCounter",
            &SamplingStatus::resetCounter,
            "Reset the counter for the given id");
    cls.def("reset",
            &SamplingStatus::reset,
            "Reset everything");
    cls.def("increaseRejectedCounter",
            &SamplingStatus::increaseRejectedCounter,
            "Register that a sample has been rejected");

    cls.def("getCount",
            py::overload_cast<>(&SamplingStatus::getCount),
            "Returns the overall number of entities");
    cls.def("getCount",
            py::overload_cast<const Id&>(&SamplingStatus::getCount),
            "Returns the number of entities for the given id");

    using namespace py::literals;
    cls.def("print",
            &SamplingStatus::print,
            "forceHeaderPrint"_a=false,
            "Print the current status");
}

} // end namespace Frackit::Python

#endif
