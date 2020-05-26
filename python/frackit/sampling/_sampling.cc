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
#include <pybind11/pybind11.h>

#include <frackit/geometry/box.hh>
#include <frackit/geometry/cylinder.hh>
#include <frackit/geometry/hollowcylinder.hh>
#include <frackit/sampling/makeuniformpointsampler.hh>

#include <frackit/python/sampling/status.hh>
#include <frackit/python/sampling/pointsampler.hh>
#include <frackit/python/sampling/boxpointsampler.hh>
#include <frackit/python/sampling/cylinderpointsampler.hh>

namespace Frackit::Python {

template<class ctype>
auto makeUniformBoxPointSampler(const Frackit::Box<ctype>& box)
-> decltype(Frackit::makeUniformPointSampler(box))
{ return Frackit::makeUniformPointSampler(box); }

template<class ctype>
auto makeUniformCylinderPointSampler(const Frackit::Cylinder<ctype>& cylinder)
-> decltype(Frackit::makeUniformPointSampler(cylinder))
{ return Frackit::makeUniformPointSampler(cylinder); }

template<class ctype>
auto makeUniformHollowCylinderPointSampler(const Frackit::HollowCylinder<ctype>& hollowCyl)
-> decltype(Frackit::makeUniformPointSampler(hollowCyl))
{ return Frackit::makeUniformPointSampler(hollowCyl); }

template<class ctype>
void registerPointSamplerCreatorFunctions(pybind11::module& module)
{
    using namespace pybind11::literals;
    module.def("makeUniformPointSampler", &makeUniformBoxPointSampler<ctype>,
               "box"_a, "returns a uniform point sampler on the given box");
    module.def("makeUniformPointSampler", &makeUniformCylinderPointSampler<ctype>,
               "cylinder"_a, "returns a uniform point sampler on the given cylinder");
    module.def("makeUniformPointSampler", &makeUniformHollowCylinderPointSampler<ctype>,
               "hollowcylinder"_a, "returns a uniform point sampler on the given hollow cylinder");
}

} // end namespace Frackit::Python

PYBIND11_MODULE(_sampling, module)
{
    // register status class
    Frackit::Python::registerSamplingStatus(module);

    // register base class for all dimensions
    Frackit::Python::registerPointSamplerInterface<double, 1>(module);
    Frackit::Python::registerPointSamplerInterface<double, 2>(module);
    Frackit::Python::registerPointSamplerInterface<double, 3>(module);

    // register point samplers on geometries
    Frackit::Python::registerBoxPointSampler<double>(module);
    Frackit::Python::registerCylinderPointSampler<double>(module);

    // register creator functions
    Frackit::Python::registerPointSamplerCreatorFunctions<double>(module);
}
