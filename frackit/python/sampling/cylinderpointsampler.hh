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
#ifndef FRACKIT_PYTHON_SAMPLING_CYLINDER_POINT_SAMPLER_HH
#define FRACKIT_PYTHON_SAMPLING_CYLINDER_POINT_SAMPLER_HH

#include <string>

#include <pybind11/pybind11.h>

#include <frackit/geometry/cylinder.hh>
#include <frackit/sampling/uniformpointsamplertraits.hh>
#include <frackit/sampling/cylinderpointsampler.hh>
#include <frackit/sampling/pointsampler.hh>

namespace Frackit::Python {

namespace py = pybind11;

template<class ctype>
void registerCylinderPointSampler(py::module& module)
{
    // TODO: Once JIT compilation is available, rewrite this
    //       such that different distributions can be used
    using Traits = Frackit::UniformPointSamplerTraits<ctype, 3>;
    using UniformSampler = Frackit::CylinderPointSampler<ctype, Traits>;

    using Point = typename UniformSampler::Point;
    using SamplerInterface = PointSampler<Point>;

    py::class_<UniformSampler, SamplerInterface> cls(module, "UniformCylinderPointSampler");

    using namespace py::literals;
    using Cylinder = Cylinder<ctype>;
    using DistroBase1 = typename Traits::DistributionBase1;
    using DistroBase2 = typename Traits::DistributionBase2;
    using DistroBase3 = typename Traits::DistributionBase3;

    // register constructor
    cls.def(py::init<const Cylinder&, const DistroBase1&, const DistroBase2&, const DistroBase3&>(),
            "cylinder"_a, "distributionR"_a, "distributionPhi"_a, "distributionH"_a);

    // register point sample function
    cls.def("__call__", &UniformSampler::operator(), "returns a randomly sampled point");
}

} // end namespace Frackit::Python

#endif
