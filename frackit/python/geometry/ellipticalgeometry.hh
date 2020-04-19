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
#ifndef FRACKIT_PYTHON_GEOMETRY_ELLIPTICALGEOMETRY_HH
#define FRACKIT_PYTHON_GEOMETRY_ELLIPTICALGEOMETRY_HH

#include <string>

#include <pybind11/pybind11.h>
#include <frackit/geometry/ellipticalgeometry.hh>

namespace Frackit::Python {

namespace py = pybind11;

namespace Detail {

    template<class ctype, int worldDim>
    void registerEllipticalGeometry(py::module& module)
    {
        using namespace py::literals;
        using EllipticalGeometry = EllipticalGeometry<ctype, worldDim>;
        using Point = typename EllipticalGeometry::Point;
        using Direction = typename EllipticalGeometry::Direction;

        // register class and define constructorss
        std::string className("EllipticalGeometry_" + std::to_string(worldDim));
        py::class_<EllipticalGeometry> cls(module, className.c_str());
        cls.def(py::init<>());
        cls.def(py::init<const Point&, const Direction&, const Direction&, ctype, ctype>(),
                "center"_a, "majorAxis"_a, "minorAxis"_a, "majorAxisLength"_a, "minorAxisLength"_a);

        // member functions
        cls.def("majorAxis", &EllipticalGeometry::majorAxis, py::return_value_policy::reference_internal, "returns the major axis");
        cls.def("minorAxis", &EllipticalGeometry::minorAxis, py::return_value_policy::reference_internal, "returns the minor axis");
        cls.def("normal", &EllipticalGeometry::normal, py::return_value_policy::reference_internal, "returns the normal of the supporting plane");
        cls.def("center", &EllipticalGeometry::center, py::return_value_policy::reference_internal, "returns the center point");

        cls.def("majorAxisLength", &EllipticalGeometry::majorAxisLength, "returns the major axis length");
        cls.def("minorAxisLength", &EllipticalGeometry::minorAxisLength, "returns the minor axis length");
        cls.def("supportingPlane", &EllipticalGeometry::supportingPlane, "returns the supporting plane");
    }

} // end namespace detail

template<class ctype>
void registerEllipticalGeometry(py::module& module)
{
    Detail::registerEllipticalGeometry<ctype, 3>(module);
}

} // end namespace Frackit::Python

#endif
