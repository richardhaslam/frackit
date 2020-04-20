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
#ifndef FRACKIT_PYTHON_GEOMETRY_ELLIPSE_HH
#define FRACKIT_PYTHON_GEOMETRY_ELLIPSE_HH

#include <string>

#include <pybind11/pybind11.h>

#include <frackit/geometry/geometry.hh>
#include <frackit/geometry/ellipticalgeometry.hh>
#include <frackit/geometry/ellipse.hh>
#include "registerdimensionproperties.hh"

namespace Frackit::Python {

namespace py = pybind11;

namespace Detail {

    template<class ctype, int worldDim>
    void registerEllipse(py::module& module)
    {
        using namespace py::literals;
        using Ellipse = Ellipse<ctype, worldDim>;
        using EllipticalGeometry = EllipticalGeometry<ctype, worldDim>;
        using Point = typename Ellipse::Point;
        using Direction = typename Ellipse::Direction;

        // register class and define constructor
        std::string className("Ellipse_" + std::to_string(worldDim));
        py::class_<Ellipse, Geometry, EllipticalGeometry> cls(module, className.c_str());
        cls.def(py::init<const Point&, const Direction&, const Direction&, ctype, ctype>(),
                "center"_a, "majorAxis"_a, "minorAxis"_a, "majorAxisLength"_a, "minorAxisLength"_a);

        // dimensionality properties
        registerDimensionProperties(cls);

        // member functions
        cls.def("name", &Ellipse::name, "name of the geometry");

        // contains queries
        cls.def("contains", py::overload_cast<const Point&, bool>(&Ellipse::contains, py::const_),
                "point"_a, "checkIfOnPlane"_a = true,
                "returns true if the given point is on the ellipse (default tolerance)");
        cls.def("contains", py::overload_cast<const Point&, ctype, bool>(&Ellipse::contains, py::const_),
                "point"_a, "eps"_a, "checkIfOnPlane"_a = true,
                "returns true if the given point is on the ellipse (given tolerance)");

        // get a point from local coordinate or angle
        cls.def("getPoint", &Ellipse::getPoint, "localCoordinate"_a,
                "return the point on the ellipse for the given local coordinate (0.0 <= localCoordinate <= 1)");
        cls.def("getPointFromAngle", &Ellipse::getPointFromAngle, "angle"_a,
                "return the point on the ellipse for the given angle w.r.t. the ellipse-local coordinate system (angle >= 0.0)");
        cls.def("getParam", &Ellipse::getParam, "point"_a, "checkIfOnEllipse"_a = true,
                "returns the local coordinate of the given point on the ellipse");
        cls.def("getAngle", &Ellipse::getAngle, "point"_a, "checkIfOnEllipse"_a = true,
                "returns the angle w.r.t. the ellipse-local coordinate system of the given point");
    }

} // end namespace detail

template<class ctype>
void registerEllipse(py::module& module)
{
    Detail::registerEllipse<ctype, 3>(module);
}

} // end namespace Frackit::Python

#endif
