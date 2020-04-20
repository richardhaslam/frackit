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
#ifndef FRACKIT_PYTHON_GEOMETRY_ELLIPSEARC_HH
#define FRACKIT_PYTHON_GEOMETRY_ELLIPSEARC_HH

#include <string>

#include <pybind11/pybind11.h>
#include <frackit/geometry/ellipsearc.hh>
#include "registerdimensionproperties.hh"

namespace Frackit::Python {

namespace py = pybind11;

namespace Detail {

    template<class ctype, int worldDim>
    void registerEllipseArc(py::module& module)
    {
        using namespace py::literals;
        using EllipseArc = EllipseArc<ctype, worldDim>;
        using Ellipse = typename EllipseArc::Ellipse;
        using Point = typename Ellipse::Point;

        // register class and define constructor
        std::string className("EllipseArc_" + std::to_string(worldDim));
        py::class_<EllipseArc, Ellipse> cls(module, className.c_str());
        cls.def(py::init<const Ellipse&, const Point&, const Point&>(),
                "supportingEllipse"_a, "source"_a, "target"_a);

        // dimensionality properties
        registerDimensionProperties(cls);

        // member functions
        cls.def("name", &EllipseArc::name, "name of the geometry");
        cls.def("source", &EllipseArc::source, "source point of the arc");
        cls.def("target", &EllipseArc::target, "target point of the arc");
        cls.def("supportingEllipse", &EllipseArc::supportingEllipse, "returns the supporting ellipse");
        cls.def("isFullEllipse", &EllipseArc::isFullEllipse, "returns true if the arc describes a full ellipse");

        // angles of source/target points
        cls.def("sourceAngleOnEllipse", &EllipseArc::sourceAngleOnEllipse,
                "returns the angle of the source point w.r.t. to the ellipse-local coordinate system");
        cls.def("targetAngleOnEllipse", &EllipseArc::targetAngleOnEllipse,
                "returns the angle of the target point w.r.t. to the ellipse-local coordinate system");

        // contains queries
        cls.def("contains", py::overload_cast<const Point&, bool>(&EllipseArc::contains, py::const_),
                "point"_a, "checkIfOnEllipse"_a = true,
                "returns true if the given point is on the ellipse arc (default tolerance)");
        cls.def("contains", py::overload_cast<const Point&, ctype, bool>(&EllipseArc::contains, py::const_),
                "point"_a, "eps"_a, "checkIfOnEllipse"_a = true,
                "returns true if the given point is on the ellipse arc (given tolerance)");

        // get a point from local coordinate or angle
        cls.def("getPoint", &EllipseArc::getPoint, "localCoordinate"_a,
                "return the point on the ellipse arc for the given local coordinate (0.0 <= localCoordinate <= 1)");
    }

} // end namespace detail

template<class ctype>
void registerEllipseArc(py::module& module)
{
    Detail::registerEllipseArc<ctype, 3>(module);
}

} // end namespace Frackit::Python

#endif
