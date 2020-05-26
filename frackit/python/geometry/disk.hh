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
#ifndef FRACKIT_PYTHON_GEOMETRY_DISK_HH
#define FRACKIT_PYTHON_GEOMETRY_DISK_HH

#include <pybind11/pybind11.h>

#include <frackit/geometry/geometry.hh>
#include <frackit/geometry/ellipticalgeometry.hh>
#include <frackit/geometry/disk.hh>
#include "registerdimensionproperties.hh"

namespace Frackit::Python {

namespace py = pybind11;

template<class ctype>
void registerDisk(py::module& module)
{
    using namespace py::literals;
    using Disk = Disk<ctype>;
    using Point = typename Disk::Point;
    using Ellipse = typename Disk::Ellipse;
    using EllipticalGeometry = EllipticalGeometry<ctype, 3>;

    // register class and define constructor
    py::class_<Disk, Geometry, EllipticalGeometry> cls(module, "Disk");
    cls.def(py::init<const Ellipse&>(), "ellipse"_a);

    // dimensionality properties
    registerDimensionProperties(cls);

    // member functions
    cls.def("name", &Disk::name, "name of the geometry");
    cls.def("area", &Disk::area, "area of the disk");
    cls.def("boundingEllipse", &Disk::boundingEllipse, "returns the ellipse that describes the disk boundary");

    // contains queries
    cls.def("contains", py::overload_cast<const Point&, bool>(&Disk::contains, py::const_),
            "point"_a, "checkIfOnEllipse"_a = true,
            "returns true if the given point is on the disk (default tolerance)");
    cls.def("contains", py::overload_cast<const Point&, ctype, bool>(&Disk::contains, py::const_),
            "point"_a, "eps"_a, "checkIfOnEllipse"_a = true,
            "returns true if the given point is on the disk (given tolerance)");

    // get a point from local coordinate or angle
    cls.def("getPoint", &Disk::getPoint, "angularFraction"_a, "radialFraction"_a,
            "return the point on the ellipse arc for the given angular and radial fractions (0.0 <= angular/radial fraction <= 1)");

    cls.def("__repr__", [&] (const Disk& d) { return "Frackit::Disk"; });
}

} // end namespace Frackit::Python

#endif
