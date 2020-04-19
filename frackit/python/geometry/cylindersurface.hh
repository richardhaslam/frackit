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
#ifndef FRACKIT_PYTHON_GEOMETRY_CYLINDER_SURFACE_HH
#define FRACKIT_PYTHON_GEOMETRY_CYLINDER_SURFACE_HH

#include <pybind11/pybind11.h>
#include <frackit/geometry/geometry.hh>
#include <frackit/geometry/cylindersurface.hh>

namespace Frackit::Python {

namespace py = pybind11;

template<class ctype>
void registerCylinderSurface(py::module& module)
{
    using namespace py::literals;
    using CylinderSurface = CylinderSurface<ctype>;
    using Circle = typename CylinderSurface::Circle;
    using Point = typename CylinderSurface::Point;

    // register class and define constructors
    py::class_<CylinderSurface, Geometry> cls(module, "CylinderSurface");
    cls.def(py::init<ctype, ctype>(), "radius"_a, "height"_a);
    cls.def(py::init<const Circle&, ctype>(), "bottomCircle"_a, "height"_a);

    // getter functions
    cls.def("name", &CylinderSurface::name, "name of the geometry class");
    cls.def("base1", &CylinderSurface::base1, "first basis vector orthogonal to center line");
    cls.def("base2", &CylinderSurface::base2, "second basis vector orthogonal to center line");
    cls.def("base3", &CylinderSurface::base3, "basis vector parallel to center line");
    cls.def("direction", &CylinderSurface::direction, "direction of the cylinder axis");

    cls.def("upperBoundingCircle", &CylinderSurface::upperBoundingCircle, "returns the upper bounding circle of the cylinder surface");
    cls.def("lowerBoundingCircle", &CylinderSurface::lowerBoundingCircle, "returns the lower bounding circle of the cylinder surface");
    cls.def("cylinder", &CylinderSurface::cylinder, "returns the cylinder this is the lateral surface of");
    cls.def("centerSegment", &CylinderSurface::centerSegment, "returns the segment describing the cylinder axis");
    cls.def("getTangentPlane", &CylinderSurface::getTangentPlane, "point"_a,
            "returns the tangent plane in the provided point p");

    cls.def("height", &CylinderSurface::height, "returns the height of the the cylinder  surface");
    cls.def("radius", &CylinderSurface::radius, "returns the raidus of the cylinder surface");
    cls.def("area", &CylinderSurface::area, "returns the area of the cylinder surface");

    // contains queries
    cls.def("contains", py::overload_cast<const Point&>(&CylinderSurface::contains, py::const_),
            "point"_a, "returns true if the given point is on the cylinder surface (default tolerance)");
    cls.def("contains", py::overload_cast<const Point&, ctype>(&CylinderSurface::contains, py::const_),
            "point"_a, "eps"_a, "returns true if the given point is on the cylinder surface (given tolerance)");
}

} // end namespace Frackit::Python

#endif
