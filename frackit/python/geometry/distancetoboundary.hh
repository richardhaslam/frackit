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
#ifndef FRACKIT_PYTHON_GEOMETRY_DISTANCE_TO_BOUNDARY_HH
#define FRACKIT_PYTHON_GEOMETRY_DISTANCE_TO_BOUNDARY_HH

#include <pybind11/pybind11.h>

#include <frackit/geometry/point.hh>
#include <frackit/geometry/disk.hh>
#include <frackit/geometry/quadrilateral.hh>
#include <frackit/geometry/polygon.hh>
#include <frackit/geometry/cylindersurface.hh>
#include <frackit/geometry/box.hh>
#include <frackit/python/geometry/brepwrapper.hh>

#include <frackit/distance/distancetoboundary.hh>

namespace Frackit::Python {

namespace py = pybind11;

namespace Detail {

template<class Geo1, class Geo2, class ctype>
ctype computeDistanceToBoundary(const Geo1& geo1, const Geo2& geo2)
{ return Frackit::computeDistance(getUnwrappedShape(geo1), getUnwrappedShape(geo2)); }

} // end namespace Detail

template<class ctype>
void registerComputeDistanceToBoundary(py::module& module)
{
    using namespace py::literals;
    using Point = Frackit::Point<ctype, 3>;
    using Quad = Frackit::Quadrilateral<ctype, 3>;
    using Polygon = Frackit::Polygon<ctype, 3>;
    using Disk = Frackit::Disk<ctype>;
    using CylinderSurface = Frackit::CylinderSurface<ctype>;
    using Box = Frackit::Box<ctype>;
    using ShapeWrapper = Frackit::Python::ShapeWrapper;
    using FaceWrapper = Frackit::Python::FaceWrapper;
    using SolidWrapper = Frackit::Python::SolidWrapper;

    // explicitly implemented overloads
    module.def("computeDistanceToBoundary",
               py::overload_cast<const Point&, const Quad&>(&Detail::computeDistanceToBoundary<Point, Quad, ctype>),
               "Returns the minimum euclidian distance of a point to the boundary of a quadrilateral");

    module.def("computeDistanceToBoundary",
               py::overload_cast<const Point&, const Polygon&>(&Detail::computeDistanceToBoundary<Point, Polygon, ctype>),
               "Returns the minimum euclidian distance of a point to the boundary of a polygon");

    module.def("computeDistanceToBoundary",
               py::overload_cast<const Point&, const Box&>(&Detail::computeDistanceToBoundary<Point, Box, ctype>),
               "Returns the minimum euclidian distance of a point to the surface of a box");

    // register overloads with shapes
    module.def("computeDistanceToBoundary",
               py::overload_cast<const ShapeWrapper&, const Disk&>(&Detail::computeDistanceToBoundary<ShapeWrapper, Disk, ctype>),
               "Returns the minimum euclidian distance of a shape to the bounding ellipse of a disk");

    module.def("computeDistanceToBoundary",
               py::overload_cast<const ShapeWrapper&, const CylinderSurface&>(&Detail::computeDistanceToBoundary<ShapeWrapper, CylinderSurface, ctype>),
               "Returns the minimum euclidian distance of a shape to the bounding circles of a cylinder surface");

    module.def("computeDistanceToBoundary",
               py::overload_cast<const ShapeWrapper&, const FaceWrapper&>(&Detail::computeDistanceToBoundary<ShapeWrapper, FaceWrapper, ctype>),
               "Returns the minimum euclidian distance of a shape to the bounding wire of a face shape");

    module.def("computeDistanceToBoundary",
               py::overload_cast<const ShapeWrapper&, const Quad&>(&Detail::computeDistanceToBoundary<ShapeWrapper, Quad, ctype>),
               "Returns the minimum euclidian distance of a shape to the bounding edges of a quadrilateral");

    module.def("computeDistanceToBoundary",
               py::overload_cast<const ShapeWrapper&, const Polygon&>(&Detail::computeDistanceToBoundary<ShapeWrapper, Polygon, ctype>),
               "Returns the minimum euclidian distance of a shape to the bounding edges of a polygon");

    module.def("computeDistanceToBoundary",
               py::overload_cast<const ShapeWrapper&, const SolidWrapper&>(&Detail::computeDistanceToBoundary<ShapeWrapper, SolidWrapper, ctype>),
               "Returns the minimum euclidian distance of a shape to the bounding shell of a solid");
}

} // end namespace Frackit::Python

#endif
