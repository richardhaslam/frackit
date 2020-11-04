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
#ifndef FRACKIT_PYTHON_GEOMETRY_DISTANCE_HH
#define FRACKIT_PYTHON_GEOMETRY_DISTANCE_HH

#include <pybind11/pybind11.h>

#include <frackit/geometry/point.hh>
#include <frackit/geometry/segment.hh>

#include <frackit/python/occutilities/brepwrapper.hh>
#include <frackit/distance/distance.hh>

namespace Frackit::Python {

namespace py = pybind11;

namespace Detail {

template<class Geo1, class Geo2, class ctype>
ctype computeDistance(const Geo1& geo1, const Geo2& geo2)
{
    return Frackit::computeDistance(OCCUtilities::getUnwrappedShape(geo1),
                                    OCCUtilities::getUnwrappedShape(geo2));
}

} // end namespace Detail

template<class ctype>
void registerComputeDistance(py::module& module)
{
    using namespace py::literals;
    using Point = Frackit::Point<ctype, 3>;
    using Segment = Frackit::Segment<ctype, 3>;
    using ShapeWrapper = OCCUtilities::ShapeWrapper;

    module.def("computeDistance",
               py::overload_cast<const Point&, const Point&>(&Detail::computeDistance<Point, Point, ctype>),
               "Returns the minimum euclidian distance between two points in 3d");

    module.def("computeDistance",
               py::overload_cast<const Segment&, const Point&>(&Detail::computeDistance<Segment, Point, ctype>),
               "Returns the minimum euclidian distance between a segment and a point in 3d",
               "segment"_a, "point"_a);

    module.def("computeDistance",
               py::overload_cast<const Point&, const Segment&>(&Detail::computeDistance<Point, Segment, ctype>),
               "Returns the minimum euclidian distance between a point and a segment in 3d",
               "point"_a, "segment"_a);

    module.def("computeDistance",
               py::overload_cast<const ShapeWrapper&, const ShapeWrapper&>(&Detail::computeDistance<ShapeWrapper, ShapeWrapper, ctype>),
               "Returns the minimum euclidian distance between two generic shapes");
}

} // end namespace Frackit::Python

#endif
