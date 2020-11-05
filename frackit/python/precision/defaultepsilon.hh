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
#ifndef FRACKIT_PYTHON_DEFAULT_EPSILON_HH
#define FRACKIT_PYTHON_DEFAULT_EPSILON_HH

#include <pybind11/pybind11.h>

// supported geometry types (so far)
#include <frackit/geometry/segment.hh>
#include <frackit/geometry/circle.hh>
#include <frackit/geometry/ellipse.hh>
#include <frackit/geometry/disk.hh>
#include <frackit/geometry/quadrilateral.hh>
#include <frackit/geometry/cylindersurface.hh>
#include <frackit/geometry/box.hh>
#include <frackit/python/geometry/brepwrapper.hh>
#include <frackit/python/geometry/ctype.hh>

#include <frackit/precision/defaultepsilon.hh>

namespace Frackit::Python {

namespace py = pybind11;

namespace Detail {

    template<class Geo>
    typename CoordinateTypeTraits<Geo>::type
    defaultEpsilon(const Geo& geo)
    { return Frackit::defaultEpsilon(getUnwrappedShape(geo)); }

    template<class Geo>
    void registerDefaultEpsilon(py::module& module)
    {
        module.def("defaultEpsilon",
                   py::overload_cast<const Geo&>(&defaultEpsilon<Geo>),
                   "returns a default tolerance value for the given geometry");
    }

} // end namespace Detail

template<class ctype>
void registerDefaultEpsilon(py::module& module)
{
    Detail::registerDefaultEpsilon<Segment<ctype, 3>>(module);
    Detail::registerDefaultEpsilon<Circle<ctype, 3>>(module);
    Detail::registerDefaultEpsilon<Ellipse<ctype, 3>>(module);
    Detail::registerDefaultEpsilon<Disk<ctype>>(module);
    Detail::registerDefaultEpsilon<Quadrilateral<ctype, 3>>(module);
    Detail::registerDefaultEpsilon<CylinderSurface<ctype>>(module);
    Detail::registerDefaultEpsilon<Box<ctype>>(module);
    Detail::registerDefaultEpsilon<ShapeWrapper>(module);
}

} // end namespace Frackit::Python

#endif // FRACKIT_PYTHON_DEFAULT_EPSILON_HH
