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
#ifndef FRACKIT_PYTHON_GEOMETRY_MAGNITUDE_HH
#define FRACKIT_PYTHON_GEOMETRY_MAGNITUDE_HH

#include <pybind11/pybind11.h>


// supported (registered) geometry types (so far)
#include <frackit/geometry/disk.hh>
#include <frackit/geometry/segment.hh>
#include <frackit/geometry/vector.hh>
#include <frackit/geometry/circle.hh>
#include <frackit/geometry/ellipse.hh>
#include <frackit/geometry/ellipsearc.hh>
#include <frackit/geometry/triangle.hh>
#include <frackit/geometry/quadrilateral.hh>
#include <frackit/geometry/polygon.hh>
#include <frackit/geometry/box.hh>
#include <frackit/geometry/sphere.hh>
#include <frackit/geometry/cylinder.hh>
#include <frackit/geometry/hollowcylinder.hh>
#include <frackit/geometry/cylindersurface.hh>
#include <frackit/python/geometry/brepwrapper.hh>
#include <frackit/python/geometry/ctype.hh>

#include <frackit/magnitude/magnitude.hh>
#include <frackit/magnitude/containedmagnitude.hh>

namespace Frackit::Python {

namespace py = pybind11;

namespace Detail {

    //! overload for computation of the magnitude of the contained part of a shape wrapper
    template<class Geo, class Domain>
    typename CoordinateTypeTraits<Geo>::type
    computeContainedMagnitude(const Geo& geo, const Domain& domain)
    { return Frackit::computeContainedMagnitude(getUnwrappedShape(geo), getUnwrappedShape(domain)); }

    //! overload for computation of the magnitude of the contained part of a shape wrapper
    template<class Geo>
    typename CoordinateTypeTraits<Geo>::type
    computeMagnitude(const Geo& geo)
    { return Frackit::computeMagnitude(getUnwrappedShape(geo)); }

    template<class Geo>
    void registerMagnitude(py::module& module)
    {
        module.def("computeMagnitude",
                   py::overload_cast<const Geo&>(&computeMagnitude<Geo>),
                   "compute the magnitude (length/area/volume) of a geometry");
    }

    template<class Geo, class Domain>
    void registerContainedMagnitude(py::module& module)
    {
        module.def("computeContainedMagnitude",
                   py::overload_cast<const Geo&, const Domain&>(&computeContainedMagnitude<Geo, Domain>),
                   "compute the magnitude of the part of an entity geometry contained in a domain geometry");
    }

} // end namespace Detail

template<class ctype>
void registerMagnitude(py::module& module)
{
    using namespace Frackit;
    Detail::registerMagnitude< Vector<ctype, 3> >(module);
    Detail::registerMagnitude< Segment<ctype, 3> >(module);
    Detail::registerMagnitude< Circle<ctype, 3> >(module);
    Detail::registerMagnitude< Ellipse<ctype, 3> >(module);
    Detail::registerMagnitude< EllipseArc<ctype, 3> >(module);
    Detail::registerMagnitude< Disk<ctype> >(module);
    Detail::registerMagnitude< Triangle<ctype, 3> >(module);
    Detail::registerMagnitude< Quadrilateral<ctype, 3> >(module);
    Detail::registerMagnitude< Polygon<ctype, 3> >(module);
    Detail::registerMagnitude< CylinderSurface<ctype> >(module);
    Detail::registerMagnitude< Cylinder<ctype> >(module);
    Detail::registerMagnitude< HollowCylinder<ctype> >(module);
    Detail::registerMagnitude< Box<ctype> >(module);
    Detail::registerMagnitude< Sphere<ctype> >(module);
    Detail::registerMagnitude< EdgeWrapper >(module);
    Detail::registerMagnitude< WireWrapper >(module);
    Detail::registerMagnitude< FaceWrapper >(module);
    Detail::registerMagnitude< ShellWrapper >(module);
    Detail::registerMagnitude< SolidWrapper >(module);
}

template<class ctype>
void registerContainedMagnitude(py::module& module)
{
    using namespace Frackit;
    Detail::registerContainedMagnitude< Segment<ctype, 3>, Box<ctype> >(module);
    Detail::registerContainedMagnitude< Segment<ctype, 3>, Cylinder<ctype> >(module);
    Detail::registerContainedMagnitude< Segment<ctype, 3>, SolidWrapper >(module);

    Detail::registerContainedMagnitude< Disk<ctype>, Box<ctype> >(module);
    Detail::registerContainedMagnitude< Disk<ctype>, Cylinder<ctype> >(module);
    Detail::registerContainedMagnitude< Disk<ctype>, SolidWrapper >(module);

    Detail::registerContainedMagnitude< Quadrilateral<ctype, 3>, Box<ctype> >(module);
    Detail::registerContainedMagnitude< Quadrilateral<ctype, 3>, Cylinder<ctype> >(module);
    Detail::registerContainedMagnitude< Quadrilateral<ctype, 3>, SolidWrapper >(module);

    Detail::registerContainedMagnitude< FaceWrapper, Box<ctype> >(module);
    Detail::registerContainedMagnitude< FaceWrapper, Cylinder<ctype> >(module);
    Detail::registerContainedMagnitude< FaceWrapper, SolidWrapper >(module);
}

} // end namespace Frackit::Python

#endif
