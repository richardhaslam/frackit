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
#ifndef FRACKIT_PYTHON_GEOMETRY_GET_BOUNDING_BOX_HH
#define FRACKIT_PYTHON_GEOMETRY_GET_BOUNDING_BOX_HH

#include <string>

#include <pybind11/pybind11.h>

#include <frackit/geometry/box.hh>
#include <frackit/geometry/circle.hh>
#include <frackit/geometry/cylinder.hh>
#include <frackit/geometry/cylindersurface.hh>
#include <frackit/geometry/disk.hh>
#include <frackit/geometry/ellipse.hh>
#include <frackit/geometry/ellipsearc.hh>
#include <frackit/geometry/polygon.hh>
#include <frackit/geometry/quadrilateral.hh>
#include <frackit/geometry/segment.hh>
#include <frackit/geometry/sphere.hh>

#include <frackit/python/geometry/ctype.hh>
#include <frackit/python/geometry/brepwrapper.hh>
#include <frackit/geometryutilities/getboundingbox.hh>

namespace Frackit::Python {

template<class Geo>
Frackit::Box<typename CoordinateTypeTraits<Geo>::type>
getBoundingBox(const Geo& geo)
{ return Frackit::getBoundingBox( getUnwrappedShape(geo) ); }

namespace Detail {

template<class Geometry>
void registerGetBoundingBox(py::module& module, const std::string& docFinish)
{
    const std::string doc = "Returns the bounding box of the given " + docFinish;
    module.def("getBoundingBox",
               py::overload_cast<const Geometry&>(&Frackit::Python::getBoundingBox<Geometry>),
               doc.c_str());

}

} // end namespace Detail

template<class ctype>
void registerGetBoundingBox(py::module& module)
{
    Detail::registerGetBoundingBox<Box<ctype>>(module, "box.");
    Detail::registerGetBoundingBox<Cylinder<ctype>>(module, "cylinder.");
    Detail::registerGetBoundingBox<CylinderSurface<ctype>>(module, "cylinder surface.");
    Detail::registerGetBoundingBox<Disk<ctype>>(module, "disk.");
    Detail::registerGetBoundingBox<Sphere<ctype>>(module, "sphere.");

    Detail::registerGetBoundingBox<Circle<ctype, 3>>(module, "circle in 3d space.");
    Detail::registerGetBoundingBox<Ellipse<ctype, 3>>(module, "ellipse in 3d space.");
    Detail::registerGetBoundingBox<EllipseArc<ctype, 3>>(module, "ellipse arc in 3d space.");
    Detail::registerGetBoundingBox<Polygon<ctype, 3>>(module, "polygon in 3d space.");
    Detail::registerGetBoundingBox<Quadrilateral<ctype, 3>>(module, "quadrilateral in 3d space.");
    Detail::registerGetBoundingBox<Segment<ctype, 3>>(module, "segment in 3d space.");

    Detail::registerGetBoundingBox<ShapeWrapper>(module, "(wrapped) OpenCascade shape.");
    Detail::registerGetBoundingBox<VertexWrapper>(module, "(wrapped) OpenCascade vertex.");
    Detail::registerGetBoundingBox<EdgeWrapper>(module, "(wrapped) OpenCascade edge.");
    Detail::registerGetBoundingBox<WireWrapper>(module, "(wrapped) OpenCascade wire.");
    Detail::registerGetBoundingBox<FaceWrapper>(module, "(wrapped) OpenCascade face.");
    Detail::registerGetBoundingBox<ShellWrapper>(module, "(wrapped) OpenCascade shell.");
    Detail::registerGetBoundingBox<SolidWrapper>(module, "(wrapped) OpenCascade solid.");
    Detail::registerGetBoundingBox<CompoundWrapper>(module, "(wrapped) OpenCascade compound.");
}

} // end namespace Frackit::Python

#endif // FRACKIT_PYTHON_GEOMETRY_GET_BOUNDING_BOX_HH
