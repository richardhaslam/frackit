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
#ifndef FRACKIT_PYTHON_BREP_UTILITIES_HH
#define FRACKIT_PYTHON_BREP_UTILITIES_HH

#include <string>
#include <type_traits>

#include <BRepTools.hxx>

#include <frackit/geometry/point.hh>
#include <frackit/geometry/segment.hh>
#include <frackit/geometry/circle.hh>
#include <frackit/geometry/ellipse.hh>
#include <frackit/geometry/ellipsearc.hh>
#include <frackit/geometry/quadrilateral.hh>
#include <frackit/geometry/disk.hh>
#include <frackit/geometry/box.hh>
#include <frackit/geometry/cylinder.hh>
#include <frackit/geometry/cylindersurface.hh>

#include <frackit/occ/breputilities.hh>
#include <frackit/occ/gputilities.hh>
#include <frackit/occ/geomutilities.hh>

#include "brepwrapper.hh"

namespace Frackit::Python {
namespace OCCUtilities {

    template<class Geometry>
    auto getShape(const Geometry& geometry)
    -> BRepWrapper< std::decay_t<decltype(Frackit::OCCUtilities::getShape(geometry))> >
    {
        using ShapeType = std::decay_t<decltype(Frackit::OCCUtilities::getShape(geometry))>;
        using WrapperType = BRepWrapper<ShapeType>;
        return WrapperType(Frackit::OCCUtilities::getShape(geometry));
    }

    ShapeWrapper readShape(const std::string& fileName)
    { return ShapeWrapper(Frackit::OCCUtilities::readShape(fileName)); }

    template<class WrappedShape>
    void write(const WrappedShape& wrappedShape, const std::string& fileName)
    {
        BRepTools::Write(wrappedShape.get(), fileName.c_str());
    }

} // end namespace OCCUtilities

template<class ctype>
void registerBRepUtilities(pybind11::module& module)
{
    namespace py = pybind11;

    // register point to shape conversions
    using P1 = Point<ctype, 1>; module.def("getShape", &OCCUtilities::getShape<P1>, "Returns the OCC BRep of a 1d point");
    using P2 = Point<ctype, 2>; module.def("getShape", &OCCUtilities::getShape<P2>, "Returns the OCC BRep of a 2d point");
    using P3 = Point<ctype, 3>; module.def("getShape", &OCCUtilities::getShape<P3>, "Returns the OCC BRep of a 3d point");

    // register segment to shape conversions
    using S1 = Segment<ctype, 1>; module.def("getShape", &OCCUtilities::getShape<S1>, "Returns the OCC BRep of a 1d segment");
    using S2 = Segment<ctype, 2>; module.def("getShape", &OCCUtilities::getShape<S2>, "Returns the OCC BRep of a 2d segment");
    using S3 = Segment<ctype, 3>; module.def("getShape", &OCCUtilities::getShape<S3>, "Returns the OCC BRep of a 3d segment");

    // register circle/ellipse to shape conversions
    module.def("getShape", &OCCUtilities::getShape<Circle<ctype, 3>>, "Returns the OCC BRep of a 3d circle");
    module.def("getShape", &OCCUtilities::getShape<Ellipse<ctype, 3>>, "Returns the OCC BRep of a 3d ellipse");
    module.def("getShape", &OCCUtilities::getShape<EllipseArc<ctype, 3>>, "Returns the OCC BRep of a 3d ellipse arc");

    // register quadrilateral to shape conversions
    module.def("getShape", &OCCUtilities::getShape<Quadrilateral<ctype, 3>>, "Returns the OCC BRep of a 3d quadrilateral");

    // register further geometry to shape conversions
    module.def("getShape", &OCCUtilities::getShape<Disk<ctype>>, "Returns the OCC BRep of a disk");
    module.def("getShape", &OCCUtilities::getShape<Box<ctype>>, "Returns the OCC BRep of a box");
    module.def("getShape", &OCCUtilities::getShape<Cylinder<ctype>>, "Returns the OCC BRep of a cylinder");
    module.def("getShape", &OCCUtilities::getShape<CylinderSurface<ctype>>, "Returns the OCC BRep of the lateral surface of a cylinder");

    // register read-in of shapes from a file
    module.def("readShape", &OCCUtilities::readShape, "Reads in the shapes from a file");

    // register write function for wrapped shapes
    using namespace OCCUtilities;
    module.def("write", &OCCUtilities::write<ShapeWrapper>, "writes a wrapped shape to a BRep file");
    module.def("write", &OCCUtilities::write<VertexWrapper>, "writes a wrapped vertex shape to a BRep file");
    module.def("write", &OCCUtilities::write<EdgeWrapper>, "writes a wrapped edge shape to a BRep file");
    module.def("write", &OCCUtilities::write<WireWrapper>, "writes a wrapped wire shape to a BRep file");
    module.def("write", &OCCUtilities::write<FaceWrapper>, "writes a wrapped face shape to a BRep file");
    module.def("write", &OCCUtilities::write<ShellWrapper>, "writes a wrapped shell shape to a BRep file");
    module.def("write", &OCCUtilities::write<SolidWrapper>, "writes a wrapped solid shape to a BRep file");
    module.def("write", &OCCUtilities::write<CompoundWrapper>, "writes a wrapped compound shape to a BRep file");
}

} // end namespace Frackit::Python

#endif
