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

#include <vector>
#include <string>
#include <type_traits>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

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

    template<class Geometry, std::enable_if_t<!IsBRepWrapper<Geometry>::value, int> = 0>
    auto getShape(const Geometry& geometry)
    -> BRepWrapper< std::decay_t<decltype(Frackit::OCCUtilities::getShape(geometry))> >
    {
        using ShapeType = std::decay_t<decltype(Frackit::OCCUtilities::getShape(geometry))>;
        using WrapperType = BRepWrapper<ShapeType>;
        return WrapperType(Frackit::OCCUtilities::getShape(geometry));
    }

    template<class Geometry, std::enable_if_t<IsBRepWrapper<Geometry>::value, int> = 0>
    Geometry getShape(const Geometry& geometry)
    { return geometry; }

    template<class ShapeWrapper>
    std::vector<VertexWrapper> getVertices(const ShapeWrapper& wrappedShape)
    {
        auto vertices = Frackit::OCCUtilities::getVertices(wrappedShape.get());

        std::vector< VertexWrapper > result;
        result.reserve(vertices.size());
        for (auto& v : vertices)
            result.emplace_back( std::move(v) );

        return result;
    }

    template<class ShapeWrapper>
    std::vector<EdgeWrapper> getEdges(const ShapeWrapper& wrappedShape)
    {
        auto edges = Frackit::OCCUtilities::getEdges(wrappedShape.get());

        std::vector< EdgeWrapper > result;
        result.reserve(edges.size());
        for (auto& e : edges)
            result.emplace_back( std::move(e) );

        return result;
    }

    template<class ShapeWrapper>
    std::vector<FaceWrapper> getFaces(const ShapeWrapper& wrappedShape)
    {
        auto faces = Frackit::OCCUtilities::getFaces(wrappedShape.get());

        std::vector< FaceWrapper > result;
        result.reserve(faces.size());
        for (auto& f : faces)
            result.emplace_back( std::move(f) );

        return result;
    }

    template<class ShapeWrapper>
    std::vector<WireWrapper> getWires(const ShapeWrapper& wrappedShape)
    {
        auto wires = Frackit::OCCUtilities::getWires(wrappedShape.get());

        std::vector< WireWrapper > result;
        result.reserve(wires.size());
        for (auto& w : wires)
            result.emplace_back( std::move(w) );

        return result;
    }

    template<class ShapeWrapper>
    std::vector<ShellWrapper> getShells(const ShapeWrapper& wrappedShape)
    {
        auto shells = Frackit::OCCUtilities::getShells(wrappedShape.get());

        std::vector< ShellWrapper > result;
        result.reserve(shells.size());
        for (auto& s : shells)
            result.emplace_back( std::move(s) );

        return result;
    }

    template<class ShapeWrapper>
    std::vector<SolidWrapper> getSolids(const ShapeWrapper& wrappedShape)
    {
        auto solids = Frackit::OCCUtilities::getSolids(wrappedShape.get());

        std::vector< SolidWrapper > result;
        result.reserve(solids.size());
        for (auto& s : solids)
            result.emplace_back( std::move(s) );

        return result;
    }

    template<class ShapeWrapper, class ctype>
    Box<ctype> getBoundingBox(const ShapeWrapper& wrappedShape)
    { return Frackit::OCCUtilities::getBoundingBox(wrappedShape.get()); }

    template<class ShapeWrapper, class ctype, int worldDim>
    ShapeWrapper translate(const ShapeWrapper& wrappedShape, const Vector<ctype, worldDim>& v)
    { return Frackit::OCCUtilities::translate(wrappedShape.get(), v); }

    template<class ShapeWrapper1, class ShapeWrapper2, class ctype>
    ShapeWrapper cut(const ShapeWrapper1& wrappedShape1,
                     const ShapeWrapper2& wrappedShape2,
                     ctype eps)
    { return Frackit::OCCUtilities::cut(wrappedShape1.get(), wrappedShape2.get(), eps); }

    template<class ShapeWrapper1, class ShapeWrapper2, class ctype>
    ShapeWrapper intersect(const ShapeWrapper1& wrappedShape1,
                           const ShapeWrapper2& wrappedShape2,
                           ctype eps)
    { return Frackit::OCCUtilities::intersect(wrappedShape1.get(), wrappedShape2.get(), eps); }

    template<class Wrapper, class ctype>
    ShapeWrapper fragment(const std::vector<Wrapper>& wrappedShapes, ctype eps)
    {
        std::vector<typename Wrapper::Shape> shapes;
        shapes.reserve(wrappedShapes.size());
        for (const auto& s : wrappedShapes)
            shapes.push_back(s.get());

        return Frackit::OCCUtilities::fragment(shapes, eps);
    }

    template<class Wrapper, class ctype>
    ShapeWrapper fuse(const std::vector<Wrapper>& wrappedShapes, ctype eps)
    {
        std::vector<typename Wrapper::Shape> shapes;
        shapes.reserve(wrappedShapes.size());
        for (const auto& s : wrappedShapes)
            shapes.push_back(s.get());

        return Frackit::OCCUtilities::fuse(shapes, eps);
    }

    ShapeWrapper readShape(const std::string& fileName)
    { return ShapeWrapper(Frackit::OCCUtilities::readShape(fileName)); }

    template<class WrappedShape>
    void write(const WrappedShape& wrappedShape, const std::string& fileName)
    {
        BRepTools::Write(wrappedShape.get(), fileName.c_str());
    }

} // end namespace OCCUtilities

namespace Detail {

    template<class WrapperType>
    void registerWrappedShapeConversions(pybind11::module& module)
    {
        module.def("getShape",
                   &OCCUtilities::getShape<WrapperType>,
                   "Returns the OCC BRep of a wrapped shape (for compatibility reasons, returns itself)");
    }

    template<class FirstType, class... WrapperTypes,
             std::enable_if_t<(sizeof...(WrapperTypes) > 0), int> = 0>
    void registerWrappedShapeConversions(pybind11::module& module)
    {
        registerWrappedShapeConversions<WrapperTypes...>(module);
        registerWrappedShapeConversions<FirstType>(module);
    }

    void registerSubShapeGetterFunctions(pybind11::module& module)
    {
        using namespace Frackit::Python::OCCUtilities;
        module.def("getVertices", &getVertices<ShapeWrapper>,    "returns the vertices of a wrapped shape");
        module.def("getVertices", &getVertices<VertexWrapper>,   "returns the vertices of a wrapped vertex");
        module.def("getVertices", &getVertices<EdgeWrapper>,     "returns the vertices of a wrapped edge");
        module.def("getVertices", &getVertices<WireWrapper>,     "returns the vertices of a wrapped wire");
        module.def("getVertices", &getVertices<FaceWrapper>,     "returns the vertices of a wrapped face");
        module.def("getVertices", &getVertices<ShellWrapper>,    "returns the vertices of a wrapped shell");
        module.def("getVertices", &getVertices<SolidWrapper>,    "returns the vertices of a wrapped solid");
        module.def("getVertices", &getVertices<CompoundWrapper>, "returns the vertices of a wrapped compound");

        module.def("getEdges", &getEdges<ShapeWrapper>,    "returns the edges of a wrapped shape");
        module.def("getEdges", &getEdges<EdgeWrapper>,     "returns the edges of a wrapped edge");
        module.def("getEdges", &getEdges<WireWrapper>,     "returns the edges of a wrapped wire");
        module.def("getEdges", &getEdges<FaceWrapper>,     "returns the edges of a wrapped face");
        module.def("getEdges", &getEdges<ShellWrapper>,    "returns the edges of a wrapped shell");
        module.def("getEdges", &getEdges<SolidWrapper>,    "returns the edges of a wrapped solid");
        module.def("getEdges", &getEdges<CompoundWrapper>, "returns the edges of a wrapped compound");

        module.def("getWires", &getWires<ShapeWrapper>,    "returns the wires of a wrapped shape");
        module.def("getWires", &getWires<WireWrapper>,     "returns the wires of a wrapped wire");
        module.def("getWires", &getWires<FaceWrapper>,     "returns the wires of a wrapped face");
        module.def("getWires", &getWires<ShellWrapper>,    "returns the wires of a wrapped shell");
        module.def("getWires", &getWires<SolidWrapper>,    "returns the wires of a wrapped solid");
        module.def("getWires", &getWires<CompoundWrapper>, "returns the wires of a wrapped compound");

        module.def("getFaces", &getFaces<ShapeWrapper>,    "returns the faces of a wrapped shape");
        module.def("getFaces", &getFaces<FaceWrapper>,     "returns the faces of a wrapped face");
        module.def("getFaces", &getFaces<ShellWrapper>,    "returns the faces of a wrapped shell");
        module.def("getFaces", &getFaces<SolidWrapper>,    "returns the faces of a wrapped solid");
        module.def("getFaces", &getFaces<CompoundWrapper>, "returns the faces of a wrapped compound");

        module.def("getShells", &getShells<ShapeWrapper>,    "returns the shells of a wrapped shape");
        module.def("getShells", &getShells<ShellWrapper>,    "returns the shells of a wrapped shell");
        module.def("getShells", &getShells<SolidWrapper>,    "returns the shells of a wrapped solid");
        module.def("getShells", &getShells<CompoundWrapper>, "returns the shells of a wrapped compound");

        module.def("getSolids", &getSolids<ShapeWrapper>,    "returns the solids of a wrapped shape");
        module.def("getSolids", &getSolids<SolidWrapper>,    "returns the solids of a wrapped solid");
        module.def("getSolids", &getSolids<CompoundWrapper>, "returns the solids of a wrapped compound");
    }

} // end namespace Detail

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

    // register "shape conversions" also for wrapped shapes
    Detail::registerWrappedShapeConversions<OCCUtilities::ShapeWrapper, OCCUtilities::VertexWrapper,
                                            OCCUtilities::EdgeWrapper,  OCCUtilities::WireWrapper,
                                            OCCUtilities::FaceWrapper,  OCCUtilities::ShellWrapper,
                                            OCCUtilities::SolidWrapper, OCCUtilities::CompoundWrapper>(module);

    // register getter functions for (sub-)shapes of a wrapped shape
    Detail::registerSubShapeGetterFunctions(module);

    // register read-in of shapes from a file
    module.def("readShape", &OCCUtilities::readShape, "Reads in the shapes from a file");

    // register bounding box computations for wrapped shapes
    using namespace OCCUtilities;
    module.def("getBoundingBox", &OCCUtilities::getBoundingBox<ShapeWrapper, ctype>,    "returns the bounding box of a wrapped shape");
    module.def("getBoundingBox", &OCCUtilities::getBoundingBox<VertexWrapper, ctype>,   "returns the bounding box of a wrapped vertex shape");
    module.def("getBoundingBox", &OCCUtilities::getBoundingBox<EdgeWrapper, ctype>,     "returns the bounding box of a wrapped edge shape");
    module.def("getBoundingBox", &OCCUtilities::getBoundingBox<WireWrapper, ctype>,     "returns the bounding box of a wrapped wire shape");
    module.def("getBoundingBox", &OCCUtilities::getBoundingBox<FaceWrapper, ctype>,     "returns the bounding box of a wrapped face shape");
    module.def("getBoundingBox", &OCCUtilities::getBoundingBox<ShellWrapper, ctype>,    "returns the bounding box of a wrapped shell shape");
    module.def("getBoundingBox", &OCCUtilities::getBoundingBox<SolidWrapper, ctype>,    "returns the bounding box of a wrapped solid shape");
    module.def("getBoundingBox", &OCCUtilities::getBoundingBox<CompoundWrapper, ctype>, "returns the bounding box of a wrapped compound shape");

    // register transformations
    module.def("translate", &OCCUtilities::translate<ShapeWrapper, ctype, 3>, "translation of a shape with a vector defined in 3d space");

    // register boolean operations for shape wrapper
    using namespace py::literals;
    module.def("cut",
               &OCCUtilities::cut<ShapeWrapper, ShapeWrapper, ctype>,
               "object"_a, "tool"_a, "tolerance"_a,
               "cuts the tool from the object");

    using namespace py::literals;
    module.def("intersect",
               &OCCUtilities::intersect<ShapeWrapper, ShapeWrapper, ctype>,
               "object1"_a, "object2"_a, "tolerance"_a,
               "returns the common part (intersection) between the two given shapes");

    using namespace py::literals;
    module.def("fragment",
               &OCCUtilities::fragment<ShapeWrapper, ctype>,
               "objects"_a, "tolerance"_a,
               "returns the fragments after intersection of the all given shapes");

    using namespace py::literals;
    module.def("fuse",
               &OCCUtilities::fuse<ShapeWrapper, ctype>,
               "objects"_a, "tolerance"_a,
               "fuse all given shapes into a single one");

    // register write function for wrapped shapes
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
