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
#ifndef FRACKIT_PYTHON_GEOMETRY_OCC_BREP_WRAPPERS_HH
#define FRACKIT_PYTHON_GEOMETRY_OCC_BREP_WRAPPERS_HH

#include <type_traits>
#include <pybind11/pybind11.h>

#include <TopoDS_Shape.hxx>

#include "registerdimensionproperties.hh"
#include "brepwrapper.hh"

namespace Frackit::Python {

namespace Detail {

    template<class Shape>
    void registerBRepWrapper(pybind11::module& module,
                             const std::string& className)
    {
        using Wrapper = BRepWrapper<Shape>;
        pybind11::class_<Wrapper> cls(module, className.c_str());
        cls.def("name", &Wrapper::name, "return the name of the wrapper");
        registerDimensionProperties(cls);

        // add constructor from other shape wrappers
        if constexpr(std::is_same_v<Shape, TopoDS_Shape>)
        {
            namespace py = pybind11;
            cls.def(py::init<ShapeWrapper>(), "construction from another brep wrapper");
            cls.def(py::init<VertexWrapper>(), "construction from another brep wrapper");
            cls.def(py::init<EdgeWrapper>(), "construction from another brep wrapper");
            cls.def(py::init<WireWrapper>(), "construction from another brep wrapper");
            cls.def(py::init<FaceWrapper>(), "construction from another brep wrapper");
            cls.def(py::init<ShellWrapper>(), "construction from another brep wrapper");
            cls.def(py::init<SolidWrapper>(), "construction from another brep wrapper");
            cls.def(py::init<CompoundWrapper>(), "construction from another brep wrapper");
        }
    }

    void registerCompoundWrapper(pybind11::module& module)
    {
        using Wrapper = CompoundWrapper;
        pybind11::class_<Wrapper> cls(module, "OCCCompoundWrapper");
        cls.def(pybind11::init<>());
        registerDimensionProperties(cls);

        // functions to add shapes to the compound
        cls.def("name", &Wrapper::name, "return the name of the wrapper");
        cls.def("add", &Wrapper::add<ShapeWrapper>, "Add a shape to the compound");
        cls.def("add", &Wrapper::add<VertexWrapper>, "Add a vertex to the compound");
        cls.def("add", &Wrapper::add<EdgeWrapper>, "Add an edge to the compound");
        cls.def("add", &Wrapper::add<WireWrapper>, "Add a wire to the compound");
        cls.def("add", &Wrapper::add<FaceWrapper>, "Add a face to the compound");
        cls.def("add", &Wrapper::add<ShellWrapper>, "Add a shell to the compound");
        cls.def("add", &Wrapper::add<SolidWrapper>, "Add a solid to the compound");
    }

} // end namespace Detail

void registerBRepWrappers(pybind11::module& module)
{
    Detail::registerBRepWrapper<TopoDS_Shape>(module, "OCCShapeWrapper");
    Detail::registerBRepWrapper<TopoDS_Vertex>(module, "OCCVertexWrapper");
    Detail::registerBRepWrapper<TopoDS_Edge>(module, "OCCEdgeWrapper");
    Detail::registerBRepWrapper<TopoDS_Wire>(module, "OCCWireWrapper");
    Detail::registerBRepWrapper<TopoDS_Face>(module, "OCCFaceWrapper");
    Detail::registerBRepWrapper<TopoDS_Shell>(module, "OCCShellWrapper");
    Detail::registerBRepWrapper<TopoDS_Solid>(module, "OCCSolidWrapper");
    Detail::registerCompoundWrapper(module);
}

} // end namespace Frackit::Python

#endif
