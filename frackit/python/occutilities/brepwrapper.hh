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
#ifndef FRACKIT_PYTHON_OCC_BREP_WRAPPER_HH
#define FRACKIT_PYTHON_OCC_BREP_WRAPPER_HH

#include <type_traits>

#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Shape.hxx>

#include <frackit/geometryutilities/name.hh>

namespace Frackit::Python::OCCUtilities {

// Wrapper class to be used around occ brep classes
template<class S>
class BRepWrapper
{
public:
    //! export underlying shape type
    using Shape = S;

    BRepWrapper(const Shape& shape) : shape_(shape) {}
    const Shape& get() const { return shape_; }
    std::string name() const { return geometryName(get()) + "_Wrapper"; }
private:
    Shape shape_;
};

// Wrapper around TopoDS_Compound
template<>
class BRepWrapper<TopoDS_Compound>
{
    using Compound = TopoDS_Compound;

public:
    //! export underlying shape type
    using Shape = Compound;

    //! Default constructor
    BRepWrapper()
    : builder_()
    {
        builder_.MakeCompound(compound_);
    }

    //! Constructor from a compound
    BRepWrapper(const Compound& c)
    : compound_(c)
    , builder_()
    {}

    //! return the compound
    const Compound& get() const
    { return compound_; }

    //! return the name of the wrapper class
    std::string name() const
    { return "TopoDS_Compound_Wrapper"; }

    //! add a wrapped shape to the compound
    template<class ShapeWrapper>
    void add(const ShapeWrapper& shape)
    {
        builder_.Add(compound_, shape.get());
    }

private:
    Compound compound_;
    BRep_Builder builder_;
};

// aliases to be used for the different wrapper classes
using ShapeWrapper = BRepWrapper<TopoDS_Shape>;
using VertexWrapper = BRepWrapper<TopoDS_Vertex>;
using EdgeWrapper = BRepWrapper<TopoDS_Edge>;
using WireWrapper = BRepWrapper<TopoDS_Wire>;
using FaceWrapper = BRepWrapper<TopoDS_Face>;
using ShellWrapper = BRepWrapper<TopoDS_Shell>;
using SolidWrapper = BRepWrapper<TopoDS_Solid>;
using CompoundWrapper = BRepWrapper<TopoDS_Compound>;

// helper struct to check if something is a BRep wrapper
template<class T>
struct IsBRepWrapper : public std::false_type {};
template<class S>
struct IsBRepWrapper<BRepWrapper<S>> : public std::true_type {};

} // end namespace Frackit::Python::OCCUtilities

#endif
