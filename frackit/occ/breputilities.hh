// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
/*!
 * \file
 * \brief Contains utility functionality for
 *        operations on objects of the BRep package
 *        as well as parsing internal geometries into
 *        BRep representations.
 */
#ifndef FRACKIT_BREP_UTILITIES_HH
#define FRACKIT_BREP_UTILITIES_HH

#include <vector>
#include <algorithm>
#include <stdexcept>

// Handle class used by OpenCascade
#include <Standard_Handle.hxx>

// objects from geometric processors package
#include <gp_Pnt.hxx>
#include <gp_Dir.hxx>
#include <gp_Elips.hxx>

// shape classes from TopoDS package
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Shape.hxx>
#include <TopTools_ListOfShape.hxx>

// explorer for shape objects
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopAbs_ShapeEnum.hxx>

// Algorithm for boolean operations on shapes
#include <BRepAlgoAPI_BuilderAlgo.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Cut.hxx>

// converter between BRep and TopoDS
#include <Bnd_Box.hxx>
#include <BRep_Tool.hxx>
#include <BRepTools.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>

// BRep primitives and operations
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeBox.hxx>

// internal geometry classes
#include <frackit/geometry/circle.hh>
#include <frackit/geometry/ellipse.hh>
#include <frackit/geometry/ellipsearc.hh>
#include <frackit/geometry/disk.hh>
#include <frackit/geometry/quadrilateral.hh>
#include <frackit/geometry/cylinder.hh>
#include <frackit/geometry/cylindersurface.hh>
#include <frackit/geometry/box.hh>

#include "gputilities.hh"
#include "geomutilities.hh"

namespace Frackit {
namespace OCCUtilities {

    //! converts a vertex shape into a point
    Point<double, 3> point(const TopoDS_Vertex& v)
    { return point(BRep_Tool::Pnt(v)); }

    //! Create an edge shape from two points
    template<class ctype, int worldDim>
    TopoDS_Edge makeEdge(const Point<ctype, worldDim>& source,
                         const Point<ctype, worldDim>& target)
    {
        // create TopoDS_Edge of the segment
        BRepBuilderAPI_MakeEdge segment(point(source), point(target));
        segment.Build();
        if (!segment.IsDone())
            throw std::runtime_error("Could not create segment edge");
        return TopoDS::Edge(segment.Shape());
    }

    //! get the BRep of a geometry (overloads provided)
    template<class Geo>
    TopoDS_Shape getShape(const Geo& p)
    { throw std::runtime_error(std::string("getShape() not implemented for " + Geo::name())); }

    //! get the BRep of a point
    template<class ctype, int worldDim>
    TopoDS_Vertex getShape(const Point<ctype, worldDim>& p)
    {
        BRepBuilderAPI_MakeVertex vertex(point(p));
        vertex.Build();
        if (!vertex.IsDone())
            throw std::runtime_error("Could not create vertex shape");
        return TopoDS::Vertex(vertex.Shape());
    }

    //! get the BRep of a segment
    template<class ctype, int worldDim>
    TopoDS_Edge getShape(const Segment<ctype, worldDim>& segment)
    { return makeEdge(segment.source(), segment.target()); }

    //! get the BRep of a cylinder
    template<class ctype>
    TopoDS_Solid getShape(const Cylinder<ctype>& cylinder)
    {
        const auto& lateral = cylinder.lateralFace();
        const auto& bottom = lateral.lowerBoundingCircle();
        auto axis = direction(bottom.normal());
        auto base1 = direction(bottom.base1());
        auto center = point(bottom.center());
        BRepPrimAPI_MakeCylinder makeCylinder(gp_Ax2(center, axis, base1),
                                              bottom.radius(),
                                              lateral.height(),
                                              2.0*M_PI);
        makeCylinder.Build();
        if (!makeCylinder.IsDone())
            throw std::runtime_error(std::string("Could not build cylinder"));
        return TopoDS::Solid(makeCylinder.Shape());
    }

    //! get the BRep of a lateral cylinder surface
    template<class ctype>
    TopoDS_Face getShape(const CylinderSurface<ctype>& cylSurface)
    {
        const auto& bottom = cylSurface.lowerBoundingCircle();
        auto axis = direction(bottom.normal());
        auto base1 = direction(bottom.base1());
        auto center = point(bottom.center());
        BRepPrimAPI_MakeCylinder makeCylinder(gp_Ax2(center, axis, base1),
                                              bottom.radius(),
                                              cylSurface.height(),
                                              2.0*M_PI);
        makeCylinder.Build();
        if (!makeCylinder.IsDone())
            throw std::runtime_error(std::string("Could not build cylinder"));
        return makeCylinder.Cylinder().LateralFace();
    }

    //! get the BRep of a 3-dimensional ellipse
    template<class ctype>
    TopoDS_Wire getShape(const Ellipse<ctype, 3>& ellipse)
    {
        gp_Dir normal = direction(ellipse.normal());
        gp_Dir majorAx = direction(ellipse.majorAxis());
        gp_Pnt center = point(ellipse.center());
        gp_Elips gpEllipse(gp_Ax2(center, normal, majorAx),
                           ellipse.majorAxisLength(),
                           ellipse.minorAxisLength());

        TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(gpEllipse);
        return BRepBuilderAPI_MakeWire(edge);
    }

    //! get the BRep of a 3-dimensional circle
    template<class ctype>
    TopoDS_Wire getShape(const Circle<ctype, 3>& circ)
    {
        return getShape(Ellipse<ctype, 3>(circ.center(),
                                          circ.base1(),
                                          circ.base2(),
                                          circ.radius(),
                                          circ.radius()));
    }

    //! get the BRep of a 3-dimensional ellipse arc
    template<class ctype>
    TopoDS_Edge getShape(const EllipseArc<ctype, 3>& arc)
    { return BRepBuilderAPI_MakeEdge( getGeomHandle(arc) ); }

    //! get the BRep of a quadrilateral in 3d space
    template<class ctype>
    TopoDS_Face getShape(const Quadrilateral<ctype, 3>& quad)
    {
        const auto v1 = getShape(quad.corner(0));
        const auto v2 = getShape(quad.corner(1));
        const auto v3 = getShape(quad.corner(2));
        const auto v4 = getShape(quad.corner(3));

        TopoDS_Edge e1 = BRepBuilderAPI_MakeEdge(v1, v2);
        TopoDS_Edge e2 = BRepBuilderAPI_MakeEdge(v2, v4);
        TopoDS_Edge e3 = BRepBuilderAPI_MakeEdge(v4, v3);
        TopoDS_Edge e4 = BRepBuilderAPI_MakeEdge(v3, v1);

        TopoDS_Wire wire = BRepBuilderAPI_MakeWire(e1, e2, e3, e4);
        return BRepBuilderAPI_MakeFace(wire);
    }

    //! get the BRep of a disk
    template<class ctype>
    TopoDS_Face getShape(const Disk<ctype>& disk)
    {
        gp_Dir normal = direction(disk.normal());
        gp_Dir majorAx = direction(disk.majorAxis());
        gp_Pnt center = point(disk.center());
        gp_Elips ellipse(gp_Ax2(center, normal, majorAx),
                         disk.majorAxisLength(),
                         disk.minorAxisLength());

        TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(ellipse);
        TopoDS_Wire wire = BRepBuilderAPI_MakeWire(edge);
        return BRepBuilderAPI_MakeFace(wire);
    }

    //! get the BRep of a box
    template<class ctype>
    TopoDS_Solid getShape(const Box<ctype>& box)
    {
        const gp_Pnt pMin(box.xMin(), box.yMin(), box.zMin());
        const gp_Pnt pMax(box.xMax(), box.yMax(), box.zMax());
        BRepPrimAPI_MakeBox makeBox(pMin, pMax);
        makeBox.Build();
        if (!makeBox.IsDone())
            throw std::runtime_error(std::string("Could not build box"));
        return TopoDS::Solid(makeBox.Shape());
    }

    //! get the BRep of shapes (for compatibilty reasons)
    const TopoDS_Vertex& getShape(const TopoDS_Vertex& v) { return v; }
    const TopoDS_Edge& getShape(const TopoDS_Edge& e) { return e; }
    const TopoDS_Face& getShape(const TopoDS_Face& f) { return f; }
    const TopoDS_Wire& getShape(const TopoDS_Wire& w) { return w; }
    const TopoDS_Solid& getShape(const TopoDS_Solid& s) { return s; }

    //! Get the vertices of a shape
    template<class Shape>
    std::vector<TopoDS_Vertex> getVertices(const Shape& shape)
    {
        TopExp_Explorer explorer;
        std::vector<TopoDS_Vertex> vertices;
        for (explorer.Init(shape, TopAbs_VERTEX); explorer.More(); explorer.Next())
            vertices.push_back(TopoDS::Vertex(explorer.Current()));
        return vertices;
    }

    //! Get the edges of a shape
    template<class Shape>
    std::vector<TopoDS_Edge> getEdges(const Shape& shape)
    {
        TopExp_Explorer explorer;
        std::vector<TopoDS_Edge> edges;
        for (explorer.Init(shape, TopAbs_EDGE); explorer.More(); explorer.Next())
            edges.push_back(TopoDS::Edge(explorer.Current()));
        return edges;
    }

    //! Get the faces of a shape
    template<class Shape>
    std::vector<TopoDS_Face> getFaces(const Shape& shape)
    {
        TopExp_Explorer explorer;
        std::vector<TopoDS_Face> faces;
        for (explorer.Init(shape, TopAbs_FACE); explorer.More(); explorer.Next())
            faces.push_back(TopoDS::Face(explorer.Current()));
        return faces;
    }

    //! Get the wires of a shape
    template<class Shape>
    std::vector<TopoDS_Wire> getWires(const Shape& shape)
    {
        std::vector<TopoDS_Wire> wires;
        for (TopExp_Explorer explorer(shape, TopAbs_WIRE); explorer.More(); explorer.Next())
            wires.push_back(TopoDS::Wire(explorer.Current()));
        return wires;
    }

    //! Get the shells of a shape
    template<class Shape>
    std::vector<TopoDS_Shell> getShells(const Shape& shape)
    {
        std::vector<TopoDS_Shell> shells;
        for (TopExp_Explorer explorer(shape, TopAbs_SHELL); explorer.More(); explorer.Next())
            shells.push_back(TopoDS::Shell(explorer.Current()));
        return shells;
    }

    //! Get the solids of a shape
    template<class Shape>
    std::vector<TopoDS_Solid> getSolids(const Shape& shape)
    {
        std::vector<TopoDS_Solid> solids;
        for (TopExp_Explorer explorer(shape, TopAbs_SOLID); explorer.More(); explorer.Next())
            solids.push_back(TopoDS::Solid(explorer.Current()));
        return solids;
    }

    //! Get bounding box
    template<class Shape, class ctype = double>
    Box<ctype> getBoundingBox(const Shape& shape)
    {
        Bnd_Box bndBox;
        BRepBndLib::Add(shape, bndBox);

        ctype xMin, yMin, zMin, xMax, yMax, zMax;
        bndBox.Get(xMin, yMin, zMin, xMax, yMax, zMax);
        return Box<ctype>(xMin, yMin, zMin, xMax, yMax, zMax);
    }

    //! Convenience function to get the shape resulting from the cut of an object with a tool
    template<class ctype>
    TopoDS_Shape cut(const TopoDS_Shape& object, const TopoDS_Shape& tool, ctype eps)
    {
        BRepAlgoAPI_Cut cut(object, tool);
        cut.SetRunParallel(false);
        cut.SetFuzzyValue(eps);
        cut.Build();
        if (!cut.IsDone())
            throw std::runtime_error(std::string("Could not perform cut operation"));

        return cut.Shape();
    }

    //! Convenience function to get the shape resulting from the intersection of two objects
    template<class ctype>
    TopoDS_Shape intersect(const TopoDS_Shape& object1, const TopoDS_Shape& object2, ctype eps)
    {
        BRepAlgoAPI_Common common(object1, object2);
        common.SetRunParallel(false);
        common.SetFuzzyValue(eps);
        common.Build();
        if (!common.IsDone())
            throw std::runtime_error(std::string("Common operation failed"));

        return common.Shape();
    }

    //! Convenience function to get the shape resulting from the fragmentation of a set of objects
    template<class ctype>
    TopoDS_Shape fragment(const std::vector<TopoDS_Shape>& objects, ctype eps)
    {
        TopTools_ListOfShape shapes;
        for (const auto& shape : objects)
            shapes.Append(shape);

        BRepAlgoAPI_BuilderAlgo fragments;
        fragments.SetRunParallel(false);
        fragments.SetArguments(shapes);
        fragments.SetFuzzyValue(eps);
        fragments.Build();
        if (!fragments.IsDone())
            throw std::runtime_error(std::string("Could not perform segment fragmentation"));

        return fragments.Shape();
    }

    /*!
     * \brief Determines the overlap edges between two sets of edges.
     * \param edges1 the first set of edges
     * \param edges2 the second set of edges
     * \note This assumes that the sets of edges don't contain any duplicates!
     * \returns A vector with the edges on which the two faces coincide.
     */
    std::vector<TopoDS_Edge> getOverlapEdges(const std::vector<TopoDS_Edge>& edges1,
                                             const std::vector<TopoDS_Edge>& edges2)
    {
        std::vector<TopoDS_Edge> overlapEdges;
        for (const auto& e : edges1)
            if ( std::count_if(edges2.begin(),
                               edges2.end(),
                               [&e] (const auto& e2) { return e2.IsSame(e); }) )
                overlapEdges.push_back(e);
        return overlapEdges;
    }

    /*!
     * \brief Determines those edges along which two faces coincide.
     * \param face1 the shape of the first face
     * \param face2 the shape of the second face
     * \returns A vector with the edges on which the two faces coincide.
     */
    std::vector<TopoDS_Edge> getOverlapEdges(const TopoDS_Face& face1,
                                             const TopoDS_Face& face2)
    { return getOverlapEdges(getEdges(face1), getEdges(face2)); }

} // end namespace OCCUtilities
} // end namespace Frackit

#endif // FRACKIT_BREP_UTILITIES_HH
