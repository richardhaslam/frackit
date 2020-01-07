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
 * \brief \todo TODO doc me.
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
#include <BRep_Tool.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>

// BRep primitives and operations
#include <BRepPrimAPI_MakeCylinder.hxx>

// internal geometry classes
#include <frackit/geometry/circle.hh>
#include <frackit/geometry/ellipse.hh>
#include <frackit/geometry/ellipsearc.hh>
#include <frackit/geometry/disk.hh>
#include <frackit/geometry/cylinder.hh>
#include <frackit/geometry/cylindersurface.hh>

#include "gputilities.hh"
#include "geomutilities.hh"

namespace Frackit {
namespace OCCUtilities {

    //! converts a vertex shape into a point
    Point<double, 3> point(const TopoDS_Vertex& v)
    { return point(BRep_Tool::Pnt(v)); }

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
    { return BRepBuilderAPI_MakeWire( getGeomHandle(arc) ); }

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

    //! Create an edge shape from two points
    template<class ctype, int dim>
    TopoDS_Edge makeEdge(const Point<ctype, dim>& source,
                         const Point<ctype, dim>& target)
    {
        // create TopoDS_Edge of the segment
        BRepBuilderAPI_MakeEdge segment(point(source), point(target));
        segment.Build();
        if (!segment.IsDone())
            throw std::runtime_error("Could not create segment edge");
        return TopoDS::Edge(segment.Shape());
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

    /*!
     * \brief Determines from a set of edges which edges are connected by a vertex.
     * \param edges vector of edges
     * \returns A vector of size edges.size() where in each entry a vector of
     *          connectivity relationships is stored. The connectivity relationships
     *          are stored in an std::array<unsignet int, 3>, in which the first entry
     *          is the index of a connected edge, the second entry is the index of the
     *          coinciding vertex on this edge (0 - first vertex, 1 - second vertex) and
     *          the third entry is the index of the coinciding vertex on the connected edge.
     */
    std::vector< std::vector< std::array<unsigned int, 3> > >
    getConnectivityInfo(const std::vector<TopoDS_Edge>& edges)
    {
        using Info = std::array<unsigned int, 3>;
        std::vector< std::vector<Info> > connectivity(edges.size());
        for (unsigned int i = 0; i < edges.size(); ++i)
        {
            for (unsigned int j = i+1; j < edges.size(); ++j)
            {
                auto& ci = connectivity[i];
                auto& cj = connectivity[j];
                const auto& vi1 = TopExp::FirstVertex(edges[i]);
                const auto& vj1 = TopExp::FirstVertex(edges[j]);
                const auto& vi2 = TopExp::LastVertex(edges[i]);
                const auto& vj2 = TopExp::LastVertex(edges[j]);
                if (vi1.IsSame(vj1)) { ci.emplace_back(Info({j, 0, 0})); cj.emplace_back(Info({i, 0, 0})); }
                else if (vi1.IsSame(vj2)) { ci.emplace_back(Info({j, 0, 1})); cj.emplace_back(Info({i, 1, 0})); }
                else if (vi2.IsSame(vj1)) { ci.emplace_back(Info({j, 1, 0})); cj.emplace_back(Info({i, 0, 1})); }
                else if (vi2.IsSame(vj2)) { ci.emplace_back(Info({j, 1, 1})); cj.emplace_back(Info({i, 1, 1})); }
            }
        }
        return connectivity;
    }

    /*!
     * \brief Determines the corner points of a set of edges, where internal vertices
     *        which connect two edges are omitted and only the corners of segments
     *        possibly consisting of several edges are extracted. This information
     *        is taken from the provided connectivity map.
     * \param edges vector of edges
     * \param connectivity connectivity information (see getConnectivityInfo())
     * \note The set of edges must not contain duplicates!
     *       Use the function getConnectivityInfo() to obtain the connectivity
     *       info of the set of edges for which you want to determine the corners.
     * \note This assumes that the set of edges describes a single connected curve.
     *       If the curve is a closed periodic curve (i.e. each edge has two neighbors)
     *       the returned start and end points are the same and are taken from the first edge.
     * \returns A pair of two-dimensional arrays, where each array contains the edge and corner
     *          index (0 - first vertex, 1 - second vertex) that compose a corner.
     */
    std::pair< std::array<unsigned int, 2>, std::array<unsigned int, 2> >
    getBoundaryVertices(const std::vector<TopoDS_Edge>& edges,
                        const std::vector<std::vector<std::array<unsigned int, 3>>>& connectivity)
    {
        using CornerInfo = std::array<unsigned int, 2>;
        if (edges.size() == 1)
            return std::make_pair( CornerInfo({0, 0}), CornerInfo({0, 1}) );

        assert(std::all_of(connectivity.begin(),
                           connectivity.end(),
                           [] (const auto& c) { return c.size() != 0; }));
        assert(std::all_of(connectivity.begin(),
                           connectivity.end(),
                           [] (const auto& c) { return c.size() <= 2; }));

        // for a closed curve, arbitrarily return the first corner as start&end
        if (std::all_of(connectivity.begin(),
                        connectivity.end(),
                        [] (const auto& c) { return c.size() == 2; }))
            return std::make_pair( CornerInfo({0, 0}), CornerInfo({0, 0}) );

        assert(std::count_if(connectivity.begin(),
                             connectivity.end(),
                             [] (const auto& c) { return c.size() == 1; }) == 2);

        std::vector<CornerInfo> corners;
        for (unsigned int i = 0; i < edges.size(); ++i)
            if (connectivity[i].size() == 1)
                connectivity[i][0][1] == 1 ? corners.emplace_back(CornerInfo({i, 0}))
                                           : corners.emplace_back(CornerInfo({i, 1}));

        assert(corners.size() == 2);
        return std::make_pair(std::move(corners[0]), std::move(corners[1]));
    }

    /*!
     * \brief Determines the corner points of a set of edges, where internal vertices
     *        which connect two edges are omitted and only the corners of segments
     *        possibly consisting of several edges are extracted. This information
     *        is taken from the provided connectivity map.
     * \param edges vector of edges
     * \note The set of edges must not contain duplicates!
     * \note This assumes that the set of edges describes a single connected curve.
     *       If the curve is a closed periodic curve (i.e. each edge has two neighbors)
     *       the returned start and end points are the same and are taken from the first edge.
     * \returns A pair of two-dimensional arrays, where each array contains the edge and corner
     *          index (0 - first vertex, 1 - second vertex) that compose a corner.
     */
    std::pair< std::array<unsigned int, 2>, std::array<unsigned int, 2> >
    getBoundaryVertices(const std::vector<TopoDS_Edge>& edges)
    { return getBoundaryVertices(edges, getConnectivityInfo(edges)); }

} // end namespace OCCUtilities
} // end namespace Frackit

#endif // FRACKIT_BREP_UTILITIES_HH
