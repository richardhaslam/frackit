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
#ifndef FRACKIT_INTERSECT_HH
#define FRACKIT_INTERSECT_HH

#include <algorithm>
#include <type_traits>
#include <stdexcept>
#include <cassert>
#include <utility>
#include <variant>
#include <vector>
#include <cmath>

// some algorithms/classes receive handles instead of objects
#include <Standard_Handle.hxx>

// classes from the geometric processors package
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <gp_Dir.hxx>
#include <gp_Ax2.hxx>
#include <gp_Elips.hxx>

// classes from the geometry package
#include <GeomAPI_IntSS.hxx>
#include <GeomAPI_IntCS.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <Geom_Plane.hxx>
#include <Geom_Surface.hxx>
#include <Geom_Line.hxx>
#include <Geom_Curve.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <Geom_Ellipse.hxx>

// builders for TopoDS_Shapes
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRep_Tool.hxx>

// BRep primitives and operations
#include <BRepPrimAPI_MakeCylinder.hxx>

// shapes to be passed to intersection algorithms
#include <TopTools_ListOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Vertex.hxx>

// class to explore the resulting shape of an intersection
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopAbs_Orientation.hxx>

// Algorithm for boolean operations on shapes
#include <BRepAlgoAPI_BuilderAlgo.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Cut.hxx>

// Geometries of intersections
#include <frackit/geometry/point.hh>
#include <frackit/geometry/line.hh>
#include <frackit/geometry/segment.hh>
#include <frackit/geometry/direction.hh>
#include <frackit/geometry/plane.hh>
#include <frackit/geometry/disk.hh>
#include <frackit/geometry/cylindricalsurface.hh>

// base tolerances for floating point comparisons
#include <frackit/geometry/precision.hh>

// utility functionality
#include <frackit/common/utilities.hh>
#include <frackit/common/math.hh>

#include "intersectiontraits.hh"

namespace Frackit {

//! \todo TODO Doc me.
template<class Geom1, class Geom2, class ctype = double>
EmptyIntersection<0> intersect(const Geom1& geo1, const Geom2& geo2, ctype eps = 0.0)
{
    std::string msg = "Intersection algorithm between \"";
    msg += geo1.name();
    msg += "\" and \"";
    msg += geo2.name();
    msg += "\" not implemented";
    throw std::runtime_error(msg);
}

//! \todo TODO Doc me.
template<class ctype, int wd>
Intersection< Segment<ctype, wd>, Segment<ctype, wd> >
intersect(const Segment<ctype, wd>& segment1,
          const Segment<ctype, wd>& segment2,
          ctype eps)
{
    using ResultType = Intersection<Segment<ctype, wd>, Segment<ctype, wd>>;
    using Segment = Frackit::Segment<ctype, wd>;
    using Direction = typename Segment::Direction;

    const auto& s1 = segment1.source();
    const auto& s2 = segment2.source();
    const auto& t1 = segment1.target();
    const auto& t2 = segment2.target();

    // if the segments are parallel, there might be an overlap
    const auto dir1 = Direction(Vector<ctype, wd>(s1, t1));
    const auto dir2 = Direction(Vector<ctype, wd>(s2, t2));
    gp_Vec v1(dir1.x(), dir1.y(), dir1.z());
    gp_Vec v2(dir2.x(), dir2.y(), dir2.z());

    if (v1.IsParallel(v2, Precision<ctype>::angular()))
    {
        // check if the lines overlap
        if (!segment1.supportingLine().contains(s2, eps))
            return ResultType( EmptyIntersection<wd>() );

        // find overlap
        const bool s2On1 = segment1.contains(s2, eps, false);
        const bool t2On1 = segment1.contains(t2, eps, false);
        if (s2On1 && t2On1) return ResultType( segment2 );

        const bool s1On2 = segment2.contains(s1, eps, false);
        const bool t1On2 = segment2.contains(t1, eps, false);
        if (s1On2 && t1On2) return ResultType( segment1 );

        if (s2On1 && s1On2) return ResultType( Segment(s2, s1) );
        if (s2On1 && t1On2) return ResultType( Segment(s2, t1) );

        if (t2On1 && s1On2) return ResultType( Segment(t2, s1) );
        if (t2On1 && t1On2) return ResultType( Segment(t2, t1) );
    }

    // Fragment the lines to find the intersection point
    std::vector<TopoDS_Shape> segmentEdges;
    segmentEdges.emplace_back(OCCUtilities::makeEdge(s1, t1));
    segmentEdges.emplace_back(OCCUtilities::makeEdge(s2, t2));

    const auto fragmentShape = OCCUtilities::fragment(segmentEdges, 0.1*eps);
    const auto edges = OCCUtilities::getEdges(fragmentShape);

    // If the result has 3 edges, a corner of one segment splits the other segment.
    // If it has 4 edges, both segments are split by the intersection point.
    // In either case, find the newly created point which is the intersection.
    const auto numEdges = edges.size();
    if (numEdges == 3 || numEdges == 4)
    {
        for (const auto& edge : edges)
        {
            for (TopExp_Explorer explorer(edge, TopAbs_VERTEX); explorer.More(); explorer.Next())
            {
                auto curPoint = OCCUtilities::point(TopoDS::Vertex(explorer.Current()));
                if (!curPoint.isEqual(s1, eps) && !curPoint.isEqual(t1, eps)
                    && !curPoint.isEqual(s2, eps) && !curPoint.isEqual(t2, eps))
                    return ResultType( std::move(curPoint) );
            }
        }

        throw std::runtime_error(std::string("Could not find intersection point"));
    }

    // a pair of corners being equal is the last possible option
    if (s1.isEqual(s2, eps)) return ResultType( s1 );
    else if (s1.isEqual(t2, eps)) return ResultType( s1 );
    else if (t1.isEqual(s2, eps)) return ResultType( t1 );
    else if (t1.isEqual(t2, eps)) return ResultType( t1 );

    return ResultType( EmptyIntersection<wd>() );
}

//! \todo TODO Doc me.
template<class ctype, int wd>
Intersection< Segment<ctype, wd>, Segment<ctype, wd> >
intersect(const Segment<ctype, wd>& segment1,
          const Segment<ctype, wd>& segment2)
{
    using std::min;
    return intersect(segment1,
                     segment2,
                     0.5*Precision<ctype>::confusion()
                     *min(segment1.length(), segment2.length()));
}

template<class ctype>
Intersection< Plane<ctype, 3>, Plane<ctype, 3> >
intersect(const Plane<ctype, 3>& plane1,
          const Plane<ctype, 3>& plane2,
          ctype eps)
{
    using ResultType = Intersection< Plane<ctype, 3>, Plane<ctype, 3> >;

    const auto d1 = OCCUtilities::direction(plane1.normal());
    const auto d2 = OCCUtilities::direction(plane2.normal());

    // if planes are identical, the result is a plane or empty
    if (gp_Vec(d1).IsParallel(gp_Vec(d2), Precision<ctype>::angular()))
    {
        if (plane1.contains(plane2.supportingPoint(), eps))
            return ResultType( plane1 );
        else
            return ResultType( EmptyIntersection<3>() );
    }

    // find the intersection line between the two planes
    const auto c1 = OCCUtilities::point(plane1.supportingPoint());
    const auto c2 = OCCUtilities::point(plane2.supportingPoint());

    Handle(Geom_Plane) gpPlane1 = new Geom_Plane(c1, d1);
    Handle(Geom_Plane) gpPlane2 = new Geom_Plane(c2, d2);
    GeomAPI_IntSS interSection(gpPlane1, gpPlane2, eps);
    if (!interSection.IsDone())
        throw std::runtime_error(std::string("Could not perform plane-plane intersection"));

    assert(interSection.NbLines() == 1);
    assert(!interSection.Line(1)->IsClosed());
    assert(interSection.Line(1)->IsCN(1));

    // The result is an infinite line.
    // The parameter range is -inf -> +inf
    // We evaluate a support point and direction for param = 0.
    const auto& sp = interSection.Line(1)->Value(0.0);
    const auto& dir = interSection.Line(1)->DN(0.0, /*derivative order*/1);

    // create Line object from support point and direction
    using Line = Frackit::Line<ctype, 3>;
    return ResultType( Line(OCCUtilities::point(sp),
                            OCCUtilities::direction(dir)) );
}

template<class ctype>
Intersection< Plane<ctype, 3>, Plane<ctype, 3> >
intersect(const Plane<ctype, 3>& plane1, const Plane<ctype, 3>& plane2)
{ return intersect(plane1, plane2, Precision<ctype>::confusion()); }

template<class ctype>
Intersection< Plane<ctype, 3>, Line<ctype, 3> >
intersect(const Plane<ctype, 3>& plane,
          const Line<ctype, 3>& line,
          ctype eps)
{
    using ResultType = Intersection< Plane<ctype, 3>, Line<ctype, 3> >;

    // if the line is not parallel to the plane, the result can only be a point
    const auto lineDir = OCCUtilities::direction(line.direction());
    const auto planeNormal = OCCUtilities::direction(plane.normal());
    if (!lineDir.IsNormal(planeNormal, Precision<ctype>::angular()))
    {
        const auto linePoint = OCCUtilities::point(line.supportingPoint());
        const auto planePoint = OCCUtilities::point(plane.supportingPoint());

        // let the the geometry package compute the intersection
        Handle(Geom_Surface) planeHandle = new Geom_Plane(planePoint, planeNormal);
        Handle(Geom_Line) lineHandle = new Geom_Line(linePoint, lineDir);
        GeomAPI_IntCS interSection(lineHandle, planeHandle);
        if (!interSection.IsDone())
            throw std::runtime_error(std::string("Could not perform disk-line intersection"));

        assert(interSection.NbSegments() == 0);
        assert(interSection.NbPoints() == 1);

        // check if point is on the disk
        const auto p = OCCUtilities::point(interSection.Point(1));
        if (plane.contains(p, eps))
            return ResultType( p );
        else
            return ResultType( EmptyIntersection<3>() );
    }

    // The line is parallel. If the distance is > eps, there is no intersection
    const auto d = Vector<ctype, 3>(line.supportingPoint(),
                                    plane.projection(line.supportingPoint()));
    if (d.length() > eps)
        return ResultType( EmptyIntersection<3>() );

    // the intersection is the line
    return ResultType( line );
}

template<class ctype>
Intersection< Plane<ctype, 3>, Line<ctype, 3> >
intersect(const Plane<ctype, 3>& plane,
          const Line<ctype, 3>& line)
{ return intersect(plane, line, Precision<ctype>::confusion()); }

template<class ctype>
Intersection< Line<ctype, 3>, Plane<ctype, 3> >
intersect(const Line<ctype, 3>& line,
          const Plane<ctype, 3>& plane,
          ctype eps)
{ return intersect(plane, line, eps); }

template<class ctype>
Intersection< Line<ctype, 3>, Plane<ctype, 3> >
intersect(const Line<ctype, 3>& line,
          const Plane<ctype, 3>& plane)
{ return intersect(plane, line); }

template<class ctype>
Intersection< Disk<ctype>, Line<ctype, 3> >
intersect(const Disk<ctype>& disk,
          const Line<ctype, 3>& line,
          ctype eps)
{
    using ResultType = Intersection< Disk<ctype>, Line<ctype, 3> >;

    // intersect the line with the plane
    const auto planeLineIS = intersect(disk.supportingPlane(), line, eps);

    // empty result
    if (std::holds_alternative<EmptyIntersection<3>>(planeLineIS))
        return ResultType( EmptyIntersection<3>() );

    // potential point intersection
    else if (std::holds_alternative<Point<ctype, 3>>(planeLineIS))
    {
        // check if point is on the disk
        const auto p = std::get<Point<ctype, 3>>(planeLineIS);
        if (disk.contains(p, Precision<ctype>::confusion(), false))
            return ResultType( p );
        else
            return ResultType( EmptyIntersection<3>() );
    }

    // potential segment intersection
    else if (std::holds_alternative<Line<ctype, 3>>(planeLineIS))
    {
        // Intersection might be a segment, a touching point on the rim or empty
        // Use the BRep algorithm package to determine the part of the line on the ellipse.
        // For this, we create a segment of the line in the neighborhood of ellipse which is
        // large enough to cover the entire disk in case they intersect.
        auto source = line.projection(disk.center());
        auto target = source;
        auto dir = Vector<ctype, 3>(line.direction());
        dir *= disk.majorAxisLength()*20.0;
        source += dir;
        target -= dir;

        // find the part of the segment on the disk
        const auto segmentShape = OCCUtilities::makeEdge(source, target);
        const auto diskShape = OCCUtilities::getShape(disk);
        const auto isShape = OCCUtilities::intersect(segmentShape, diskShape, 0.1*eps);
        const auto isEdges = OCCUtilities::getEdges(isShape);

        assert(isEdges.size() <= 1);
        if (isEdges.size() == 1)
        {
            const auto vertices = OCCUtilities::getVertices(isEdges[0]);
            assert(vertices.size() == 2);
            return ResultType( Segment<ctype, 3>(OCCUtilities::point(vertices[0]),
                                                 OCCUtilities::point(vertices[1])) );
        }

        // the line can still touch the ellipse on the rim
        const auto cutShape = OCCUtilities::cut(segmentShape, diskShape, 0.1*eps);
        const auto vertices = OCCUtilities::getVertices(cutShape);

        assert(vertices.size() == 2 || vertices.size() == 4);
        if (vertices.size() == 2)
            return ResultType( EmptyIntersection<3>() );
        else
        {
            // find the new point
            for (const auto& v : vertices)
            {
                const auto curPoint = OCCUtilities::point(v);
                if (!curPoint.isEqual(source, eps) && !curPoint.isEqual(target, eps))
                    return ResultType( curPoint );
            }
        }

        throw std::runtime_error(std::string("Unexpected code behaviour"));
    }
    else
        throw std::runtime_error(std::string("Unexpected Plane-Line intersection result"));
}

template<class ctype>
Intersection< Disk<ctype>, Line<ctype, 3> >
intersect(const Disk<ctype>& disk,
          const Line<ctype, 3>& line)
{
    return intersect(disk,
                     line,
                     Precision<ctype>::confusion()*disk.minorAxisLength());
}

template<class ctype>
Intersection< Line<ctype, 3>, Disk<ctype> >
intersect(const Line<ctype, 3>& line,
          const Disk<ctype>& disk,
          ctype eps)
{ return intersect(disk, line, eps); }

template<class ctype>
Intersection< Line<ctype, 3>, Disk<ctype> >
intersect(const Line<ctype, 3>& line,
          const Disk<ctype>& disk)
{ return intersect(disk, line); }

template<class ctype>
Intersection< Disk<ctype>, Disk<ctype> >
intersect(const Disk<ctype>& disk1,
          const Disk<ctype>& disk2,
          ctype eps)
{
    using ResultType  = Intersection< Disk<ctype>, Disk<ctype> >;
    static constexpr int worldDim = 3;

    // first intersect the supporting planes
    std::string isGeomName;
    const auto planeIS = intersect(disk1.supportingPlane(), disk2.supportingPlane());
    std::visit([&] (auto&& is) { isGeomName = is.name(); }, planeIS);

    // if the planes don't intersect, the disks don't either
    if (std::holds_alternative<EmptyIntersection<worldDim>>(planeIS))
        return ResultType( EmptyIntersection<worldDim>() );

    // if the result is a line, find the possible segment or point intersection
    else if (std::holds_alternative<Line<ctype, worldDim>>(planeIS))
    {
        const auto& isLine = std::get<Line<ctype, worldDim>>(planeIS);
        const auto is1 = intersect(disk1, isLine);
        const auto is2 = intersect(disk2, isLine);

        // both disks intersect the support plane of the other disks
        if (std::holds_alternative<Segment<ctype, worldDim>>(is1)
            && std::holds_alternative<Segment<ctype, worldDim>>(is2))
        {
            const auto& seg1 = std::get<Segment<ctype, worldDim>>(is1);
            const auto& seg2 = std::get<Segment<ctype, worldDim>>(is2);
            const auto segmentIS = intersect(seg1, seg2);

            if (std::holds_alternative<Segment<ctype, worldDim>>(segmentIS))
                return ResultType( std::get<Segment<ctype, worldDim>>(segmentIS) );
            if (std::holds_alternative<Point<ctype, worldDim>>(segmentIS))
                return ResultType( std::get<Point<ctype, worldDim>>(segmentIS) );
            else
                return ResultType( EmptyIntersection<worldDim>() );
        }

        // one of the disks might still touch the other disk
        if (std::holds_alternative<Segment<ctype, worldDim>>(is1)
            && std::holds_alternative<Point<ctype, worldDim>>(is2))
        {
            const auto& seg1 = std::get<Segment<ctype, worldDim>>(is1);
            const auto& p2 = std::get<Point<ctype, worldDim>>(is2);
            if (seg1.contains(p2, eps))
                return ResultType( p2 );
        }
        else if (std::holds_alternative<Point<ctype, worldDim>>(is1)
            && std::holds_alternative<Segment<ctype, worldDim>>(is2))
        {
            const auto& p1 = std::get<Point<ctype, worldDim>>(is1);
            const auto& seg2 = std::get<Segment<ctype, worldDim>>(is2);
            if (seg2.contains(p1, eps))
                return ResultType( p1 );
        }

        // intersection is empty
        return ResultType( EmptyIntersection<worldDim>() );
    }

    // if the result is a plane, find the intersection surface
    else if (std::holds_alternative<Plane<ctype, worldDim>>(planeIS))
        throw std::runtime_error(std::string("NotImplemented: planar disk-disk intersections"));

    throw std::runtime_error(std::string("Unexpected plane-plane intersection result"));
}

template<class ctype>
Intersection< Disk<ctype>, Disk<ctype> >
intersect(const Disk<ctype>& disk1,
          const Disk<ctype>& disk2)
{
    using std::min;
    auto eps = min(disk1.minorAxisLength(), disk2.minorAxisLength());
    eps *= Precision<ctype>::confusion();
    return intersect(disk1, disk2, eps);
}

template<class ctype>
Intersection< CylindricalSurface<ctype>, Disk<ctype> >
intersect(const CylindricalSurface<ctype>& cylSurface,
          const Disk<ctype>& disk,
          ctype eps)
{
    // possible result geometries
    using ResultType = Intersection< CylindricalSurface<ctype>, Disk<ctype> >;
    using Point = Frackit::Point<ctype, 3>;
    using Segment = Frackit::Segment<ctype, 3>;
    using EllipseArc = Frackit::EllipseArc<ctype, 3>;
    using Ellipse = Frackit::Ellipse<ctype, 3>;
    using Direction = typename Ellipse::Direction;
    using Vector = typename Direction::Vector;

    // Intersect the disk with the cylinder
    const auto diskFace = OCCUtilities::getShape(disk);
    const auto cylinder = OCCUtilities::getShape(cylSurface.cylinder());
    const auto cylLateralFace = OCCUtilities::getShape(cylSurface);

    const auto containedFaceShape = OCCUtilities::intersect(diskFace, cylinder, 0.1*eps);
    auto containedFaces = OCCUtilities::getFaces(containedFaceShape);
    assert(containedFaces.size() <= 1);

    // intersect the disk boundary with the cylinder surface (detect touching points)
    const auto ellipseShape = OCCUtilities::getShape(disk.boundingEllipse());
    const auto ellipseCut = OCCUtilities::cut(ellipseShape, cylLateralFace, 0.1*eps);
    const auto ellipseCutVertices = OCCUtilities::getVertices(ellipseCut);

    // only touching points are possible
    if (containedFaces.size() == 0)
    {
        // find vertices on surface (avoid duplicates)
        std::vector<Point> touchingPoints;
        for (const auto& v : ellipseCutVertices)
        {
            auto p = OCCUtilities::point(v);
            if (cylSurface.contains(p, eps)
                && std::none_of(touchingPoints.begin(),
                                touchingPoints.end(),
                                [&p, eps] (const auto& pnt) { return pnt.isEqual(p, eps); }))
                touchingPoints.emplace_back(std::move(p));
        }

        return ResultType({touchingPoints.begin(), touchingPoints.end()});
    }

    // there are intersection edges, get orientation of the geometries
    gp_Vec dn(OCCUtilities::direction(disk.normal()));
    gp_Vec ca(OCCUtilities::direction(cylSurface.direction()));
    const bool diskIsParallel = dn.IsNormal(ca, Precision<ctype>::angular());
    const bool diskIsOrthogonal = dn.IsParallel(ca, Precision<ctype>::angular());

    // ... and (maybe) the ellipse on which the results live
    Ellipse infEllipse;
    if (diskIsOrthogonal)
        infEllipse = Ellipse(disk.supportingPlane().projection(cylSurface.centerSegment().source()),
                             cylSurface.base1(), cylSurface.base2(),
                             cylSurface.radius(), cylSurface.radius());
    else if (!diskIsParallel)
    {
        Vector cylDir(cylSurface.direction());
        cylDir *= disk.majorAxisLength();

        const auto& diskPlane = disk.supportingPlane();
        Direction majAxis(Vector(disk.center(), diskPlane.projection(disk.center() + cylDir)));
        Direction minAxis(crossProduct(Vector(disk.normal()), Vector(majAxis)));
        if (!isRightHandSystem(Vector(majAxis), Vector(minAxis), Vector(disk.normal())))
            majAxis.invert();

        using std::cos;
        const ctype majAxisLength = cylSurface.radius()/cos(dn.Angle(ca));
        const auto& cylAxisLine = cylSurface.centerSegment().supportingLine();
        const auto center = std::get<Point>(intersect(diskPlane, cylAxisLine, eps));
        infEllipse = Ellipse(center, majAxis, minAxis, majAxisLength, cylSurface.radius());
    }

    // intersect the wire of the contained disk with cylinder surface
    const auto containedFaceWires = OCCUtilities::getWires(containedFaceShape);
    assert(containedFaceWires.size() == 1);

    const auto wireIntersection = OCCUtilities::intersect(containedFaceWires[0], cylLateralFace, 0.1*eps);
    const auto wiresOnCylSurface = OCCUtilities::getWires(wireIntersection);

    // parse each wire into its basic geometry
    std::vector<EllipseArc> isArcs;
    std::vector<Segment> isSegments;
    for (const auto& wire : wiresOnCylSurface)
    {
        const auto wireEdges = OCCUtilities::getEdges(wire);
        const auto corners = OCCUtilities::getBoundaryVertices(wireEdges);

        const auto& c1 = corners.first[1] == 0 ? TopExp::FirstVertex(wireEdges[corners.first[0]])
                                               : TopExp::LastVertex(wireEdges[corners.first[0]]);
        const auto& c2 = corners.second[1] == 0 ? TopExp::FirstVertex(wireEdges[corners.second[0]])
                                                : TopExp::LastVertex(wireEdges[corners.second[0]]);
        const auto p1 = OCCUtilities::point(c1);
        const auto p2 = OCCUtilities::point(c2);

        // if both corners are the same the full ellipse is the intersection
        if (p1.isEqual(p2, eps))
            return ResultType({ infEllipse });

        if (!diskIsParallel)
        {
            // get corners of the wire
            EllipseArc arc1(infEllipse, p1, p2);
            EllipseArc arc2(infEllipse, p2, p1);
            const auto center1 = OCCUtilities::point(arc1.getPoint(0.5));
            const auto center2 = OCCUtilities::point(arc2.getPoint(0.5));

            // select the arc whose center is on the set of given edges
            unsigned int resultArcIndex = 0;
            for (const auto& edge : wireEdges)
            {
                // get unbounded curve and parameter bounds (umin, umax), then trim
                Standard_Real uMin, uMax;
                Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, uMin, uMax);
                curve = new Geom_TrimmedCurve(curve, uMin, uMax);
                GeomAPI_ProjectPointOnCurve c1OnCurve(center1, curve);
                GeomAPI_ProjectPointOnCurve c2OnCurve(center2, curve);
                if (c1OnCurve.LowerDistance() < eps) { resultArcIndex = 1; break; }
                if (c2OnCurve.LowerDistance() < eps) { resultArcIndex = 2; break; }
            }

            assert(resultArcIndex != 0);
            if (resultArcIndex == 1) isArcs.emplace_back( std::move(arc1) );
            else isArcs.emplace_back( std::move(arc2) );
        }
        else
            isSegments.emplace_back(p1, p2);
    }

    // There might still be intersection points (avoid duplicates!)
    std::vector<Point> isPoints;
    if (wiresOnCylSurface.size() < 2)
    {
        for (const auto& v : ellipseCutVertices)
        {
            auto p = OCCUtilities::point(v);
            if (cylSurface.contains(p, eps))
            {
                const auto onSeg = std::any_of(isSegments.begin(),
                                               isSegments.end(),
                                               [&p, eps] (const auto& seg) { return seg.contains(p, eps); });
                const auto onArc = std::any_of(isArcs.begin(),
                                               isArcs.end(),
                                               [&p, eps] (const auto& arc) { return arc.contains(p, eps); });
                if (!onSeg && !onArc && std::none_of(isPoints.begin(),
                                                     isPoints.end(),
                                                     [&p, eps] (const auto isP) { return isP.isEqual(p, eps); }))
                        isPoints.emplace_back(std::move(p));
            }
        }
    }

    ResultType result;
    for (auto&& p : isPoints)
        result.emplace_back(std::move(p));
    for (auto&& s : isSegments)
        result.emplace_back(std::move(s));
    for (auto&& arc : isArcs)
        result.emplace_back(std::move(arc));

    return result.empty() ? ResultType({ EmptyIntersection<3>() }) : result;
}

template<class ctype>
Intersection< CylindricalSurface<ctype>, Disk<ctype> >
intersect(const CylindricalSurface<ctype>& cylSurface, const Disk<ctype>& disk)
{
    using std::min;
    auto eps = min(disk.minorAxisLength(), cylSurface.radius());
    eps = min(eps, cylSurface.height());
    eps *= Precision<ctype>::confusion();
    return intersect(cylSurface, disk, eps);
}

template<class ctype>
Intersection< Disk<ctype>, CylindricalSurface<ctype> >
intersect(const Disk<ctype>& disk, const CylindricalSurface<ctype>& cylSurface, ctype eps)
{ return intersect(cylSurface, disk, eps); }

template<class ctype>
Intersection< Disk<ctype>, CylindricalSurface<ctype> >
intersect(const Disk<ctype>& disk, const CylindricalSurface<ctype>& cylSurface)
{
    using std::min;
    auto eps = min(disk.minorAxisLength(), cylSurface.radius());
    eps = min(eps, cylSurface.height());
    eps *= Precision<ctype>::confusion();
    return intersect(cylSurface, disk, eps);
}

} // end namespace Frackit

#endif // FRACKIT_INTERSECT_HH
