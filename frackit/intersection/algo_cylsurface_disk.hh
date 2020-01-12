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
 * \brief Contains the intersection algorithm
 *        between a lateral cylinder surface and a disk.
 */
#ifndef FRACKIT_CYLINDERSURFACE_DISK_INTERSECTION_HH
#define FRACKIT_CYLINDERSURFACE_DISK_INTERSECTION_HH

#include <vector>
#include <algorithm>

#include <BRep_Tool.hxx>
#include <TopExp.hxx>
#include <TopoDS_Vertex.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>

#include <frackit/occ/breputilities.hh>
#include <frackit/occ/geomutilities.hh>
#include <frackit/occ/gputilities.hh>

#include <frackit/geometry/disk.hh>
#include <frackit/geometry/cylindersurface.hh>

#include "intersectiontraits.hh"
#include "emptyintersection.hh"

namespace Frackit {
namespace IntersectionAlgorithms {

//! Intersect a lateral cylinder surface and a disk
//! The result can be:
//! - an ellipse
//! - ellipse arc(s)
//! - segment(s)
//! - touching points
template<class ctype>
Intersection< CylinderSurface<ctype>, Disk<ctype> >
intersect_cylinderSurface_disk(const CylinderSurface<ctype>& cylSurface,
                               const Disk<ctype>& disk,
                               ctype eps)
{
    // possible result geometries
    using ResultType = Intersection< CylinderSurface<ctype>, Disk<ctype> >;
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

    // detect possible touching points
    const auto diskWires = OCCUtilities::getWires(diskFace);

    assert(diskWires.size() == 1);
    const auto diskWireCut = OCCUtilities::cut(diskWires[0], cylLateralFace, eps*0.1);
    const auto diskWireCutVertices = OCCUtilities::getVertices(diskWireCut);

    std::vector<Point> touchCandidates;
    touchCandidates.reserve(diskWireCutVertices.size());
    for (const auto& v : diskWireCutVertices)
    {
        const auto p = OCCUtilities::point(v);
        if (cylSurface.contains(p, eps))
            if (std::none_of(touchCandidates.begin(),
                             touchCandidates.end(),
                             [&p, eps] (const auto& tp) { return tp.isEqual(p, eps); }))
                touchCandidates.emplace_back(std::move(p));
    }

    // only touching points are possible
    if (containedFaces.empty())
    {
        if (touchCandidates.empty())
            return ResultType({ EmptyIntersection<3>() });

        ResultType result;
        for (auto& p : touchCandidates)
            result.emplace_back(std::move(p));
        return result;
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
        if (BRep_Tool::IsClosed(wire))
            return ResultType({ infEllipse });

        // find corners of the arc that this wire describes
        TopoDS_Vertex v1, v2;
        TopExp::Vertices(wire, v1, v2);
        assert(!v1.IsNull());
        assert(!v2.IsNull());
        assert(!v1.IsSame(v2));

        const auto p1 = OCCUtilities::point(v1);
        const auto p2 = OCCUtilities::point(v2);

        if (!diskIsParallel)
        {
            // select the arc whose center is on the set of given edges
            EllipseArc arc1(infEllipse, p1, p2);
            EllipseArc arc2(infEllipse, p2, p1);
            const auto center1 = OCCUtilities::point(arc1.getPoint(0.5));
            const auto center2 = OCCUtilities::point(arc2.getPoint(0.5));

            unsigned int resultArcIndex = 0;
            for (const auto& edge : OCCUtilities::getEdges(wire))
            {
                auto curve = OCCUtilities::getGeomHandle(edge);
                GeomAPI_ProjectPointOnCurve c1OnCurve(center1, curve);
                GeomAPI_ProjectPointOnCurve c2OnCurve(center2, curve);

                if (c1OnCurve.NbPoints() > 0 && c1OnCurve.LowerDistance() < eps) { resultArcIndex = 1; break; }
                if (c2OnCurve.NbPoints() > 0 && c2OnCurve.LowerDistance() < eps) { resultArcIndex = 2; break; }
            }

            assert(resultArcIndex != 0);
            if (resultArcIndex == 1) isArcs.emplace_back( std::move(arc1) );
            else isArcs.emplace_back( std::move(arc2) );
        }
        else
            isSegments.emplace_back(p1, p2);
    }

    // (Maybe) add touching points
    ResultType result;
    for (auto& p : touchCandidates)
    {
        const auto onSeg = std::any_of(isSegments.begin(),
                                       isSegments.end(),
                                       [&p, eps] (const auto& seg) { return seg.contains(p, eps); });
        const auto onArc = std::any_of(isArcs.begin(),
                                       isArcs.end(),
                                       [&p, eps] (const auto& arc) { return arc.contains(p, eps); });
        if (!onSeg && !onArc)
            result.emplace_back(std::move(p));
    }

    // add segments and arcs
    for (auto&& s : isSegments)
        result.emplace_back(std::move(s));
    for (auto&& arc : isArcs)
        result.emplace_back(std::move(arc));

    return result.empty() ? ResultType({ EmptyIntersection<3>() }) : result;
}

} // end namespace IntersectionAlgorithms
} // end namespace Frackit

#endif // FRACKIT_CYLINDERSURFACE_DISK_INTERSECTION_HH
