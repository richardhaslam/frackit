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

#include <frackit/common/utilities.hh>
#include <frackit/intersection/intersectiontraits.hh>
#include <frackit/intersection/emptyintersection.hh>

namespace Frackit {
namespace IntersectionAlgorithms {

//! Intersect a lateral cylinder surface and a disk
//! The result can be:
//! - an ellipse
//! - ellipse arc(s)
//! - segment(s)
//! - touching points
template<class ctype>
Intersection< CylindricalSurface<ctype>, Disk<ctype> >
intersect_cylinderSurface_disk(const CylindricalSurface<ctype>& cylSurface,
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

} // end namespace IntersectionAlgorithms
} // end namespace Frackit

#endif // FRACKIT_CYLINDERSURFACE_DISK_INTERSECTION_HH
