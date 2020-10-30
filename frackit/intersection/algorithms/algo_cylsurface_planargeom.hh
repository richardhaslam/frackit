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
/*!
 * \file
 * \brief Contains the intersection algorithm
 *        between a lateral cylinder surface
 *        and a planar (2d) geometry in 3d space.
 */
#ifndef FRACKIT_CYLINDERSURFACE_PLANAR_GEOMETRY_INTERSECTION_HH
#define FRACKIT_CYLINDERSURFACE_PLANAR_GEOMETRY_INTERSECTION_HH

#include <vector>
#include <algorithm>

#include <BRep_Tool.hxx>
#include <TopExp.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>

#include <frackit/occ/breputilities.hh>
#include <frackit/occ/geomutilities.hh>
#include <frackit/occ/gputilities.hh>

#include <frackit/geometry/cylindersurface.hh>
#include <frackit/common/extractdimension.hh>

#include <frackit/intersection/intersectiontraits.hh>
#include <frackit/intersection/emptyintersection.hh>

namespace Frackit {
namespace IntersectionAlgorithms {
namespace Detail {

    // gets the start and end point of the given wire
    template<class ctype>
    std::pair<Point<ctype, 3>, Point<ctype, 3>> getWireTips(const TopoDS_Wire& wire)
    {
        // find corners of the arc that this wire describes
        TopoDS_Vertex v1, v2;
        TopExp::Vertices(wire, v1, v2);
        assert(!v1.IsNull());
        assert(!v2.IsNull());
        assert(!v1.IsSame(v2));

        return std::make_pair(OCCUtilities::point(v1), OCCUtilities::point(v2));
    }

    // parses the wire into the ellipse arc living on the given ellipse
    template<class ctype>
    EllipseArc<ctype, 3> makeArc(const TopoDS_Wire& wire,
                                 const Ellipse<ctype, 3>& ellipse,
                                 const ctype eps)
    {
        using EllipseArc = EllipseArc<ctype, 3>;

        const auto [p1, p2] = getWireTips<ctype>(wire);

        // select the arc whose center is on the set of given edges
        EllipseArc arc1(ellipse, p1, p2);
        EllipseArc arc2(ellipse, p2, p1);
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
        return resultArcIndex == 1 ? arc1 : arc2;
    }

} // end namespace Detail

/*!
 * \brief Intersect a lateral cylinder surface and
 *        a planar geometry in three-dimensional space.
 *        The result can be:
 *        - an ellipse
 *        - ellipse arc(s)
 *        - segment(s)
 *        - touching points
 * \param cylSurface The lateral surface of a cylinder
 * \param faceGeom The planar geometry
 * \param charLength A characteristic length scale of faceGeom
 * \param eps Tolerance to be used for boolean operations
 */
template<class ctype, class PlanarGeometry>
Intersection< CylinderSurface<ctype>, PlanarGeometry >
intersect_cylinderSurface_planarGeometry(const CylinderSurface<ctype>& cylSurface,
                                         const PlanarGeometry& faceGeom,
                                         ctype charLength,
                                         ctype eps)
{
    static_assert( (DimensionalityTraits<PlanarGeometry>::geometryDimension() == 2 &&
                    DimensionalityTraits<PlanarGeometry>::worldDimension() == 3),
                  "This algorithm expects planar, two-dimensional geometries in 3d space");

    // possible result geometries
    using ResultType = Intersection< CylinderSurface<ctype>, PlanarGeometry >;
    using Point = Frackit::Point<ctype, 3>;
    using Segment = Frackit::Segment<ctype, 3>;
    using EllipseArc = Frackit::EllipseArc<ctype, 3>;
    using Ellipse = Frackit::Ellipse<ctype, 3>;
    using Direction = typename Ellipse::Direction;
    using Vector = typename Direction::Vector;

    // determine the orientation of the geometries
    const auto& faceGeomPlane = faceGeom.supportingPlane();
    gp_Vec dn(OCCUtilities::direction(faceGeomPlane.normal()));
    gp_Vec ca(OCCUtilities::direction(cylSurface.direction()));
    const bool faceGeomIsParallel = dn.IsNormal(ca, Precision<ctype>::angular());
    const bool faceGeomIsOrthogonal = dn.IsParallel(ca, Precision<ctype>::angular());
    const bool resultIsElliptical = faceGeomIsOrthogonal || !faceGeomIsParallel;

    // Depending on the orientation, the results might live on an ellipse.
    Ellipse infEllipse;
    if (faceGeomIsOrthogonal)
        infEllipse = Ellipse(faceGeomPlane.projection(cylSurface.centerSegment().source()),
                             cylSurface.base1(), cylSurface.base2(),
                             cylSurface.radius(), cylSurface.radius());
    else if (!faceGeomIsParallel)
    {
        Vector cylDir(cylSurface.direction());
        cylDir *= charLength;

        Direction majAxis(Vector(faceGeom.center(), faceGeomPlane.projection(faceGeom.center() + cylDir)));
        Direction minAxis(crossProduct(Vector(faceGeomPlane.normal()), Vector(majAxis)));
        if (!isRightHandSystem(Vector(majAxis), Vector(minAxis), Vector(faceGeomPlane.normal())))
            majAxis.invert();

        using std::cos;
        using std::abs;
        const ctype majAxisLength = abs(cylSurface.radius()/cos(dn.Angle(ca)));
        const auto& cylAxisLine = cylSurface.centerSegment().supportingLine();
        const auto center = std::get<Point>(intersect(faceGeomPlane, cylAxisLine, eps));
        infEllipse = Ellipse(center, majAxis, minAxis, majAxisLength, cylSurface.radius());
    }

    // Determine the part of the planar geometry inside the cylinder
    const auto faceShape = OCCUtilities::getShape(faceGeom);
    const auto cylinder = OCCUtilities::getShape(cylSurface.cylinder());
    const auto cylLateralFace = OCCUtilities::getShape(cylSurface);
    const auto containedFaceShape = OCCUtilities::intersect(faceShape, cylinder, 0.1*eps);
    const auto containedFaces = OCCUtilities::getFaces(containedFaceShape);
    assert(containedFaces.size() <= 1);

    // if the face is completely outside, there might still be touching edges or corners
    // first, check for wires on the cylinder surfaces
    const auto faceGeomWires = OCCUtilities::getWires(faceShape);
    assert(faceGeomWires.size() == 1);
    if (containedFaces.empty())
    {
        const auto is = OCCUtilities::intersect(faceGeomWires[0], cylLateralFace, 0.1*eps);
        const auto wires = OCCUtilities::getWires(is);

        if (!wires.empty())
        {
            ResultType result;
            if (resultIsElliptical)
            {
                for (const auto& w : wires)
                    result.emplace_back(Detail::makeArc(w, infEllipse, eps));
            }
            else
            {
                for (const auto& w : wires)
                {
                    const auto [p1, p2] = Detail::getWireTips<ctype>(w);
                    result.emplace_back(Segment(p1, p2));
                }
            }
            return result;
        }
    }

    // detect possible touching points
    const auto faceGeomWireCut = OCCUtilities::cut(faceGeomWires[0], cylLateralFace, eps*0.1);
    const auto faceGeomWireCutVertices = OCCUtilities::getVertices(faceGeomWireCut);

    std::vector<Point> touchCandidates;
    touchCandidates.reserve(faceGeomWireCutVertices.size());
    for (const auto& v : faceGeomWireCutVertices)
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

    // intersect the wire of the contained face shape with cylinder surface
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

        if (faceGeomIsOrthogonal || !faceGeomIsParallel)
            isArcs.emplace_back(Detail::makeArc(wire, infEllipse, eps));
        else
        {
            const auto [p1, p2] = Detail::getWireTips<ctype>(wire);
            isSegments.emplace_back(p1, p2);
        }
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

#endif // FRACKIT_CYLINDERSURFACE_PLANAR_GEOMETRY_INTERSECTION_HH
