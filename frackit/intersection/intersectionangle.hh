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
 * \ingroup Intersection
 * \brief Class that contains functionality to determine
 *        the angle in which two geometries intersect.
 */
#ifndef FRACKIT_INTERSECTION_ANGLE_HH
#define FRACKIT_INTERSECTION_ANGLE_HH

#include <cmath>
#include <algorithm>
#include <cassert>
#include <variant>
#include <array>
#include <vector>
#include <stdexcept>
#include <limits>
#include <type_traits>

#include <gp_Pnt2d.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>

#include <frackit/precision/defaultepsilon.hh>

#include <frackit/geometry/point.hh>
#include <frackit/geometry/segment.hh>
#include <frackit/geometry/ellipse.hh>
#include <frackit/geometry/ellipsearc.hh>
#include <frackit/geometry/plane.hh>
#include <frackit/geometry/disk.hh>
#include <frackit/geometry/line.hh>
#include <frackit/geometry/cylindersurface.hh>

#include <frackit/geometryutilities/name.hh>
#include <frackit/geometryutilities/isplanar.hh>
#include <frackit/common/extractdimension.hh>

#include <frackit/occ/gputilities.hh>
#include <frackit/occ/geomutilities.hh>
#include <frackit/occ/breputilities.hh>

#include "emptyintersection.hh"
#include "intersect.hh"

namespace Frackit {

/*!
 * \ingroup Intersection
 * \brief Class that contains functions to compute the
 *        intersection angle between two geometries,
 *        intersecting in a given intersection geometry.
 */
template<class ctype = double>
class IntersectionAngle
{

public:

    /*!
     * \brief Returns the angle in which two geometries intersect
     * \param geo1 The first geometry
     * \param geo2 The second geometry
     * \param isGeom The geometry of the intersection between geo1 and geo2
     * \note This overload is active when no specialization is available
     */
    template<class Geo1, class Geo2, class IsGeometry>
    ctype operator() (const Geo1& geo1,
                      const Geo2& geo2,
                      const IsGeometry& isGeom) const
    {
        std::string msg = "Intersection angle not implemented for ";
        msg += "\"" + geometryName(geo1) + "\"";
        msg += " and ";
        msg += "\"" + geometryName(geo2) + "\" and the intersection geometry ";
        msg += "\"" + geometryName(isGeom) + "\"";
        throw std::runtime_error( msg );
    }

    /*!
     * \brief Returns the angle in which two planes intersect in a line
     * \param plane1 The first plane
     * \param plane2 The second plane
     * \param isLine The line in which the two planes intersect
     */
    template<int wd>
    ctype operator() (const Plane<ctype, wd>& plane1,
                      const Plane<ctype, wd>& plane2,
                      const Line<ctype, wd>& isLine) const
    {
        using std::abs;
        using std::acos;
        using Vector = Vector<ctype, wd>;
        return acos( abs(Vector(plane1.normal())*Vector(plane2.normal())) );
    }

    /*!
     * \brief Returns the angle in which two planes intersect
     * \param plane1 The first plane
     * \param plane2 The second plane
     * \note This overload is used when one knows that the planes
     *       intersect but does not want to compute the intersection.
     */
    template<int wd>
    ctype operator() (const Plane<ctype, wd>& plane1,
                      const Plane<ctype, wd>& plane2) const
    {
        assert( !isEmptyIntersection(intersect(plane1, plane2)) );

        using std::abs;
        using std::acos;
        using Vector = Vector<ctype, wd>;
        return acos( abs(Vector(plane1.normal())*Vector(plane2.normal())) );
    }

    /*!
     * \brief Returns the angle in which two planar 2-dimensional
     *        geometries embedded in 3d space touch in a point
     * \param geo1 The first planar geometry
     * \param geo2 The first planar geometry
     * \param isPoint The touching point
     */
    template<class Geo1, class Geo2,
             std::enable_if_t<IsPlanarGeometry<Geo1>::value
                              && IsPlanarGeometry<Geo2>::value, int> = 0>
    ctype operator() (const Geo1& geo1,
                      const Geo2& geo2,
                      const Point<ctype, 3>& isPoint) const
    { return (*this)(geo1.supportingPlane(), geo2.supportingPlane()); }

    /*!
     * \brief Returns the angle in which two planar 2-dimensional
     *        geometries embedded in 3d space touch intersect in a segment
     * \param geo1 The first planar geometry
     * \param geo2 The first planar geometry
     * \param isSeg The intersection segment
     */
    template<class Geo1, class Geo2,
             std::enable_if_t<IsPlanarGeometry<Geo1>::value
                              && IsPlanarGeometry<Geo2>::value, int> = 0>
    ctype operator() (const Geo1& geo1,
                      const Geo2& geo2,
                      const Segment<ctype, 3>& isSeg) const
    { return (*this)(geo1.supportingPlane(), geo2.supportingPlane()); }

    /*!
     * \brief Returns the angle in which a planar
     *        2-dimensional geometry in 2d space
     *        and a cylinder surface touch in a point.
     * \param geo The planar geometry
     * \param cylSurface The cylinder surface
     * \param isPoint The touching point
     */
    template<class Geo, std::enable_if_t<IsPlanarGeometry<Geo>::value, int> = 0>
    ctype operator() (const Geo& geo,
                      const CylinderSurface<ctype>& cylSurface,
                      const Point<ctype, 3>& isPoint) const
    { return (*this)(geo.supportingPlane(), cylSurface.getTangentPlane(isPoint)); }

    /*!
     * \brief Returns the angle in which a cylinder
     *        surface and planar 2-dimensional geometry
     *        in 2d space touch in a point.
     * \param cylSurface The cylinder surface
     * \param geo The planar geometry
     * \param isPoint The touching point
     */
    template<class Geo, std::enable_if_t<IsPlanarGeometry<Geo>::value, int> = 0>
    ctype operator() (const CylinderSurface<ctype>& cylSurface,
                      const Geo& geo,
                      const Point<ctype, 3>& isPoint) const
    { return(*this)(geo, cylSurface, isPoint); }

    /*!
     * \brief Returns the angle in which a planar 2-dimensional geometry
     *        and a cylinder surface intersect in a segment.
     * \param geo The planar geometry
     * \param cylSurface The cylinder surface
     * \param isSeg The intersection segment
     * \note An intersection segment on the cylinder surface can only
     *       occur if the intersection is parallel to the cylinder axis.
     *       Thus, the tangent plane on the cylinder is the same along
     *       the intersection segment, and we compute the angle for an
     *       arbitrary point on the segment; here: the first corner.
     */
    template<class Geo, std::enable_if_t<IsPlanarGeometry<Geo>::value, int> = 0>
    ctype operator() (const Geo& geo,
                      const CylinderSurface<ctype>& cylSurface,
                      const Segment<ctype, 3>& isSeg) const
    { return (*this)(geo.supportingPlane(), cylSurface.getTangentPlane(isSeg.source())); }

    /*!
     * \brief Returns the angle in which a cylinder surface and
     *        a planar 2-dimensional geometry intersect in a segment.
     * \param cylSurface The cylinder surface
     * \param geo The planar geometry
     * \param isSeg The intersection segment
     * \note An intersection segment on the cylinder surface can only
     *       occur if the intersection is parallel to the cylinder axis.
     *       Thus, the tangent plane on the cylinder is the same along
     *       the intersection segment, and we compute the angle for an
     *       arbitrary point on the segment; here: the first corner.
     */
    template<class Geo, std::enable_if_t<IsPlanarGeometry<Geo>::value, int> = 0>
    ctype operator() (const CylinderSurface<ctype>& cylSurface,
                      const Geo& geo,
                      const Segment<ctype, 3>& isSeg) const
    { return (*this)(geo, cylSurface, isSeg); }

    /*!
     * \brief Returns the angle in which a planar 2-dimensional
     *        geometry and a cylinder surface intersect in an ellipse arc.
     * \param geo The planar geometry
     * \param cylSurface The cylinder surface
     * \param isArc The intersection arc
     */
    template<class Geo, std::enable_if_t<IsPlanarGeometry<Geo>::value, int> = 0>
    ctype operator() (const Geo& geo,
                      const CylinderSurface<ctype>& cylSurface,
                      const EllipseArc<ctype, 3>& isArc) const
    {
        // use the minimum angle between the geometry plane and the
        // tangent plane on the surface at four sample points
        constexpr std::array<ctype, 5> params({0.01, 0.25, 0.5, 0.75, 0.99});
        ctype resultAngle = std::numeric_limits<ctype>::max();
        const auto& geoPlane = geo.supportingPlane();

        using std::min;
        for (auto param : params)
        {
            const auto p = isArc.getPoint(param);
            const auto tangentPlane = cylSurface.getTangentPlane(p);
            resultAngle = min(resultAngle, (*this)(geoPlane, tangentPlane));
        }

        return resultAngle;
    }

    /*!
     * \brief Returns the angle in which a cyliner surface and
     *        a planar 2-dimensional geometry intersect in an ellipse arc.
     * \param cylSurface The cylinder surface
     * \param geo The planar geometry
     * \param isArc The intersection arc
     */
    template<class Geo, std::enable_if_t<IsPlanarGeometry<Geo>::value, int> = 0>
    ctype operator() (const CylinderSurface<ctype>& cylSurface,
                      const Geo& geo,
                      const EllipseArc<ctype, 3>& isArc) const
    { return (*this)(geo, cylSurface, isArc); }

    /*!
     * \brief Returns the angle in which a planar 2-dimensional geometry
     *        and a cylinder surface intersect in an ellipse.
     * \param cylSurface The cylinder surface
     * \param geo The planar geometry
     * \param isEllipse The intersection ellipse
     */
    template<class Geo, std::enable_if_t<IsPlanarGeometry<Geo>::value, int> = 0>
    ctype operator() (const Geo& geo,
                      const CylinderSurface<ctype>& cylSurface,
                      const Ellipse<ctype, 3>& isEllipse) const
    {
        // use the minimum angle between the geometry plane and the
        // tangent plane on the surface at eight sample points
        const auto& geoPlane = geo.supportingPlane();
        ctype resultAngle = std::numeric_limits<ctype>::max();
        constexpr std::array<ctype, 9> params({0.001, 0.125, 0.25, 0.375,
                                              0.5,   0.625, 0.75, 0.875,
                                              0.999});

        using std::min;
        for (auto param : params)
        {
            const auto p = isEllipse.getPoint(param);
            const auto tangentPlane = cylSurface.getTangentPlane(p);
            resultAngle = min(resultAngle, (*this)(geoPlane, tangentPlane));
        }

        return resultAngle;
    }

    /*!
     * \brief Returns the angle in which a cylinder surface
     *        and a planar 2-dimensional geometry intersect in an ellipse.
     * \param geo The planar geometry
     * \param cylSurface The cylinder surface
     * \param isEllipse The intersection ellipse
     */
    template<class Geo, std::enable_if_t<IsPlanarGeometry<Geo>::value, int> = 0>
    ctype operator() (const CylinderSurface<ctype>& cylSurface,
                      const Geo& geo,
                      const Ellipse<ctype, 3>& isEllipse) const
    { return (*this)(geo, cylSurface, isEllipse); }

    /*!
     * \brief Returns the angle in which a planar 2d geometry
     *        and a face shape touch in a point.
     * \param geo The planar geometry
     * \param face The face shape
     * \param isPoint The touching point
     */
    template<class Geo, std::enable_if_t<IsPlanarGeometry<Geo>::value, int> = 0>
    ctype operator() (const Geo& geo,
                      const TopoDS_Face& face,
                      const Point<ctype, 3>& isPoint) const
    {
        // get the parameters of this point on the face via orthogonal projection
        const auto geomSurface = OCCUtilities::getGeomHandle(face);
        GeomAPI_ProjectPointOnSurf projection(OCCUtilities::point(isPoint), geomSurface);

        // since the intersection point should have been on the face, distance < eps!
        assert(projection.LowerDistance() < defaultEpsilon(face));

        ctype paramU, paramV;
        projection.LowerDistanceParameters(paramU, paramV);

        // construct the tangent plane of the face in the point
        gp_Pnt p;
        gp_Vec baseVec1, baseVec2;
        geomSurface->D1(projection, paramV, p, baseVec1, baseVec2);

        const auto base1 = OCCUtilities::vector(baseVec1);
        const auto base2 = OCCUtilities::vector(baseVec2);
        const Direction<ctype, 3> normal(crossProduct(base1, base2));
        const Plane<ctype, 3> tangentPlane(isPoint, normal);
        return (*this)(geo.supportingPlane(), tangentPlane);
    }

    /*!
     * \brief Returns the angle in which two face shapes touch in a point.
     * \param face1 The first face shape
     * \param face2 The second face shape
     * \param isPoint The touching point
     */
    ctype operator() (const TopoDS_Face& face1,
                      const TopoDS_Face& face2,
                      const Point<ctype, 3>& isPoint) const
    {
        // get the parameters of this point on the face via orthogonal projection
        const auto geomSurface1 = OCCUtilities::getGeomHandle(face1);
        const auto geomSurface2 = OCCUtilities::getGeomHandle(face2);
        GeomAPI_ProjectPointOnSurf projection1(OCCUtilities::point(isPoint), geomSurface1);
        GeomAPI_ProjectPointOnSurf projection2(OCCUtilities::point(isPoint), geomSurface2);

        // since the intersection point should have been on the face, distance < eps!
        assert(projection1.LowerDistance() < defaultEpsilon(face1));
        assert(projection2.LowerDistance() < defaultEpsilon(face2));

        ctype paramU1, paramV1, paramU2, paramV2;
        projection1.LowerDistanceParameters(paramU1, paramV1);
        projection2.LowerDistanceParameters(paramU2, paramV2);

        // construct the tangent planes of the faces in the point
        gp_Pnt p1, p2;
        gp_Vec baseVec11, baseVec12, baseVec21, baseVec22;
        geomSurface1->D1(projection1, paramV1, p1, baseVec11, baseVec12);
        geomSurface2->D1(projection2, paramV2, p2, baseVec21, baseVec22);

        using Direction = Direction<ctype, 3>;
        const auto base11 = OCCUtilities::vector(baseVec11);
        const auto base12 = OCCUtilities::vector(baseVec12);
        const Plane<ctype, 3> tangentPlane1(isPoint, Direction(crossProduct(base11, base12)));

        const auto base21 = OCCUtilities::vector(baseVec21);
        const auto base22 = OCCUtilities::vector(baseVec22);
        const Plane<ctype, 3> tangentPlane2(isPoint, Direction(crossProduct(base21, base22)));

        return (*this)(tangentPlane1, tangentPlane2);
    }

    /*!
     * \brief Returns the angle in which a face shape
     *        and a planar 2d geometry intersect in a point.
     * \param face The face shape
     * \param geo The planar geometry
     * \param isPoint The touching point
     */
    template<class Geo,
             std::enable_if_t<IsPlanarGeometry<Geo>::value, int> = 0>
    ctype operator() (const TopoDS_Face& face,
                      const Geo& geo,
                      const Point<ctype, 3>& isPoint) const
    { return (*this)(geo, face, isPoint); }

    /*!
     * \brief Returns the angle in which a planar 2d geometry
     *        and a face shape intersect in an edge.
     * \param geo The planar geometry
     * \param face The face shape
     * \param isEdge The intersection edge
     */
    template<class Geo, std::enable_if_t<IsPlanarGeometry<Geo>::value, int> = 0>
    ctype operator() (const Geo& geo,
                      const TopoDS_Face& face,
                      const TopoDS_Edge& isEdge) const
    {
        // compute the angle at several sample points along the edge and take minimum
        const auto edgeHandle = OCCUtilities::getGeomHandle(isEdge);
        const auto deltaParam = edgeHandle->LastParameter() - edgeHandle->FirstParameter();

        ctype resultAngle = std::numeric_limits<ctype>::max();
        constexpr std::array<ctype, 5> paramFactors({0.01, 0.25, 0.5, 0.75, 0.99});

        using std::min;
        for (auto f : paramFactors)
        {
            const auto param = edgeHandle->FirstParameter() + f*deltaParam;
            const auto isPoint = OCCUtilities::point(edgeHandle->Value(param));
            resultAngle = min(resultAngle, (*this)(geo, face, isPoint));
        }

        return resultAngle;
    }

    /*!
     * \brief Returns the angle in which two face shapes intersect in an edge.
     * \param face1 The first face shape
     * \param face2 The second face shape
     * \param isEdge The intersection edge
     */
    ctype operator() (const TopoDS_Face& face1,
                      const TopoDS_Face& face2,
                      const TopoDS_Edge& isEdge) const
    {
        // compute the angle at several sample points along the edge and take minimum
        const auto edgeHandle = OCCUtilities::getGeomHandle(isEdge);
        const auto deltaParam = edgeHandle->LastParameter() - edgeHandle->FirstParameter();

        ctype resultAngle = std::numeric_limits<ctype>::max();
        constexpr std::array<ctype, 5> paramFactors({0.01, 0.25, 0.5, 0.75, 0.99});

        using std::min;
        for (auto f : paramFactors)
        {
            const auto param = edgeHandle->FirstParameter() + f*deltaParam;
            const auto isPoint = OCCUtilities::point(edgeHandle->Value(param));
            resultAngle = min(resultAngle, (*this)(face1, face2, isPoint));
        }

        return resultAngle;
    }

    /*!
     * \brief Returns the angle in which a face shape
     *        and a planar 2d geometry intersect in an edge.
     * \param face The face shape
     * \param geo The planar geometry
     * \param isEdge The intersection edge
     */
    template<class Geo,
             std::enable_if_t<IsPlanarGeometry<Geo>::value
                              && !(std::is_same_v<Geo, TopoDS_Face>), int> = 0>
    ctype operator() (const TopoDS_Face& face,
                      const Geo& geo,
                      const TopoDS_Edge& isEdge) const
    { return (*this)(geo, face, isEdge); }

    /*!
     * \brief Returns the angle in which a geometry
     *        intersects a face shape in a face.
     * \param geo The planar geometry
     * \param face The face shape
     * \param isFace The intersection face
     * \note if the intersection object is a face itself,
     *       the angle in which the geometries intersect is always zero.
     */
    template<class Geo>
    ctype operator() (const Geo& geo,
                      const TopoDS_Face& face,
                      const TopoDS_Face& isFace) const
    { return 0.0; }

    /*!
     * \brief Returns the angle in which a face shape
     *        intersects a planar 2d geometry in a face.
     * \param face The face shape
     * \param geo The planar geometry
     * \param isFace The intersection face
     */
    template<class Geo,
             std::enable_if_t<IsPlanarGeometry<Geo>::value
                              && !(std::is_same_v<Geo, TopoDS_Face>), int> = 0>
    ctype operator() (const TopoDS_Face& face,
                      const Geo& geo,
                      const TopoDS_Face& isFace) const
    { return (*this)(geo, face, isFace); }

    /*!
     * \brief Overload for an intersection variant
     * \param geo1 The first geometry
     * \param geo2 The second geometry
     * \param intersection The intersection variant
     */
    template<class Geo1, class Geo2, class... T>
    ctype operator() (const Geo1& geo1,
                      const Geo2& geo2,
                      const std::variant<T...>& intersection) const
    { return std::visit([&] (auto&& is) { return (*this)(geo1, geo2, is); }, intersection); }

    /*!
     * \brief Overload for general intersections
     *        containing possibly various types.
     * \note We take the minimum ocurring angle here
     */
    template<class Geo1, class Geo2, class T>
    ctype operator() (const Geo1& geo1,
                      const Geo2& geo2,
                      const std::vector<T>& intersections) const
    {
        using std::min;
        ctype result = std::numeric_limits<ctype>::max();
        for (const auto& is : intersections)
            result = min(result, (*this)(geo1, geo2, is));

        return result;
    }
};

} // end namespace Frackit

#endif // FRACKIT_INTERSECTION_ANGLE_HH
