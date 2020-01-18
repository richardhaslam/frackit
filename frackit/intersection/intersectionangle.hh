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

#include <gp_Pnt2d.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <TopoDS_Face.hxx>
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
#include <frackit/geometry/name.hh>

#include <frackit/occ/gputilities.hh>
#include <frackit/occ/geomutilities.hh>
#include <frackit/occ/breputilities.hh>

#include "emptyintersection.hh"
#include "intersect.hh"

namespace Frackit {

/*!
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
                      const IsGeometry& isGeom)
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
                      const Line<ctype, wd>& isLine)
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
                      const Plane<ctype, wd>& plane2)
    {
        assert( !isEmptyIntersection(intersect(plane1, plane2)) );

        using std::abs;
        using std::acos;
        using Vector = Vector<ctype, wd>;
        return acos( abs(Vector(plane1.normal())*Vector(plane2.normal())) );
    }

    /*!
     * \brief Returns the angle in which two disks touch in a point
     * \param disk1 The first plane
     * \param disk2 The second plane
     * \param isPoint The touching point
     */
    ctype operator() (const Disk<ctype>& disk1,
                      const Disk<ctype>& disk2,
                      const Point<ctype, 3>& isPoint)
    { return (*this)(disk1.supportingPlane(), disk2.supportingPlane()); }

    /*!
     * \brief Returns the angle in which two disks intersect in a segment
     * \param disk1 The first plane
     * \param disk2 The second plane
     * \param isSeg The intersection segment
     */
    ctype operator() (const Disk<ctype>& disk1,
                      const Disk<ctype>& disk2,
                      const Segment<ctype, 3>& isSeg)
    { return (*this)(disk1.supportingPlane(), disk2.supportingPlane()); }

    /*!
     * \brief Returns the angle in which a disk and a
     *        cylinder surface touch in a point.
     * \param disk The disk
     * \param cylSurface The cylinder surface
     * \param isPoint The touching point
     */
    ctype operator() (const Disk<ctype>& disk,
                      const CylinderSurface<ctype>& cylSurface,
                      const Point<ctype, 3>& isPoint)
    { return (*this)(disk.supportingPlane(), cylSurface.getTangentPlane(isPoint)); }

    /*!
     * \brief Returns the angle in which a cylinder
     *        surface and a disk touch in a point.
     * \param cylSurface The cylinder surface
     * \param disk The disk
     * \param isPoint The touching point
     */
    ctype operator() (const CylinderSurface<ctype>& cylSurface,
                      const Disk<ctype>& disk,
                      const Point<ctype, 3>& isPoint)
    { return(*this)(disk, cylSurface, isPoint); }

    /*!
     * \brief Returns the angle in which a disk
     *        and a cylinder surface intersect in a segment.
     * \param disk The disk
     * \param cylSurface The cylinder surface
     * \param isSeg The intersection segment
     * \note An intersection segment on the cylinder surface can only
     *       occur if the intersection is parallel to the cylinder axis.
     *       Thus, the tangent plane on the cylinder is the same along
     *       the intersection segment, and we compute the angle for an
     *       arbitrary point on the segment; here: the first corner.
     */
    ctype operator() (const Disk<ctype>& disk,
                      const CylinderSurface<ctype>& cylSurface,
                      const Segment<ctype, 3>& isSeg)
    { return (*this)(disk.supportingPlane(), cylSurface.getTangentPlane(isSeg.source())); }

    /*!
     * \brief Returns the angle in which a cylinder surface
     *        and a disk intersect in a segment.
     * \param disk The disk
     * \param cylSurface The cylinder surface
     * \param isSeg The intersection segment
     * \note An intersection segment on the cylinder surface can only
     *       occur if the intersection is parallel to the cylinder axis.
     *       Thus, the tangent plane on the cylinder is the same along
     *       the intersection segment, and we compute the angle for an
     *       arbitrary point on the segment; here: the first corner.
     */
    ctype operator() (const CylinderSurface<ctype>& cylSurface,
                      const Disk<ctype>& disk,
                      const Segment<ctype, 3>& isSeg)
    { return (*this)(disk, cylSurface, isSeg); }

    /*!
     * \brief Returns the angle in which a cylinder surface
     *        and a disk intersect in an ellipse arc.
     * \param disk The disk
     * \param cylSurface The cylinder surface
     * \param isArc The intersection arc
     */
    ctype operator() (const Disk<ctype>& disk,
                      const CylinderSurface<ctype>& cylSurface,
                      const EllipseArc<ctype, 3>& isArc)
    {
        // use the minimum angle between the disk plane and the
        // tangent plane on the surface at four sample points
        std::array<ctype, 5> params({0.0, 0.25, 0.5, 0.75, 1.0});
        ctype resultAngle = std::numeric_limits<ctype>::max();
        const auto& diskPlane = disk.supportingPlane();

        using std::min;
        for (auto param : params)
        {
            const auto p = isArc.getPoint(param);
            const auto tangentPlane = cylSurface.getTangentPlane(p);
            resultAngle = min(resultAngle, (*this)(diskPlane, tangentPlane));
        }

        return resultAngle;
    }

    /*!
     * \brief Returns the angle in which a disk
     *        and a cylinder surface intersect in an ellipse arc.
     * \param cylSurface The cylinder surface
     * \param disk The disk
     * \param isArc The intersection arc
     */
    ctype operator() (const CylinderSurface<ctype>& cylSurface,
                      const Disk<ctype>& disk,
                      const EllipseArc<ctype, 3>& isArc)
    { return (*this)(disk, cylSurface, isArc); }

    /*!
     * \brief Returns the angle in which a disk
     *        and a cylinder surface intersect in an ellipse.
     * \param cylSurface The cylinder surface
     * \param disk The disk
     * \param isEllipse The intersection ellipse
     */
    ctype operator() (const Disk<ctype>& disk,
                      const CylinderSurface<ctype>& cylSurface,
                      const Ellipse<ctype, 3>& isEllipse)
    {
        // use the minimum angle between the disk plane and the
        // tangent plane on the surface at eight sample points
        const auto& diskPlane = disk.supportingPlane();
        ctype resultAngle = std::numeric_limits<ctype>::max();
        std::array<ctype, 9> params({0.0, 0.125, 0.25, 0.375,
                                     0.5, 0.625, 0.75, 0.875,
                                     1.0});

        using std::min;
        for (auto param : params)
        {
            const auto p = isEllipse.getPoint(param);
            const auto tangentPlane = cylSurface.getTangentPlane(p);
            resultAngle = min(resultAngle, (*this)(diskPlane, tangentPlane));
        }

        return resultAngle;
    }

    /*!
     * \brief Returns the angle in which a cylinder surface
     *        and a disk intersect in an ellipse.
     * \param disk The disk
     * \param cylSurface The cylinder surface
     * \param isEllipse The intersection ellipse
     */
    ctype operator() (const CylinderSurface<ctype>& cylSurface,
                      const Disk<ctype>& disk,
                      const Ellipse<ctype, 3>& isEllipse)
    { return (*this)(disk, cylSurface, isEllipse); }

    /*!
     * \brief Returns the angle in which a disk
     *        and a face shape touch in a point.
     * \param disk The disk
     * \param face The face shape
     * \param isPoint The touching point
     */
    ctype operator() (const Disk<ctype>& disk,
                      const TopoDS_Face& face,
                      const Point<ctype, 3>& isPoint)
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
        return (*this)(disk.supportingPlane(), tangentPlane);
    }

    /*!
     * \brief Returns the angle in which a face shape
     *        and a disk intersect in a point.
     * \param face The face shape
     * \param disk The disk
     * \param isPoint The touching point
     */
    ctype operator() (const TopoDS_Face& face,
                      const Disk<ctype>& disk,
                      const Point<ctype, 3>& isPoint)
    { return (*this)(disk, face, isPoint); }

    /*!
     * \brief Returns the angle in which a disk
     *        and a face shape intersect in an edge.
     * \param disk The disk
     * \param face The face shape
     * \param isEdge The intersection edge
     */
    ctype operator() (const Disk<ctype>& disk,
                      const TopoDS_Face& face,
                      const TopoDS_Edge& isEdge)
    {
        // compute the angle at several sample points along the edge and take minimum
        const auto edgeHandle = OCCUtilities::getGeomHandle(isEdge);
        const auto deltaParam = edgeHandle->LastParameter() - edgeHandle->FirstParameter();

        ctype resultAngle = std::numeric_limits<ctype>::max();
        std::array<ctype, 5> paramFactors({0.0, 0.25, 0.5, 0.75, 1.0});

        using std::min;
        for (auto f : paramFactors)
        {
            const auto param = edgeHandle->FirstParameter() + f*deltaParam;
            const auto isPoint = OCCUtilities::point(edgeHandle->Value(param));
            resultAngle = min(resultAngle, (*this)(disk, face, isPoint));
        }

        return resultAngle;
    }

    /*!
     * \brief Returns the angle in which a face shape
     *        and a disk intersect in an edge.
     * \param face The face shape
     * \param disk The disk
     * \param isEdge The intersection edge
     */
    ctype operator() (const TopoDS_Face& face,
                      const Disk<ctype>& disk,
                      const TopoDS_Edge& isEdge)
    { return (*this)(disk, face, isEdge); }

    /*!
     * \brief Returns the angle in which a disk intersects
     *        a face shape in a face.
     * \param disk The disk
     * \param face The face shape
     * \param isFace The intersection face
     */
    ctype operator() (const Disk<ctype>& disk,
                      const TopoDS_Face& face,
                      const TopoDS_Face& isFace)
    { return 0.0; }

    /*!
     * \brief Returns the angle in which a face shape
     *        intersects a disk in a face.
     * \param face The face shape
     * \param disk The disk
     * \param isFace The intersection face
     */
    ctype operator() (const TopoDS_Face& face,
                      const Disk<ctype>& disk,
                      const TopoDS_Face& isFace)
    { return (*this)(disk, face, isFace); }

    /*!
     * \brief Overload for an intersection variant
     * \param geo1 The first geometry
     * \param geo2 The second geometry
     * \param intersection The intersection variant
     */
    template<class Geo1, class Geo2, class... T>
    ctype operator() (const Geo1& geo1,
                      const Geo2& geo2,
                      const std::variant<T...>& intersection)
    { return std::visit([&] (auto&& is) { return (*this)(geo1, geo2, is); }, intersection); }

    /*!
     * \brief Overload for general intersections
     *        containing possibly various types.
     * \note We take the minimum ocurring angle here
     */
    template<class Geo1, class Geo2, class T>
    ctype operator() (const Geo1& geo1,
                      const Geo2& geo2,
                      const std::vector<T>& intersections)
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
