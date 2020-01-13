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
 * \brief Implementation details on distance to boundary constraints.
 */
#ifndef FRACKIT_ENTITYNETWORK_CONSTRAINT_IMPL_ADMISSIBLE_BOUNDARY_DISTANCE_HH
#define FRACKIT_ENTITYNETWORK_CONSTRAINT_IMPL_ADMISSIBLE_BOUNDARY_DISTANCE_HH

#include <stdexcept>
#include <variant>
#include <vector>

#include <TopExp.hxx>

#include <frackit/distance/distance.hh>
#include <frackit/distance/distancetoboundary.hh>
#include <frackit/distance/pointonboundary.hh>

#include <frackit/magnitude/magnitude.hh>
#include <frackit/intersection/intersectiontraits.hh>

#include <frackit/geometry/point.hh>
#include <frackit/geometry/segment.hh>
#include <frackit/geometry/ellipsearc.hh>
#include <frackit/geometry/ellipse.hh>
#include <frackit/geometry/name.hh>

namespace Frackit {
namespace ConstraintImpl {

    /*!
     * \brief Evaluates if the distance of an intersection geometry to
     *        the boundary of an entity is above a given threshold.
     */
    template<class IsGeom, class Geo, class ctype>
    bool isAdmissibleDistanceToBoundary(const IsGeom& is, const Geo& entity, ctype threshold)
    {
        std::string msg = "Distance to boundary not implemented for intersection geometry ";
        msg += "\"" + geometryName(is) + "\"";
        msg += " and entity geometry ";
        msg += "\"" + geometryName(entity) + "\"";
        throw std::runtime_error(std::string(msg));
    }

    /*!
     * \brief Overload for empty intersections. These fulfill the distance constraint.
     */
    template<int wd, class Geo, class ctype>
    bool isAdmissibleDistanceToBoundary(const EmptyIntersection<wd>& is,
                                        const Geo& entity,
                                        ctype threshold)
    { return true; }

    /*!
     * \brief Overload for point intersections.
     *        Points that lie on the boundary of the entity are considered
     *        to fulfill the distance constraint.
     */
    template<class ctype, int wd, class Geo, class ctype2>
    bool isAdmissibleDistanceToBoundary(const Point<ctype, wd>& is,
                                        const Geo& entity,
                                        ctype2 threshold)
    {
        if (pointOnGeometryBoundary(is, entity))
            return true;
        return computeDistanceToBoundary(is, entity) >= threshold;
    }

    /*!
     * \brief Overload for segment intersections.
     */
    template<class ctype, int wd, class Geo, class ctype2>
    bool isAdmissibleDistanceToBoundary(const Segment<ctype, wd>& seg,
                                        const Geo& entity,
                                        ctype2 threshold)
    {
        if (!isAdmissibleDistanceToBoundary(seg.source(), entity, threshold))
            return false;
        return isAdmissibleDistanceToBoundary(seg.target(), entity, threshold);
    }

    /*!
     * \brief Overload for ellipse arc intersections.
     */
    template<class ctype, int wd, class Geo, class ctype2>
    bool isAdmissibleDistanceToBoundary(const EllipseArc<ctype, wd>& arc,
                                        const Geo& entity,
                                        ctype2 threshold)
    {
        const auto sourceOnBound = pointOnGeometryBoundary(arc.source(), entity);
        const auto targetOnBound = pointOnGeometryBoundary(arc.target(), entity);

        // if any of the two rouches the boundary, only check source or target
        if (sourceOnBound || targetOnBound)
        {
            if (!sourceOnBound
                && computeDistanceToBoundary(arc.source(), entity) < threshold)
                return false;
            if (!targetOnBound
                && computeDistanceToBoundary(arc.target(), entity) < threshold)
                return false;
            return true;
        }

        // otherwise check a number of points along the arc
        std::vector<ctype> params({0.0, 0.25, 0.5, 0.75, 1.0});
        return std::all_of(params.begin(),
                           params.end(),
                           [&] (auto param) { return computeDistanceToBoundary(arc.getPoint(param), entity) >= threshold; });
    }

    /*!
     * \brief Overload for ellipse intersections on cylinder surfaces.
     */
    template<class ctype, int wd, class ctype2, class ctype3>
    bool isAdmissibleDistanceToBoundary(const Ellipse<ctype, wd>& arc,
                                        const CylinderSurface<ctype2>& entity,
                                        ctype3 threshold)
    {
        // for this we would first have to intersect the ellipse with the
        // upper and lower bounding circles to make sure there is no intersection
        // point!
        throw std::runtime_error(std::string("TODO: IMPLEMENT"));
    }

    /*!
     * \brief Overload for edge intersections on disks.
     */
    template<class ctype, class ctype2>
    bool isAdmissibleDistanceToBoundary(const TopoDS_Edge& edge,
                                        const Disk<ctype>& entity,
                                        ctype2 threshold)
    {
        const auto p1 = OCCUtilities::point( TopExp::FirstVertex(edge) );
        const auto p2 = OCCUtilities::point( TopExp::LastVertex(edge) );
        if (!isAdmissibleDistanceToBoundary(p1, entity, threshold))
            return false;
        return isAdmissibleDistanceToBoundary(p2, entity, threshold);
    }

    /*!
     * \brief Overload for face intersections on disks.
     */
    template<class ctype, class ctype2>
    bool isAdmissibleDistanceToBoundary(const TopoDS_Face& face,
                                        const Disk<ctype>& entity,
                                        ctype2 threshold)
    { throw std::runtime_error(std::string("NotImplemented: face-entity boundary distance")); }

    /*!
     * \brief Overload for edge intersections on faces.
     */
    template<class ctype>
    bool isAdmissibleDistanceToBoundary(const TopoDS_Edge& edge,
                                        const TopoDS_Face& entity,
                                        ctype threshold)
    {
        const auto p1 = OCCUtilities::point( TopExp::FirstVertex(edge) );
        const auto p2 = OCCUtilities::point( TopExp::LastVertex(edge) );
        if (!isAdmissibleDistanceToBoundary(p1, entity, threshold))
            return false;
        return isAdmissibleDistanceToBoundary(p2, entity, threshold);
    }

    /*!
     * \brief Overload for face intersections on faces.
     */
    template<class ctype>
    bool isAdmissibleDistanceToBoundary(const TopoDS_Face& face,
                                        const TopoDS_Face& entity,
                                        ctype threshold)
    { throw std::runtime_error(std::string("NotImplemented: face-entity boundary distance")); }

    /*!
     * \brief Overload for std::variant.
     */
    template<class... T, class Geo, class ctype>
    bool isAdmissibleDistanceToBoundary(const std::variant<T...>& intersection,
                                        const Geo& entity,
                                        ctype threshold)
    {
        return std::visit([&] (auto&& is)
                          { return isAdmissibleDistanceToBoundary(is, entity, threshold); },
                          intersection);
    }

    /*!
     * \brief Overload for a vector of intersections.
     */
    template<class T, class Geo, class ctype>
    bool isAdmissibleDistanceToBoundary(const std::vector<T>& intersection,
                                        const Geo& entity,
                                        ctype threshold)
    {
        return std::all_of(intersection.begin(),
                           intersection.end(),
                           [&] (const auto& is)
                           { return isAdmissibleDistanceToBoundary(is, entity, threshold); });
    }

} // end namespace ConstraintImpl
} // end namespace Frackit

#endif // FRACKIT_ENTITYNETWORK_CONSTRAINT_IMPL_ADMISSIBLE_BOUNDARY_DISTANCE_HH
