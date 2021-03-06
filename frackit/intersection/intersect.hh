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
 * \brief Contains functionality to compute the
 *        intersections between geometries.
 */
#ifndef FRACKIT_INTERSECT_HH
#define FRACKIT_INTERSECT_HH

#include <cmath>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shell.hxx>

// Geometries of intersections
#include <frackit/geometry/line.hh>
#include <frackit/geometry/segment.hh>
#include <frackit/geometry/plane.hh>
#include <frackit/geometry/disk.hh>
#include <frackit/geometry/quadrilateral.hh>
#include <frackit/geometry/polygon.hh>
#include <frackit/geometry/cylindersurface.hh>
#include <frackit/geometry/ctype.hh>

#include <frackit/geometryutilities/name.hh>
#include <frackit/geometryutilities/getboundingbox.hh>

#include <frackit/occ/breputilities.hh>
#include <frackit/precision/defaultepsilon.hh>

#include "intersectiontraits.hh"
#include "bboxintersection.hh"
#include "emptyintersection.hh"

#include "algorithms/algo_segment_segment.hh"
#include "algorithms/algo_plane_plane.hh"
#include "algorithms/algo_plane_line.hh"
#include "algorithms/algo_disk_line.hh"
#include "algorithms/algo_quadrilateral_line.hh"
#include "algorithms/algo_quadrilateral_quadrilateral.hh"
#include "algorithms/algo_quadrilateral_disk.hh"
#include "algorithms/algo_polygon_disk.hh"
#include "algorithms/algo_polygon_polygon.hh"
#include "algorithms/algo_disk_disk.hh"
#include "algorithms/algo_cylsurface_disk.hh"
#include "algorithms/algo_cylsurface_face.hh"
#include "algorithms/algo_cylsurface_quadrilateral.hh"
#include "algorithms/algo_cylsurface_polygon.hh"
#include "algorithms/algo_face_face_3d.hh"

namespace Frackit {

// Forward declaration of the interface without given tolerance
template<class Geom1, class Geom2>
auto intersect(const Geom1& geo1, const Geom2& geo2);

/*!
 * \ingroup Intersection
 * \brief Interface for intersecting two geometries.
 * \param geo1 The first geometry
 * \param geo2 The second geometry
 * \param eps Tolerance to be used for floating point comparison
 */
template<class Geom1, class Geom2>
EmptyIntersection<0> intersect(const Geom1& geo1,
                               const Geom2& geo2,
                               typename CoordinateTypeTraits<Geom1>::type eps)
{
    std::string msg = "Intersection algorithm between \"";
    msg += geometryName(geo1);
    msg += "\" and \"";
    msg += geometryName(geo2);
    msg += "\" not implemented";
    throw std::runtime_error(msg);
}

/*!
 * \ingroup Intersection
 * \brief Intersect two segments.
 * \param segment1 The first segment
 * \param segment2 The second segment
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype, int wd>
Intersection< Segment<ctype, wd>, Segment<ctype, wd> >
intersect(const Segment<ctype, wd>& segment1, const Segment<ctype, wd>& segment2, ctype eps)
{
    if (!doIntersect(getBoundingBox(segment1), getBoundingBox(segment2), eps))
        return {EmptyIntersection<wd, ctype>()};
    return IntersectionAlgorithms::intersect_segment_segment(segment1, segment2, eps);
}

/*!
 * \ingroup Intersection
 * \brief Intersect two planes.
 * \param plane1 The first plane
 * \param plane2 The second plane
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Plane<ctype, 3>, Plane<ctype, 3> >
intersect(const Plane<ctype, 3>& plane1, const Plane<ctype, 3>& plane2, ctype eps)
{ return IntersectionAlgorithms::intersect_plane_plane(plane1, plane2, eps); }

/*!
 * \ingroup Intersection
 * \brief Intersect a plane and a line.
 * \param plane The plane
 * \param line The line
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Plane<ctype, 3>, Line<ctype, 3> >
intersect(const Plane<ctype, 3>& plane, const Line<ctype, 3>& line, ctype eps)
{ return IntersectionAlgorithms::intersect_plane_line(plane, line, eps); }

/*!
 * \ingroup Intersection
 * \brief Intersect a line and a plane.
 * \param line The plane
 * \param plane The line
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Line<ctype, 3>, Plane<ctype, 3> >
intersect(const Line<ctype, 3>& line, const Plane<ctype, 3>& plane, ctype eps)
{ return intersect(plane, line, eps); }

/*!
 * \ingroup Intersection
 * \brief Intersect a disk and a line.
 * \param disk The disk
 * \param line The line
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Disk<ctype>, Line<ctype, 3> >
intersect(const Disk<ctype>& disk, const Line<ctype, 3>& line, ctype eps)
{ return IntersectionAlgorithms::intersect_disk_line(disk, line, eps); }

/*!
 * \ingroup Intersection
 * \brief Intersect a line and a disk.
 * \param line The line
 * \param disk The disk
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Line<ctype, 3>, Disk<ctype> >
intersect(const Line<ctype, 3>& line, const Disk<ctype>& disk, ctype eps)
{ return intersect(disk, line, eps); }

/*!
 * \ingroup Intersection
 * \brief Intersect a three-dimensional quadrilateral and a line.
 * \param quad The quadrilateral
 * \param line The line
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Quadrilateral<ctype, 3>, Line<ctype, 3> >
intersect(const Quadrilateral<ctype, 3>& quad, const Line<ctype, 3>& line, ctype eps)
{ return IntersectionAlgorithms::intersect_quadrilateral_line(quad, line, eps); }


/*!
 * \ingroup Intersection
 * \brief Intersect a line and a three-dimensional quadrilateral.
 * \param line The line
 * \param quad The quadrilateral
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Line<ctype, 3>, Quadrilateral<ctype, 3>>
intersect(const Line<ctype, 3>& line, const Quadrilateral<ctype, 3>& quad, ctype eps)
{ return intersect(quad, line, eps); }

/*!
 * \ingroup Intersection
 * \brief Intersect two disks.
 * \param disk1 The first disk
 * \param disk2 The second disk
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Disk<ctype>, Disk<ctype> >
intersect(const Disk<ctype>& disk1, const Disk<ctype>& disk2, ctype eps)
{
    if (!doIntersect(getBoundingBox(disk1), getBoundingBox(disk2), eps))
        return {EmptyIntersection<3, ctype>()};
    return IntersectionAlgorithms::intersect_disk_disk(disk1, disk2, eps);
}

/*!
 * \ingroup Intersection
 * \brief Intersect two quadrilaterals in 3d space.
 * \param quad1 The first quadrilateral
 * \param quad2 The second quadrilateral
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Quadrilateral<ctype, 3>, Quadrilateral<ctype, 3> >
intersect(const Quadrilateral<ctype, 3>& quad1, const Quadrilateral<ctype, 3>& quad2, ctype eps)
{
    if (!doIntersect(getBoundingBox(quad1), getBoundingBox(quad2), eps))
        return {EmptyIntersection<3, ctype>()};
    return IntersectionAlgorithms::intersect_quadrilateral_quadrilateral(quad1, quad2, eps);
}

/*!
 * \ingroup Intersection
 * \brief Intersect two polygons in 3d space.
 * \param polygon1 The first polygon
 * \param polygon2 The second polygon
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Polygon<ctype, 3>, Polygon<ctype, 3> >
intersect(const Polygon<ctype, 3>& polygon1, const Polygon<ctype, 3>& polygon2, ctype eps)
{
    if (!doIntersect(getBoundingBox(polygon1), getBoundingBox(polygon2), eps))
        return {EmptyIntersection<3, ctype>()};
    return IntersectionAlgorithms::intersect_polygon_polygon(polygon1, polygon2, eps);
}

/*!
 * \ingroup Intersection
 * \brief Intersect a quadrilateral and a disk in 3d space.
 * \param quad The quadrilateral
 * \param disk The disk
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Quadrilateral<ctype, 3>, Disk<ctype> >
intersect(const Quadrilateral<ctype, 3>& quad, const Disk<ctype>& disk, ctype eps)
{
    if (!doIntersect(getBoundingBox(quad), getBoundingBox(disk), eps))
        return {EmptyIntersection<3, ctype>()};
    return IntersectionAlgorithms::intersect_quadrilateral_disk(quad, disk, eps);
}

/*!
 * \ingroup Intersection
 * \brief Intersect a disk and a quadrilateral in 3d space.
 * \param disk The disk
 * \param quad The quadrilateral
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Disk<ctype>, Quadrilateral<ctype, 3> >
intersect(const Disk<ctype>& disk, const Quadrilateral<ctype, 3>& quad, ctype eps)
{ return intersect(quad, disk, eps); }

/*!
 * \ingroup Intersection
 * \brief Intersect a polygon and a disk in 3d space.
 * \param polygon The polygon
 * \param disk The disk
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Polygon<ctype, 3>, Disk<ctype> >
intersect(const Polygon<ctype, 3>& polygon, const Disk<ctype>& disk, ctype eps)
{
    if (!doIntersect(getBoundingBox(polygon), getBoundingBox(disk), eps))
        return {EmptyIntersection<3, ctype>()};
    return IntersectionAlgorithms::intersect_polygon_disk(polygon, disk, eps);
}

/*!
 * \ingroup Intersection
 * \brief Intersect a disk and a polygon in 3d space.
 * \param polygon The polygon
 * \param quad The quadrilateral
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Disk<ctype>, Polygon<ctype, 3> >
intersect(const Disk<ctype>& disk, const Polygon<ctype, 3>& polygon, ctype eps)
{ return intersect(polygon, disk, eps); }

/*!
 * \ingroup Intersection
 * \brief Intersect a lateral cylinder surface and a disk.
 * \param cylSurface The lateral cylinder surface
 * \param disk The disk
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< CylinderSurface<ctype>, Disk<ctype> >
intersect(const CylinderSurface<ctype>& cylSurface, const Disk<ctype>& disk, ctype eps)
{
    if (!doIntersect(getBoundingBox(cylSurface), getBoundingBox(disk), eps))
        return {EmptyIntersection<3, ctype>()};
    return IntersectionAlgorithms::intersect_cylinderSurface_disk(cylSurface, disk, eps);
}

/*!
 * \ingroup Intersection
 * \brief Intersect a disk and a lateral cylinder surface.
 * \param disk The disk
 * \param cylSurface The lateral cylinder surface
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Disk<ctype>, CylinderSurface<ctype> >
intersect(const Disk<ctype>& disk, const CylinderSurface<ctype>& cylSurface, ctype eps)
{ return intersect(cylSurface, disk, eps); }

/*!
 * \ingroup Intersection
 * \brief Intersect a lateral cylinder surface and a quadrilateral.
 * \param cylSurface The lateral cylinder surface
 * \param quad The quadrilateral
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< CylinderSurface<ctype>, Quadrilateral<ctype, 3> >
intersect(const CylinderSurface<ctype>& cylSurface, const Quadrilateral<ctype, 3>& quad, ctype eps)
{
    if (!doIntersect(getBoundingBox(cylSurface), getBoundingBox(quad), eps))
        return {EmptyIntersection<3, ctype>()};
    return IntersectionAlgorithms::intersect_cylinderSurface_quadrilateral(cylSurface, quad, eps);
}

/*!
 * \ingroup Intersection
 * \brief Intersect a quadrilateral and a lateral cylinder surface.
 * \param quad The quadrilateral
 * \param cylSurface The lateral cylinder surface
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Quadrilateral<ctype, 3>, CylinderSurface<ctype> >
intersect(const Quadrilateral<ctype, 3>& quad, const CylinderSurface<ctype>& cylSurface, ctype eps)
{ return intersect(cylSurface, quad, eps); }

/*!
 * \ingroup Intersection
 * \brief Intersect a lateral cylinder surface and a polygon.
 * \param cylSurface The lateral cylinder surface
 * \param polygon The polygon
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< CylinderSurface<ctype>, Polygon<ctype, 3> >
intersect(const CylinderSurface<ctype>& cylSurface, const Polygon<ctype, 3>& polygon, ctype eps)
{
    if (!doIntersect(getBoundingBox(cylSurface), getBoundingBox(polygon), eps))
        return {EmptyIntersection<3, ctype>()};
    return IntersectionAlgorithms::intersect_cylinderSurface_polygon(cylSurface, polygon, eps);
}

/*!
 * \ingroup Intersection
 * \brief Intersect a polygon and a lateral cylinder surface.
 * \param polygon The polygon
 * \param cylSurface The lateral cylinder surface
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Polygon<ctype, 3>, CylinderSurface<ctype> >
intersect(const Polygon<ctype, 3>& polygon, const CylinderSurface<ctype>& cylSurface, ctype eps)
{ return intersect(cylSurface, polygon, eps); }

/*!
 * \ingroup Intersection
 * \brief Intersect a lateral cylinder surface and a TopoDS_Face.
 * \param cylSurface The lateral cylinder surface
 * \param face The face shape
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< CylinderSurface<ctype>, TopoDS_Face >
intersect(const CylinderSurface<ctype>& cylSurface, const TopoDS_Face& face, ctype eps)
{
    if (!doIntersect(getBoundingBox(cylSurface), getBoundingBox(face), eps))
        return {EmptyIntersection<3, ctype>()};
    return IntersectionAlgorithms::intersect_cylinderSurface_face(cylSurface, face, eps);
}

/*!
 * \ingroup Intersection
 * \brief Intersect a TopoDS_Face and a lateral cylinder surface.
 * \param face The face shape
 * \param cylSurface The lateral cylinder surface
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< TopoDS_Face, CylinderSurface<ctype> >
intersect(const TopoDS_Face& face, const CylinderSurface<ctype>& cylSurface, ctype eps)
{ return intersect(cylSurface, face, eps); }

/*!
 * \ingroup Intersection
 * \brief Intersect two face shapes.
 * \param face1 The first face shape
 * \param face2 The second face shape
 * \param eps Tolerance to be used for floating point comparisons
 * \todo TODO: How to distinguish here if 2d setups are considered?
 */
template<class ctype>
Intersection< TopoDS_Face, TopoDS_Face >
intersect(const TopoDS_Face& face1, const TopoDS_Face& face2, ctype eps)
{
    if (!doIntersect(getBoundingBox(face1), getBoundingBox(face2), eps))
        return {EmptyIntersection<3, ctype>()};
    return IntersectionAlgorithms::intersect_face_face_3d(face1, face2, eps);
}

/*!
 * \ingroup Intersection
 * \brief Intersect a disk and a face shape.
 * \param disk The disk
 * \param face The face shape
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Disk<ctype>, TopoDS_Face >
intersect(const Disk<ctype>& disk, const TopoDS_Face& face, ctype eps)
{ return intersect(OCCUtilities::getShape(disk), face, eps); }

/*!
 * \ingroup Intersection
 * \brief Intersect a face shape and a disk.
 * \param face The face shape
 * \param disk The disk
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< TopoDS_Face, Disk<ctype> >
intersect(const TopoDS_Face& face, const Disk<ctype>& disk, ctype eps)
{ return intersect(disk, face, eps); }

/*!
 * \ingroup Intersection
 * \brief Intersect a quadrilateral and a face shape.
 * \param quad The quadrilateral
 * \param face The face shape
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Quadrilateral<ctype, 3>, TopoDS_Face >
intersect(const Quadrilateral<ctype, 3>& quad, const TopoDS_Face& face, ctype eps)
{ return intersect(OCCUtilities::getShape(quad), face, eps); }

/*!
 * \ingroup Intersection
 * \brief Intersect a face shape and a quadrilateral.
 * \param face The face shape
 * \param quad The quadrilateral
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< TopoDS_Face, Quadrilateral<ctype, 3> >
intersect(const TopoDS_Face& face, const Quadrilateral<ctype, 3>& quad, ctype eps)
{ return intersect(quad, face, eps); }

/*!
 * \ingroup Intersection
 * \brief Intersect a polygon and a face shape.
 * \param polygon The polygon
 * \param face The face shape
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< Polygon<ctype, 3>, TopoDS_Face >
intersect(const Polygon<ctype, 3>& polygon, const TopoDS_Face& face, ctype eps)
{ return intersect(OCCUtilities::getShape(polygon), face, eps); }

/*!
 * \ingroup Intersection
 * \brief Intersect a face shape and a polygon.
 * \param face The face shape
 * \param polygon The polygon
 * \param eps Tolerance to be used for floating point comparisons
 */
template<class ctype>
Intersection< TopoDS_Face, Polygon<ctype, 3> >
intersect(const TopoDS_Face& face, const Polygon<ctype, 3>& polygon, ctype eps)
{ return intersect(polygon, face, eps); }

/*!
 * \ingroup Intersection
 * \brief Interface for intersecting two geometries.
 * \param geo1 The first geometry
 * \param geo2 The second geometry
 * \note This overload selects a default tolerance
 */
template<class Geom1, class Geom2>
auto intersect(const Geom1& geo1, const Geom2& geo2)
{
    using std::min;
    const auto eps1 = defaultEpsilon(geo1);
    const auto eps2 = defaultEpsilon(geo2);
    return intersect(geo1, geo2, min(eps1, eps2));
}

} // end namespace Frackit

#endif // FRACKIT_INTERSECT_HH
