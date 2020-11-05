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
 * \ingroup Distance
 * \brief Contains functionality for computing the distance
 *        of a geometry to the boundary of another geometry.
 *        The difference to the distance computation functions
 *        is that those return a distance of zero if the geometries
 *        intersect, while here, it is only zero if it intersects the boundary.
 *        For instance, the distance of a point to a box is zero if the point
 *        is inside the box, but the distance to the box surface may be non-zero.
 */
#ifndef FRACKIT_DISTANCE_TO_BOUNDARY_HH
#define FRACKIT_DISTANCE_TO_BOUNDARY_HH

#include <cmath>
#include <limits>
#include <algorithm>

#include <Extrema_ExtAlgo.hxx>
#include <Extrema_ExtFlag.hxx>
#include <BRepTools.hxx>
#include <BRepClass3d.hxx>

#include <frackit/precision/precision.hh>
#include <frackit/geometry/disk.hh>
#include <frackit/geometry/quadrilateral.hh>
#include <frackit/geometry/polygon.hh>
#include <frackit/geometry/cylindersurface.hh>

#include <frackit/geometryutilities/name.hh>

#include "distance.hh"

namespace Frackit {

/*!
 * \ingroup Distance
 * \brief Compute the distance of a geometry
 *        to the boundary of another geometry.
 * \param geo1 The first geometry
 * \param geo2 The second geometry
 * \param deflection The epsilon used in the BrepExtrema command
 * \param extFlag The flag passed to the BrepExtrema command (MIN/MAX/MINMAX)
 * \param extAlgo The algorithm passed to the BrepExtrema command (TREE/GRAD)
 * \note This is the default overload throwing an error. Overloads have to be
 *       implemented for pairs of geometries.
 */
template<class Geo1, class Geo2>
Impl::PCT<Geo1, Geo2>
computeDistanceToBoundary(const Geo1& geo1,
                          const Geo2& geo2,
                          Impl::PCT<Geo1, Geo2> deflection = Precision<Impl::PCT<Geo1, Geo2>>::confusion(),
                          Extrema_ExtFlag extFlag = Extrema_ExtFlag_MINMAX,
                          Extrema_ExtAlgo extAlgo = Extrema_ExtAlgo_Grad)
{
    std::string msg = "Distance to boundary computation not implemented for ";
    msg += "\"" + geometryName(geo1) + "\"";
    msg += " and ";
    msg += "\"" + geometryName(geo2)  + "\"";
    throw std::runtime_error(msg);
}

/*!
 * \ingroup Distance
 * \brief Compute the distance of a geometry
 *        to the bounding ellipse of a disk.
 * \param geo The geometry
 * \param disk The disk
 * \param deflection The epsilon used in the BrepExtrema command
 * \param extFlag The flag passed to the BrepExtrema command (MIN/MAX/MINMAX)
 * \param extAlgo The algorithm passed to the BrepExtrema command (TREE/GRAD)
 */
template<class Geo, class ctype>
Impl::PCT<Geo, Disk<ctype>>
computeDistanceToBoundary(const Geo& geo,
                          const Disk<ctype>& disk,
                          Impl::PCT<Geo, Disk<ctype>> deflection
                            = Precision<Impl::PCT<Geo, Disk<ctype>>>::confusion(),
                          Extrema_ExtFlag extFlag = Extrema_ExtFlag_MINMAX,
                          Extrema_ExtAlgo extAlgo = Extrema_ExtAlgo_Grad)
{
    return computeDistance(OCCUtilities::getShape(geo),
                           OCCUtilities::getShape(disk.boundingEllipse()),
                           deflection,
                           extFlag,
                           extAlgo);
}

/*!
 * \ingroup Distance
 * \brief Compute the distance of a point
 *        to the bounding of a quadrilateral.
 * \param p The point
 * \param quad The quadrilateral
 */
template<class ctype1, class ctype2>
PromotedType<ctype1, ctype2>
computeDistanceToBoundary(const Point<ctype1, 3>& p,
                          const Quadrilateral<ctype2, 3>& quad)
{
    using std::min;
    using ctype = PromotedType<ctype1, ctype2>;

    ctype distance = std::numeric_limits<ctype>::max();
    for (unsigned int i = 0; i < quad.numEdges(); ++i)
        distance = min(distance, computeDistance(quad.edge(i), p));
    return distance;
}

/*!
 * \ingroup Distance
 * \brief Compute the distance of a point
 *        to the bounding of a polygon.
 * \param p The point
 * \param polygon The polygon
 */
template<class ctype1, class ctype2>
PromotedType<ctype1, ctype2>
computeDistanceToBoundary(const Point<ctype1, 3>& p,
                          const Polygon<ctype2, 3>& polygon)
{
    using std::min;
    using ctype = PromotedType<ctype1, ctype2>;

    ctype distance = std::numeric_limits<ctype>::max();
    for (unsigned int i = 0; i < polygon.numEdges(); ++i)
        distance = min(distance, computeDistance(polygon.edge(i), p));
    return distance;
}

/*!
 * \ingroup Distance
 * \brief Compute the distance of a geometry
 *        to the bounding circles of a cylinder surface.
 * \param geo The geometry
 * \param cylSurface The cylinder surface
 * \param deflection The epsilon used in the BrepExtrema command
 * \param extFlag The flag passed to the BrepExtrema command (MIN/MAX/MINMAX)
 * \param extAlgo The algorithm passed to the BrepExtrema command (TREE/GRAD)
 */
template<class Geo, class ctype>
Impl::PCT<Geo, CylinderSurface<ctype>>
computeDistanceToBoundary(const Geo& geo,
                          const CylinderSurface<ctype>& cylSurface,
                          Impl::PCT<Geo, CylinderSurface<ctype>> deflection
                            = Precision<Impl::PCT<Geo, CylinderSurface<ctype>>>::confusion(),
                          Extrema_ExtFlag extFlag = Extrema_ExtFlag_MINMAX,
                          Extrema_ExtAlgo extAlgo = Extrema_ExtAlgo_Grad)
{
    using std::min;
    return min(computeDistance(OCCUtilities::getShape(geo),
                               OCCUtilities::getShape(cylSurface.lowerBoundingCircle()),
                               deflection,
                               extFlag,
                               extAlgo),
               computeDistance(OCCUtilities::getShape(geo),
                               OCCUtilities::getShape(cylSurface.upperBoundingCircle()),
                               deflection,
                               extFlag,
                               extAlgo));
}

/*!
 * \ingroup Distance
 * \brief Compute the distance of a shape
 *        to the bounding ellipse of a disk.
 * \param shape The shape
 * \param disk The disk
 * \param deflection The epsilon used in the BrepExtrema command
 * \param extFlag The flag passed to the BrepExtrema command (MIN/MAX/MINMAX)
 * \param extAlgo The algorithm passed to the BrepExtrema command (TREE/GRAD)
 */
template<class ctype>
ctype computeDistanceToBoundary(const TopoDS_Shape& shape,
                                const Disk<ctype>& disk,
                                ctype deflection = Precision<ctype>::confusion(),
                                Extrema_ExtFlag extFlag = Extrema_ExtFlag_MINMAX,
                                Extrema_ExtAlgo extAlgo = Extrema_ExtAlgo_Grad)
{
    return computeDistance(shape,
                           OCCUtilities::getShape(disk.boundingEllipse()),
                           deflection,
                           extFlag,
                           extAlgo);
}

/*!
 * \ingroup Distance
 * \brief Compute the distance of a shape
 *        to the boundary of a TopoDS_Face.
 * \param shape The shape
 * \param face The TopoDS_Face
 * \param deflection The epsilon used in the BrepExtrema command
 * \param extFlag The flag passed to the BrepExtrema command (MIN/MAX/MINMAX)
 * \param extAlgo The algorithm passed to the BrepExtrema command (TREE/GRAD)
 * \note This only computes the distance to the outer wire of the face.
 *       Thus, this assumes faces that do not have a whole inside.
 */
template<class ctype = double>
ctype computeDistanceToBoundary(const TopoDS_Shape& shape,
                                const TopoDS_Face& face,
                                ctype deflection = Precision<ctype>::confusion(),
                                Extrema_ExtFlag extFlag = Extrema_ExtFlag_MINMAX,
                                Extrema_ExtAlgo extAlgo = Extrema_ExtAlgo_Grad)
{
    return computeDistance(shape,
                           BRepTools::OuterWire(face),
                           deflection,
                           extFlag,
                           extAlgo);
}

/*!
 * \ingroup Distance
 * \brief Compute the distance of a shape
 *        to the bounding wire of a quadrilateral.
 * \param shape The shape
 * \param quad The quadrilateral
 * \param deflection The epsilon used in the BrepExtrema command
 * \param extFlag The flag passed to the BrepExtrema command (MIN/MAX/MINMAX)
 * \param extAlgo The algorithm passed to the BrepExtrema command (TREE/GRAD)
 */
template<class ctype>
ctype computeDistanceToBoundary(const TopoDS_Shape& shape,
                                const Quadrilateral<ctype, 3>& quad,
                                ctype deflection = Precision<ctype>::confusion(),
                                Extrema_ExtFlag extFlag = Extrema_ExtFlag_MINMAX,
                                Extrema_ExtAlgo extAlgo = Extrema_ExtAlgo_Grad)
{
    using std::min;
    ctype minDist = std::numeric_limits<ctype>::max();
    for (unsigned int i = 0; i < quad.numCorners(); ++i)
        minDist = min(minDist, computeDistance(shape,
                                               quad.edge(i),
                                               deflection, extFlag, extAlgo));
    return minDist;
}

/*!
 * \ingroup Distance
 * \brief Compute the distance of a shape
 *        to the bounding wire of a polygon.
 * \param shape The shape
 * \param polygon The polygon
 * \param deflection The epsilon used in the BrepExtrema command
 * \param extFlag The flag passed to the BrepExtrema command (MIN/MAX/MINMAX)
 * \param extAlgo The algorithm passed to the BrepExtrema command (TREE/GRAD)
 */
template<class ctype>
ctype computeDistanceToBoundary(const TopoDS_Shape& shape,
                                const Polygon<ctype, 3>& polygon,
                                ctype deflection = Precision<ctype>::confusion(),
                                Extrema_ExtFlag extFlag = Extrema_ExtFlag_MINMAX,
                                Extrema_ExtAlgo extAlgo = Extrema_ExtAlgo_Grad)
{
    using std::min;
    ctype minDist = std::numeric_limits<ctype>::max();
    for (unsigned int i = 0; i < polygon.numCorners(); ++i)
        minDist = min(minDist, computeDistance(shape,
                                               polygon.edge(i),
                                               deflection, extFlag, extAlgo));
    return minDist;
}

/*!
 * \ingroup Distance
 * \brief Compute the distance of a geometry
 *        to the boundary of a TopoDS_Face.
 * \param geo The geometry
 * \param face The TopoDS_Face
 * \param deflection The epsilon used in the BrepExtrema command
 * \param extFlag The flag passed to the BrepExtrema command (MIN/MAX/MINMAX)
 * \param extAlgo The algorithm passed to the BrepExtrema command (TREE/GRAD)
 */
template<class Geo>
typename Geo::ctype
computeDistanceToBoundary(const Geo& geo,
                          const TopoDS_Face& face,
                          typename Geo::ctype deflection
                          = Precision<typename Geo::ctype>::confusion(),
                          Extrema_ExtFlag extFlag = Extrema_ExtFlag_MINMAX,
                          Extrema_ExtAlgo extAlgo = Extrema_ExtAlgo_Grad)
{
    return computeDistance(OCCUtilities::getShape(geo),
                           BRepTools::OuterWire(face),
                           deflection,
                           extFlag,
                           extAlgo);
}

/*!
 * \ingroup Distance
 * \brief Compute the distance of a geometry
 *        to the boundary of a TopoDS_Solid.
 * \param geo The geometry
 * \param solid The TopoDS_Solid
 * \param deflection The epsilon used in the BrepExtrema command
 * \param extFlag The flag passed to the BrepExtrema command (MIN/MAX/MINMAX)
 * \param extAlgo The algorithm passed to the BrepExtrema command (TREE/GRAD)
 */
template<class Geo>
typename Geo::ctype
computeDistanceToBoundary(const Geo& geo,
                          const TopoDS_Solid& solid,
                          typename Geo::ctype deflection
                          = Precision<typename Geo::ctype>::confusion(),
                          Extrema_ExtFlag extFlag = Extrema_ExtFlag_MINMAX,
                          Extrema_ExtAlgo extAlgo = Extrema_ExtAlgo_Grad)
{
    return computeDistance(OCCUtilities::getShape(geo),
                           BRepClass3d::OuterShell(solid),
                           deflection,
                           extFlag,
                           extAlgo);
}

} // end namespace Frackit

#endif // FRACKIT_DISTANCE_TO_BOUNDARY_HH
