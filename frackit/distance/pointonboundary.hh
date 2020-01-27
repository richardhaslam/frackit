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
 * \ingroup Distance
 * \brief Contains functionality to evaluate if points are contained
 *        on the boundary of geometries, i.e. if the distance to the
 *        boundary is below a given numerical threshold.
 */
#ifndef FRACKIT_POINT_ON_GEOMETRY_BOUNDARY_HH
#define FRACKIT_POINT_ON_GEOMETRY_BOUNDARY_HH

#include <string>

#include <BRepTools.hxx>
#include <TopoDS_Face.hxx>

#include <frackit/common/extractdimension.hh>
#include <frackit/precision/defaultepsilon.hh>

#include <frackit/geometry/point.hh>
#include <frackit/geometry/disk.hh>
#include <frackit/geometry/cylindersurface.hh>

#include <frackit/geometryutilities/name.hh>

#include "pointongeometry.hh"

namespace Frackit {

/*!
 * \ingroup Distance
 * \brief Evaluate if a point lies on the boundary of a geometry.
 * \param p The point
 * \param geo The geometry
 * \param eps Epsilon value to be used for the check
 * \note This is the default overload throwing an error. Overloads
 *       have to be implemented for the different geometries.
 */
template<class ctype1, int wd, class Geo, class ctype2>
bool pointOnGeometryBoundary(const Point<ctype1, wd>& p, const Geo& geo, ctype2 eps)
{
    static_assert(wd == DimensionalityTraits<Geo>::worldDimension(), "World dimension mismatch");
    std::string msg = "Point on boundary not implemented for ";
    msg += "\"" + geometryName(geo) + "\" \n";
    throw std::runtime_error(msg);
}

/*!
 * \ingroup Distance
 * \brief Evaluate if a point lies on the boundary of a disk.
 * \param p The point
 * \param disk The disk
 * \param eps Epsilon value to be used for the check
 */
template<class ctype1, int wd, class ctype2>
bool pointOnGeometryBoundary(const Point<ctype1, wd>& p,
                             const Disk<ctype2>& disk,
                             ctype2 eps)
{ return pointOnGeometry(p, disk.boundingEllipse(), eps); }

/*!
 * \ingroup Distance
 * \brief Evaluate if a point lies on the boundary of
 *        a quadrilateral in 3d space.
 * \param p The point
 * \param quad The quadrilateral
 * \param eps Epsilon value to be used for the check
 */
template<class ctype1, class ctype2, class ctype3>
bool pointOnGeometryBoundary(const Point<ctype1, 3>& p,
                             const Quadrilateral<ctype2, 3>& quad,
                             ctype3 eps)
{
    for (unsigned int i = 0; i < quad.numEdges(); ++i)
        if (pointOnGeometry(p, quad.edge(i), eps))
            return true;
    return false;
}

/*!
 * \ingroup Distance
 * \brief Evaluate if a point lies on the boundary of a cylinder surface.
 * \param p The point
 * \param cylSurface The cylinder surface
 * \param eps Epsilon value to be used for the check
 */
template<class ctype1, int wd, class ctype2>
bool pointOnGeometryBoundary(const Point<ctype1, wd>& p,
                             const CylinderSurface<ctype2>& cylSurface,
                             ctype2 eps)
{
    if (pointOnGeometry(p, cylSurface.upperBoundingCircle(), eps))
        return true;
    return pointOnGeometry(p, cylSurface.lowerBoundingCircle(), eps);
}

/*!
 * \ingroup Distance
 * \brief Evaluate if a point lies on the boundary of a face shape.
 * \param p The point
 * \param face The face shape
 * \param eps Epsilon value to be used for the check
 */
template<class ctype1, class ctype2>
bool pointOnGeometryBoundary(const Point<ctype1, 3>& p,
                             const TopoDS_Face& face,
                             ctype2 eps)
{ return pointOnGeometry(p, BRepTools::OuterWire(face)); }

/*!
 * \ingroup Distance
 * \brief Evaluate if a point lies on the boundary of a geometry.
 * \param p The point
 * \param geo The geometry
 * \note This overload selects the default epsilon
 */
template<class ctype1, int wd, class Geo>
bool pointOnGeometryBoundary(const Point<ctype1, wd>& p, const Geo& geo)
{ return pointOnGeometryBoundary(p, geo, defaultEpsilon(geo)); }

} // end namespace Frackit

#endif // FRACKIT_POINT_ON_GEOMETRY_BOUNDARY_HH
