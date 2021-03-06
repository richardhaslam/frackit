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
 * \ingroup Magnitude
 * \brief Contains functionality for computing the
 *        lengths of one-dimensional geometries.
 */
#ifndef FRACKIT_MAGNITUDE_LENGTH_HH
#define FRACKIT_MAGNITUDE_LENGTH_HH

#include <cmath>
#include <numeric>

#include <Standard_Handle.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Wire.hxx>

#include <Geom_Curve.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <GCPnts_AbscissaPoint.hxx>

#include <frackit/occ/geomutilities.hh>
#include <frackit/occ/breputilities.hh>

#include <frackit/geometry/circle.hh>
#include <frackit/geometry/ellipse.hh>
#include <frackit/geometry/ellipsearc.hh>

namespace Frackit {

/*!
 * \ingroup Magnitude
 * \brief Returns the length of a geometry.
 * \note This is the default overload trying to get
 *       the length from the geometry itself.
 */
template<class Geometry>
typename Geometry::ctype computeLength(const Geometry& geom)
{ return geom.length(); }

/*!
 * \ingroup Magnitude
 * \brief Returns the length of a curve of the Geom package.
 */
template<class ctype = double>
ctype computeLength(const Handle(Geom_Curve)& curve)
{
    const auto uMin = curve->FirstParameter();
    const auto uMax = curve->LastParameter();
    GeomAdaptor_Curve adaptorCurve(curve, uMin, uMax);
    return GCPnts_AbscissaPoint::Length(adaptorCurve, uMin, uMax);
}

/*!
 * \ingroup Magnitude
 * \brief Returns the length of a BRep edge.
 */
template<class ctype = double>
ctype computeLength(const TopoDS_Edge& edge)
{ return computeLength(OCCUtilities::getGeomHandle(edge)); }

/*!
 * \ingroup Magnitude
 * \brief Returns the length of a BRep wire.
 */
template<class ctype = double>
ctype computeLength(const TopoDS_Wire& wire)
{
    const auto edges = OCCUtilities::getEdges(wire);
    auto addLength = [] (const ctype& cur, const TopoDS_Edge& edge)
                        { return cur + computeLength(edge); };
    return std::accumulate(edges.begin(), edges.end(), ctype(0.0), addLength);
}

/*!
 * \ingroup Magnitude
 * \brief Returns the length of a circle (circumference).
 */
template<class ctype, int worldDim>
ctype computeLength(const Circle<ctype, worldDim>& circle)
{ return 2.0*M_PI*circle.radius(); }

/*!
 * \ingroup Magnitude
 * \brief Returns the length of an ellipse.
 */
template<class ctype, int worldDim>
ctype computeLength(const Ellipse<ctype, worldDim>& ellipse)
{ return M_PI*(ellipse.majorAxisLength() + ellipse.minorAxisLength()); }

/*!
 * \ingroup Magnitude
 * \brief Returns the length of an elliptical arc.
 */
template<class ctype>
ctype computeLength(const EllipseArc<ctype, 3>& arc)
{
    if (arc.isFullEllipse())
         return M_PI*(arc.majorAxisLength() + arc.minorAxisLength());
    return computeLength(OCCUtilities::getGeomHandle(arc));
}

} // end namespace Frackit

#endif // FRACKIT_MAGNITUDE_LENGTH_HH
