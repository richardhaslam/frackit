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
 * \brief \todo TODO doc me.
 */
#ifndef FRACKIT_GEOM_UTILITIES_HH
#define FRACKIT_GEOM_UTILITIES_HH

// Handle class used by OpenCascade
#include <Standard_Handle.hxx>

// Real number type used by OpenCascade
#include <Standard_Real.hxx>

// objects from geometric processors package
#include <gp_Ax2.hxx>
#include <gp_Elips.hxx>

// objects from the geometry package
#include <Geom_Ellipse.hxx>
#include <Geom_Curve.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <Geom_Surface.hxx>

// Shapes from the TopoDS package
#include <TopoDS_Edge.hxx>
#include <BRep_Tool.hxx>

// internal geometry classes
#include <frackit/geometry/ellipse.hh>
#include <frackit/geometry/ellipsearc.hh>

#include "gputilities.hh"

namespace Frackit {
namespace OCCUtilities {

    //! returns a Geom_Curve handle for an ellipse
    template<class ctype>
    Handle(Geom_Curve) getGeomHandle(const Ellipse<ctype, 3>& ellipse)
    {
        gp_Dir normal = direction(ellipse.normal());
        gp_Dir majorAx = direction(ellipse.majorAxis());
        gp_Pnt center = point(ellipse.center());
        gp_Elips gpEllipse(gp_Ax2(center, normal, majorAx),
                           ellipse.majorAxisLength(),
                           ellipse.minorAxisLength());
        return new Geom_Ellipse(gpEllipse);
    }

    //! returns a Geom_Curve handle for an ellipse arc
    template<class ctype>
    Handle(Geom_Curve) getGeomHandle(const EllipseArc<ctype, 3>& ellipseArc)
    {
        const auto& ellipse = ellipseArc.supportingEllipse();
        const auto geomEllipseHandle = getGeomHandle(ellipse);
        const auto angleSource = ellipseArc.sourceAngleOnEllipse();
        auto angleTarget = ellipseArc.targetAngleOnEllipse();
        if (angleTarget < angleSource)
            angleTarget += 2.0*M_PI;

        return new Geom_TrimmedCurve(geomEllipseHandle, angleSource, angleTarget);
    }

    //! returns a Geom_Curve handle for an edge
    Handle(Geom_Curve) getGeomHandle(const TopoDS_Edge& edge)
    {
        // get unbounded curve and parameter bounds (umin, umax)
        Standard_Real uMin, uMax;
        Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, uMin, uMax);

        // return trimmed curve
        return new Geom_TrimmedCurve(curve, uMin, uMax);
    }

    //! returns a Geom_Surface handle for a face
    Handle(Geom_Surface) getGeomHandle(const TopoDS_Face& face)
    { return BRep_Tool::Surface(face); }

} // end namespace OCCUtilities
} // end namespace Frackit

#endif // FRACKIT_GEOM_UTILITIES_HH
