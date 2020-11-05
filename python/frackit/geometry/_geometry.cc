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
#include <pybind11/pybind11.h>

#include <frackit/python/geometry/geometry.hh>
#include <frackit/python/geometry/point.hh>
#include <frackit/python/geometry/direction.hh>
#include <frackit/python/geometry/plane.hh>
#include <frackit/python/geometry/vector.hh>
#include <frackit/python/geometry/segment.hh>
#include <frackit/python/geometry/line.hh>
#include <frackit/python/geometry/triangle.hh>
#include <frackit/python/geometry/quadrilateral.hh>
#include <frackit/python/geometry/polygon.hh>
#include <frackit/python/geometry/ellipticalgeometry.hh>
#include <frackit/python/geometry/ellipse.hh>
#include <frackit/python/geometry/ellipsearc.hh>
#include <frackit/python/geometry/disk.hh>
#include <frackit/python/geometry/circle.hh>
#include <frackit/python/geometry/cylinder.hh>
#include <frackit/python/geometry/hollowcylinder.hh>
#include <frackit/python/geometry/cylindersurface.hh>
#include <frackit/python/geometry/box.hh>
#include <frackit/python/geometry/sphere.hh>
#include <frackit/python/geometry/brepwrappers.hh>
#include <frackit/python/geometry/emptyintersection.hh>

#include <frackit/python/geometry/distance.hh>
#include <frackit/python/geometry/distancetoboundary.hh>
#include <frackit/python/geometry/intersect.hh>

PYBIND11_MODULE(_geometry, module)
{
    Frackit::Python::registerGeometry(module);

    Frackit::Python::registerPoint<double>(module);

    Frackit::Python::registerVector<double>(module);
    Frackit::Python::registerDirection<double>(module);
    Frackit::Python::registerSegment<double>(module);
    Frackit::Python::registerLine<double>(module);

    Frackit::Python::registerEllipticalGeometry<double>(module);
    Frackit::Python::registerEllipse<double>(module);
    Frackit::Python::registerEllipseArc<double>(module);
    Frackit::Python::registerCircle<double>(module);

    Frackit::Python::registerPlane<double>(module);
    Frackit::Python::registerTriangle<double>(module);
    Frackit::Python::registerQuadrilateral<double>(module);
    Frackit::Python::registerPolygon<double>(module);
    Frackit::Python::registerDisk<double>(module);
    Frackit::Python::registerCylinderSurface<double>(module);

    Frackit::Python::registerCylinder<double>(module);
    Frackit::Python::registerHollowCylinder<double>(module);
    Frackit::Python::registerBox<double>(module);
    Frackit::Python::registerSphere<double>(module);

    // wrapper classes for OpenCascade BRep shapes
    Frackit::Python::registerBRepWrappers(module);

    // register empty intersection type as geometry class
    Frackit::Python::registerEmptyIntersection<double>(module);

    // distance queries
    Frackit::Python::registerComputeDistance<double>(module);
    Frackit::Python::registerComputeDistanceToBoundary<double>(module);

    // intersection functions
    Frackit::Python::registerIntersectionFunctions<double>(module);
}
