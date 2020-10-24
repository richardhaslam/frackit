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
 * \ingroup Geometry
 * \brief Classes that implement polygons in n-dimensional space.
 */
#ifndef FRACKIT_GEOMETRY_POLYGON_HH
#define FRACKIT_GEOMETRY_POLYGON_HH

#include <cassert>
#include <cmath>
#include <string>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <iterator>

#include <frackit/common/math.hh>
#include <frackit/precision/precision.hh>

#include "geometry.hh"
#include "point.hh"
#include "vector.hh"
#include "segment.hh"
#include "triangle.hh"
#include "plane.hh"
#include "triangulation.hh"

namespace Frackit {

/*!
 * \ingroup Geometry
 * \brief Class that implements polygons
 *        in a coordinate space with the dimension worldDim.
 * \tparam CT The type used for coordinates
 * \tparam worldDim The dimension of the coordinate space
 */
template<class CT, int worldDim>
class Polygon;

/*!
 * \ingroup Geometry
 * \brief Implementation of polygons in 3d space.
 * \tparam CT The type used for coordinates
 */
template<class CT>
class Polygon<CT, 3> : public Geometry
{
    using Triangle = Frackit::Triangle<CT, 3>;
    using Vector = Frackit::Vector<CT, 3>;

public:
    //! export dimensionality
    static constexpr int myDimension() { return 2; };
    static constexpr int worldDimension() { return 3; };

    //! export type used for coordinates
    using ctype = CT;

    //! export underlying geometry types
    using Point = Frackit::Point<ctype, 3>;
    using Segment = Frackit::Segment<ctype, 3>;
    using Plane = Frackit::Plane<ctype, 3>;

    //! Default constructor
    Polygon() = default;

    /*!
     * \brief The constructor.
     * \param corners The corner points of the polygon
     * \note We expect the corners to be sorted, that is,
     *       that consecutive corners in the provided container
     *       are neighboring corners of the polygon. Moreover,
     *       we expect the corners to be convertible to PointType.
     */
    template<class CornerStorage>
    Polygon(const CornerStorage& corners)
    {
        std::copy(corners.begin(), corners.end(), std::back_inserter(corners_));
        if (corners_.size() < 3)
            throw std::runtime_error("At least 3 corners are required!");

        // determine if polygon is convex
        determineConvexity_();

        // compute area & center
        const auto triangulation = triangulate(corners_);
        supportPlane_ = Plane(corners_[0], corners_[1], corners_[2]);
        center_ = Point(0.0, 0.0, 0.0);
        area_ = 0.0;

        for (const auto& t : triangulation)
        {
            auto d = Vector(Point(0.0, 0.0, 0.0), t.center());
            d *= t.area();

            center_ += d;
            area_ += t.area();
        }
        center_ /= area_;

        // check if polygon is degenerated (e.g. because points are unsorted)
        assert( !isDegenerated_(triangulation) );
        assert( std::all_of(corners_.begin(), corners_.end(),
        [&] (const auto& c)
        {
            auto eps = Precision<ctype>::confusion();
            eps *= Vector(center_, corners_[0]).length();
            return supportPlane_.contains(c, eps);
        }) );
    }

    //! Return the name of this geometry
    std::string name() const override { return "Polygon_3d"; }
    //! Return the number of corners
    std::size_t numCorners() const { return corners_.size(); }
    //! Return the i-th corner
    const Point& corner(unsigned int i) const
    {
        assert(i < numCorners());
        return corners_[i];
    }

    //! Return the number of edges
    std::size_t numEdges() const { return corners_.size(); }
    //! Return the i-th edge
    Segment edge(unsigned int i) const
    {
        assert(i < numEdges());
        return Segment(corner(i), corner((i+1)%numEdges()));
    }

    //! Return the center of the polygon
    const Point& center() const { return center_; }
    //! Return the plane this polygon is embedded in
    const Plane& supportingPlane() const { return supportPlane_; }
    //! returns true if the polygon is convex
    bool isConvex() const { return isConvex_; }
    //! Return the area of the polygon
    ctype area() const
    {
        if (!isConvex())
            throw std::runtime_error("area() for non-convex polygons not implemented");
        return area_;
    }

    /*!
     * \brief Returns true if a point is on the polygon
     * \param p The point to be checked
     * \param eps The tolerance to be used
     * \param checkIfOnPlane This can be set to false in case one
     *                       knows that the point lies on the supporting
     *                       plane in order to skip this check at this point.
     */
    bool contains(const Point& p, ctype eps, bool checkIfOnPlane = true) const
    {
        if (isConvex_)
            return convexContains(p, eps, checkIfOnPlane);
        else
            throw std::runtime_error("contains() query for non-convex polygons");
    }

    /*!
     * \brief Returns true if a point is on the polygon
     * \param p The point to be checked
     * \param checkIfOnPlane This can be set to false in case one
     *                       knows that the point lies on the supporting
     *                       plane in order to skip this check at this point.
     * \note This overload uses a default epsilon
     */
    bool contains(const Point& p,  bool checkIfOnPlane = true) const
    {
        auto eps = Precision<ctype>::confusion();
        eps *= Vector(corner(0), center_).length();
        return contains(p, eps, checkIfOnPlane);
    }

private:
    /*!
     * \brief Contains query for convex polygons
     * \todo TODO: There must be an easier way to achieve this?
     */
    bool convexContains(const Point& p, ctype eps, bool checkIfOnPlane) const
    {
        if (checkIfOnPlane)
            if (!supportingPlane().contains(p, eps))
                return false;

        // sum up the areas of sub-triangles with the given point
        ctype otherArea(0.0);
        for (std::size_t i = 0; i < numCorners()-1; ++i)
            otherArea += Triangle(corner(i), corner(i+1), p).area();
        otherArea += Triangle(corner(numCorners()-1), corner(0), p).area();

        // if area is equal to the polygon area, point is on it
        using std::abs;
        return abs(area() - otherArea) < eps*eps;
    }

    //! returns true if the polygon, represented by its triangulation, is degenerated
    template<class Triangulation>
    bool isDegenerated_(const Triangulation& triangulation) const
    {
        // TODO: implement
        return false;
    }

    //! determines if the polygon is convex or not
    void determineConvexity_()
    {
        isConvex_ = true;

        // use normal over first corner as test vector
        const auto testVector = crossProduct(Vector(edge(numEdges()-1).direction()),
                                             Vector(edge(0).direction()));

        // check if normals around the other edges point in the same direction
        for (std::size_t i = 1; i < numEdges()-1; ++i)
        {
            const auto n = crossProduct(Vector(edge(i).direction()),
                                        Vector(edge(i+1).direction()));

            // angle changed, polygon is concave
            if (std::signbit(n*testVector))
            { isConvex_ = false; return; }
        }
    }

    std::vector<Point> corners_;
    Plane supportPlane_;
    Point center_;
    ctype area_;
    bool isConvex_;
};

} // end namespace Frackit

#endif // FRACKIT_GEOMETRY_POLYGON_HH
