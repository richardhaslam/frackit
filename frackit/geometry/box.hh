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
 * \brief Class that implements axis-aligned boxes in 3d space.
 */
#ifndef FRACKIT_GEOMETRY_BOX_HH
#define FRACKIT_GEOMETRY_BOX_HH

#include <cassert>
#include <stdexcept>
#include <string>

#include <frackit/precision/precision.hh>
#include "geometry.hh"
#include "point.hh"
#include "segment.hh"
#include "quadrilateral.hh"
#include "vector.hh"

namespace Frackit {

/*!
 * \ingroup Geometry
 * \brief Class that implements axis-aligned boxes in 3d space.
 * \tparam CT The type used for coordinates
 */
template<class CT>
class Box : public Geometry
{
    using Vector = Frackit::Vector<CT, 3>;

public:
    //! export dimensionality
    static constexpr int myDimension() { return 3; };
    static constexpr int worldDimension() { return 3; };

    //! export type used for coordinates
    using ctype = CT;

    //! export underlying geometry types
    using Point = Frackit::Point<ctype, 3>;
    using Quadrilateral = Frackit::Quadrilateral<ctype, 3>;
    using Segment = Frackit::Segment<ctype, 3>;

    /*!
     * \brief Construct a box from the minimum/maximum coordinates.
     */
    Box(ctype xmin, ctype ymin, ctype zmin,
        ctype xmax, ctype ymax, ctype zmax)
    : xMin_(xmin), xMax_(xmax)
    , yMin_(ymin), yMax_(ymax)
    , zMin_(zmin), zMax_(zmax)
    {
        if (xMax_ < xMin_) throw std::runtime_error( std::string("Could not construct box: xMax < xMin!") );
        if (yMax_ < yMin_) throw std::runtime_error( std::string("Could not construct box: ymax < ymin!") );
        if (zMax_ < zMin_) throw std::runtime_error( std::string("Could not construct box: zmax < zMin!") );
    }

    /*!
     * \brief Construct a box from the minimum/maximum corner points.
     */
    Box(const Point& pMin, const Point& pMax)
    : Box(pMin.x(), pMin.y(), pMin.z(),
          pMax.x(), pMax.y(), pMax.z())
    {}

    //! Return the name of this geometry.
    std::string name() const override { return "Box"; }
    //! Return the x-coordinate of the first corner.
    ctype xMin() const { return xMin_; }
    //! Return the y-coordinate of the first corner.
    ctype yMin() const { return yMin_; }
    //! Return the z-coordinate of the first corner.
    ctype zMin() const { return zMin_; }
    //! Return the x-coordinate of the last corner.
    ctype xMax() const { return xMax_; }
    //! Return the y-coordinate of the last corner.
    ctype yMax() const { return yMax_; }
    //! Return the z-coordinate of the last corner.
    ctype zMax() const { return zMax_; }

    //! Return the volume of the box
    ctype volume() const
    {
        return (xMax() - xMin())
               *(yMax() - yMin())
               *(zMax() - zMin());
    }

    //! Return the number of corners
    static constexpr std::size_t numCorners()
    { return 8; }

    //! Return the corner for the given index
    Point corner(unsigned int cornerIdx) const
    {
        assert(cornerIdx < numCorners());
        switch (cornerIdx)
        {
            case 0: return Point({xMin_, yMin_, zMin_});
            case 1: return Point({xMax_, yMin_, zMin_});
            case 2: return Point({xMin_, yMax_, zMin_});
            case 3: return Point({xMax_, yMax_, zMin_});
            case 4: return Point({xMin_, yMin_, zMax_});
            case 5: return Point({xMax_, yMin_, zMax_});
            case 6: return Point({xMin_, yMax_, zMax_});
            case 7: return Point({xMax_, yMax_, zMax_});
            default: throw std::runtime_error( std::string("Box: Invalid corner index!") );
        }
    }

    //! Return the number of edges.
    static constexpr std::size_t numEdges()
    { return 12; }

    //! Return the edge for the given index
    Segment edge(unsigned int edgeIdx) const
    {
        assert(edgeIdx < numEdges());
        switch (edgeIdx)
        {
            case 0: return Segment(corner(0), corner(1));
            case 1: return Segment(corner(0), corner(2));
            case 2: return Segment(corner(1), corner(3));
            case 3: return Segment(corner(2), corner(3));
            case 4: return Segment(corner(0), corner(4));
            case 5: return Segment(corner(1), corner(5));
            case 6: return Segment(corner(2), corner(6));
            case 7: return Segment(corner(3), corner(7));
            case 8: return Segment(corner(4), corner(5));
            case 9: return Segment(corner(4), corner(6));
            case 10: return Segment(corner(5), corner(7));
            case 11: return Segment(corner(6), corner(7));
            default: throw std::runtime_error( std::string("Box: Invalid edge index!") );
        }
    }

    //! Return the number of faces
    static constexpr std::size_t numFaces()
    { return 6; }

    //! Return the face for the given index
    Quadrilateral face(unsigned int faceIdx) const
    {
        assert(faceIdx < numFaces());
        switch (faceIdx)
        {
            case 0: return Quadrilateral(corner(0), corner(1), corner(2), corner(3)); // bottom
            case 1: return Quadrilateral(corner(0), corner(1), corner(4), corner(5)); // front
            case 2: return Quadrilateral(corner(0), corner(2), corner(4), corner(6)); // left
            case 3: return Quadrilateral(corner(1), corner(3), corner(5), corner(7)); // right
            case 4: return Quadrilateral(corner(2), corner(3), corner(6), corner(7)); // back
            case 5: return Quadrilateral(corner(4), corner(5), corner(6), corner(7)); // top
            default: throw std::runtime_error( std::string("Box: Invalid face index!") );
        }
    }

    /*!
     * \brief Returns true if a point lies within the box.
     * \param p The point to be checked
     * \param eps The epsilon (tolerance) value to be used
     */
    bool contains(const Point& p, ctype eps) const
    {
        return p.x() > xMin() - eps && p.x() < xMax() + eps
               && p.y() > yMin() - eps && p.y() < yMax() + eps
               && p.z() > zMin() - eps && p.z() < zMax() + eps;
    }

    /*!
     * \brief Returns true if a point lies within the box.
     * \param p The point to be checked
     * \note This overload uses a default epsilon.
     */
    bool contains(const Point& p) const
    {
        auto eps = Precision<ctype>::confusion();
        eps *= Vector(corner(0), corner(7)).length();
        return contains(p, eps);
    }

private:
    ctype xMin_; ctype xMax_;
    ctype yMin_; ctype yMax_;
    ctype zMin_; ctype zMax_;
};

} // end namespace Frackit

#endif // FRACKIT_GEOMETRY_BOX_HH
