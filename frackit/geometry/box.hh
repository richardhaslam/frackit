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
#ifndef FRACKIT_GEOMETRY_BOX_HH
#define FRACKIT_GEOMETRY_BOX_HH

#include <cassert>
#include <stdexcept>

#include "point.hh"
#include "segment.hh"
#include "quadrilateral.hh"

namespace Frackit {

/*!
 * \brief \todo TODO doc me.
 */
template<class CT>
class Box
{

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
     * \brief \todo TODO doc me.
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

    //! \todo TODO doc me.
    static std::string name() { return "Box"; }
    //! \todo TODO doc me.
    ctype xMin() const { return xMin_; }
    //! \todo TODO doc me.
    ctype yMin() const { return yMin_; }
    //! \todo TODO doc me.
    ctype zMin() const { return zMin_; }
    //! \todo TODO doc me.
    ctype xMax() const { return xMax_; }
    //! \todo TODO doc me.
    ctype yMax() const { return yMax_; }
    //! \todo TODO doc me.
    ctype zMax() const { return zMax_; }

    //! \todo TODO doc me.
    ctype volume() const
    {
        return (xMax() - xMin())
               *(yMax() - yMin())
               *(zMax() - zMin());
    }

    //! \todo TODO doc me.
    static constexpr std::size_t numCorners()
    { return 8; }

    //! \todo TODO doc me.
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

    //! \todo TODO doc me.
    static constexpr std::size_t numEdges()
    { return 12; }

    //! \todo TODO doc me.
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

    //! \todo TODO doc me.
    static constexpr std::size_t numFaces()
    { return 6; }

    //! \todo TODO doc me.
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

    //! Returns true if a point lies within the box
    bool contains(const Point& p, ctype eps) const
    {
        return p.x() > xMin() - eps && p.x() < xMax() + eps
               && p.y() > yMin() - eps && p.y() < yMax() + eps
               && p.z() > zMin() - eps && p.z() < zMax() + eps;
    }

    //! Returns true if a point lies within the box
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
