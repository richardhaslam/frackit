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
 * \brief Class that implements directions (unit vectors).
 */
#ifndef FRACKIT_GEOMETRY_DIRECTION_HH
#define FRACKIT_GEOMETRY_DIRECTION_HH

#include <string>

#include "geometry.hh"
#include "vector.hh"

namespace Frackit {

// Forward declarations
template<class CT, int wd> class Vector;

/*!
 * \ingroup Geometry
 * \brief Class representing directions (unit vectors).
 * \tparam CT The type used for coordinates
 * \tparam wd The dimension of the space
 */
template<class CT, int wd>
class Direction;

/*!
 * \ingroup Geometry
 * \brief Direction class implementation in 1d space.
 * \tparam CT The type used for coordinates
 */
template<class CT>
class Direction<CT, 1> : public Geometry
{
public:
    //! export underlying vector type
    using Vector = typename Frackit::Vector<CT, 1>;

    //! export dimensionality
    static constexpr int myDimension() { return 1; };
    static constexpr int worldDimension() { return 1; };

    //! default constructor
    Direction() = default;

    /*!
     * \brief Constructs a direction from a vector.
     * \param v The vector to construct a direction from
     */
    template<class ctype>
    Direction(const Frackit::Vector<ctype, 1>& v)
    {
        auto vNorm = v;
        vNorm /= v.length();
        coordinates_[0] = vNorm.x();
    }

    //! Return the name of this geometry
    std::string name() const override { return "Direction_1d"; }

    //! Return the x-coordinate of the direction
    CT x() const { return coordinates_[0]; }

    //! Return the vector a*this
    Vector operator* (CT a) const
    {
        Vector v(x());
        v *= a;
        return v;
    }

    //! Invert the direction
    void invert()
    {
        for (auto& c : coordinates_)
            c *= -1.0;
    }

private:
    std::array<CT, 1> coordinates_;
};

/*!
 * \ingroup Geometry
 * \brief Direction class implementation in 2d space.
 * \tparam CT The type used for coordinates
 */
template<class CT>
class Direction<CT, 2> : public Geometry
{
public:
    //! export underlying vector type
    using Vector = typename Frackit::Vector<CT, 2>;

    //! export dimensionality
    static constexpr int myDimension() { return 1; };
    static constexpr int worldDimension() { return 2; };

    //! default constructor
    Direction() = default;

    /*!
     * \brief Constructs a direction from a vector.
     * \param v The vector to construct a direction from
     */
    template<class ctype>
    Direction(const Frackit::Vector<ctype, 2>& v)
    {
        auto vNorm = v;
        vNorm /= v.length();
        coordinates_[0] = vNorm.x();
        coordinates_[1] = vNorm.y();
    }

    //! Return the name of this geometry
    std::string name() const override { return "Direction_2d"; }

    //! Return the x-coordinate of the direction
    CT x() const { return coordinates_[0]; }

    //! Return the y-coordinate of the direction
    CT y() const { return coordinates_[1]; }

    //! Return the vector a*this
    Vector operator* (CT a) const
    {
        Vector v(x(), y());
        v *= a;
        return v;
    }

    //! Invert the direction
    void invert()
    {
        for (auto& c : coordinates_)
            c *= -1.0;
    }

private:
    std::array<CT, 2> coordinates_;
};

/*!
 * \ingroup Geometry
 * \brief Direction class implementation in 3d space.
 * \tparam CT The type used for coordinates
 */
template<class CT>
class Direction<CT, 3> : public Geometry
{
public:
    //! export underlying vector type
    using Vector = typename Frackit::Vector<CT, 3>;

    //! export dimensionality
    static constexpr int myDimension() { return 1; };
    static constexpr int worldDimension() { return 3; };

    //! default constructor
    Direction() = default;

    /*!
     * \brief Constructs a direction from a vector.
     * \param v The vector to construct a direction from
     */
    template<class ctype>
    Direction(const Frackit::Vector<ctype, 3>& v)
    {
        auto vNorm = v;
        vNorm /= v.length();
        coordinates_[0] = vNorm.x();
        coordinates_[1] = vNorm.y();
        coordinates_[2] = vNorm.z();
    }

    //! Return the name of this geometry
    std::string name() const override { return "Direction_3d"; }

    //! Return the x-coordinate of the direction
    CT x() const { return coordinates_[0]; }

    //! Return the y-coordinate of the direction
    CT y() const { return coordinates_[1]; }

    //! Return the y-coordinate of the direction
    CT z() const { return coordinates_[2]; }

    //! Return the vector a*this
    Vector operator* (CT a) const
    {
        Vector v(x(), y(), z());
        v *= a;
        return v;
    }

    //! Invert the direction
    void invert()
    {
        for (auto& c : coordinates_)
            c *= -1.0;
    }

private:
    std::array<CT, 3> coordinates_;
};

/*!
 * \ingroup Geometry
 * \brief Writes out a direction in 1d space.
 * \relates Direction
 *
 * \param s std::ostream to write to
 * \param d Direction to write
 *
 * \returns the output stream (s)
 */
template<class CT>
std::ostream& operator<< (std::ostream& s, const Direction<CT, 1>& d)
{
    s << d.x();
    return s;
}

/*!
 * \ingroup Geometry
 * \brief Writes out a direction in 2d space.
 * \relates Direction
 *
 * \param s std::ostream to write to
 * \param d Direction to write
 *
 * \returns the output stream (s)
 */
template<class CT>
std::ostream& operator<< (std::ostream& s, const Direction<CT, 2>& d)
{
    s << d.x() << " " << d.y();
    return s;
}

/*!
 * \ingroup Geometry
 * \brief Writes out a direction in 2d space.
 * \relates Direction
 *
 * \param s std::ostream to write to
 * \param d Direction to write
 *
 * \returns the output stream (s)
 */
template<class CT>
std::ostream& operator<< (std::ostream& s, const Direction<CT, 3>& d)
{
    s << d.x() << " " << d.y() << " " << d.z();
    return s;
}

} // end namespace Frackit

#endif // FRACKIT_GEOMETRY_DIRECTION_HH
