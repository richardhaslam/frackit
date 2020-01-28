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
 * \brief Classes that implement points in n-dimensional space.
 */
#ifndef FRACKIT_GEOMETRY_POINT_HH
#define FRACKIT_GEOMETRY_POINT_HH

#include <array>
#include <string>
#include <cassert>
#include <algorithm>
#include <initializer_list>

#include <frackit/precision/precision.hh>
#include "geometry.hh"
#include "vector.hh"

namespace Frackit {

// Forwad declaration
template<class CT, int wd> class Vector;

/*!
 * \ingroup Geometry
 * \brief Base class for Point implementations.
 * \tparam Impl The actual vector class implementation
 * \tparam CT The type used for coordinates
 * \tparam wd The dimension of the space
 */
template<class Impl, class CT, int wd>
class PointBase : public Geometry
{

public:
    //! export dimensionality
    static constexpr int myDimension() { return 0; };
    static constexpr int worldDimension() { return wd; };

    //! export type used for coordinates
    using ctype = CT;

    /*!
     * \brief Default constructor. Creates a point at the origin.
     */
    PointBase()
    {
        std::fill(coordinates_.begin(), coordinates_.end(), 0.0);
    }

    /*!
     * \brief Construction from an initializer list.
     * \param l The initializer list containing the coordinates
     */
    PointBase(const std::initializer_list<CT>& l)
    {
        const auto minDim = std::min(static_cast<std::size_t>(wd), l.size());
        std::copy_n( l.begin(), minDim, coordinates_.begin());
    }

    //! Return the name of this geometry class
    std::string name() const override { return "Point_" + std::to_string(wd) + "d"; }

    //! Returns true if the given point is equal to this one
    bool isEqual(const Impl& other, ctype eps) const
    { return Vector<ctype, wd>(asImp_(), other).squaredLength() < eps*eps; }

    //! Returns true if the given point is equal to this one (default eps)
    bool isEqual(const Impl& other) const
    {
        const auto scale = Vector<ctype, wd>(Impl(), asImp_()).length();
        const auto eps = Precision<ctype>::confusion()*scale;
        return isEqual(other, eps);
    }

protected:
    //! Provide access to the underlying coordinates
    ctype operator[] (unsigned int i) const
    {
        assert(i < coordinates_.size());
        return coordinates_[i];
    }

    //! Provide access to the coordinate storage
    std::array<ctype, wd>& coordinates()
    { return coordinates_; }

    //! Returns the implementation (static polymorphism)
    Impl &asImp_() { return *static_cast<Impl*>(this); }
    //! \copydoc asImp_()
    const Impl &asImp_() const { return *static_cast<const Impl*>(this); }

private:
    std::array<ctype, wd> coordinates_;
};

/*!
 * \ingroup Geometry
 * \brief Point class implementation.
 * \tparam CT The type used for coordinates
 * \tparam wd The dimension of the space
 */
template<class CT, int wd>
class Point;

/*!
 * \ingroup Geometry
 * \brief Point class implementation in 1d space.
 */
template<class CT>
class Point<CT, 1>
: public PointBase< Point<CT, 1>, CT, 1 >
{
    using ThisType = Point<CT, 1>;
    using ParentType = PointBase<ThisType, CT, 1>;

public:
    //! pull up base class' constructors
    using ParentType::ParentType;

    //! Construction from coordinates
    Point(CT x)
    : ParentType({x})
    {}

    //! Return the x-coordinate of the vector
    CT x() const { return (*this)[0]; }

    //! Move this point with the vector v
    template<class C>
    ThisType& operator+= (const Vector<C, 1>& v)
    {
        this->coordinates()[0] += v.x();
        return *this;
    }

    //! Move this point with the vector v
    template<class C>
    ThisType operator+ (const Vector<C, 1>& v) const
    {
        auto tmp(*this);
        tmp += v;
        return tmp;
    }

    //! Move this point with the vector v
    template<class C>
    ThisType& operator-= (const Vector<C, 1>& v)
    {
        this->coordinates()[0] -= v.x();
        return *this;
    }

    //! Move this point with the vector v
    template<class C>
    ThisType operator- (const Vector<C, 1>& v) const
    {
        auto tmp(*this);
        tmp -= v;
        return tmp;
    }
};

/*!
 * \ingroup Geometry
 * \brief Point class implementation in 2d space.
 */
template<class CT>
class Point<CT, 2>
: public PointBase< Point<CT, 2>, CT, 2 >
{
    using ThisType = Point<CT, 2>;
    using ParentType = PointBase<ThisType, CT, 2>;

public:
    //! pull up base class' constructors
    using ParentType::ParentType;

    //! Construction from coordinates
    Point(CT x, CT y)
    : ParentType({x, y})
    {}

    //! Return the x-coordinate of the vector
    CT x() const { return (*this)[0]; }

    //! Return the y-coordinate of the vector
    CT y() const { return (*this)[1]; }

    //! Move this point with the vector v
    template<class C>
    ThisType& operator+= (const Vector<C, 2>& v)
    {
        this->coordinates()[0] += v.x();
        this->coordinates()[1] += v.y();
        return *this;
    }

    //! Move this point with the vector v
    template<class C>
    ThisType operator+ (const Vector<C, 2>& v) const
    {
        auto tmp(*this);
        tmp += v;
        return tmp;
    }

    //! Move this point with the vector v
    template<class C>
    ThisType& operator-= (const Vector<C, 2>& v)
    {
        this->coordinates()[0] -= v.x();
        this->coordinates()[1] -= v.y();
        return *this;
    }

    //! Move this point with the vector v
    template<class C>
    ThisType operator- (const Vector<C, 2>& v) const
    {
        auto tmp(*this);
        tmp -= v;
        return tmp;
    }
};

/*!
 * \ingroup Geometry
 * \brief Point class implementation in 3d space.
 */
template<class CT>
class Point<CT, 3>
: public PointBase< Point<CT, 3>, CT, 3 >
{
    using ThisType = Point<CT, 3>;
    using ParentType = PointBase<ThisType, CT, 3>;

public:
    //! pull up base class' constructors
    using ParentType::ParentType;

    //! Construction from coordinates
    Point(CT x, CT y, CT z)
    : ParentType({x, y, z})
    {}

    //! Return the x-coordinate of the vector
    CT x() const { return (*this)[0]; }

    //! Return the y-coordinate of the vector
    CT y() const { return (*this)[1]; }

    //! Return the z-coordinate of the vector
    CT z() const { return (*this)[2]; }

    //! Move this point with the vector v
    template<class C>
    ThisType& operator+= (const Vector<C, 3>& v)
    {
        this->coordinates()[0] += v.x();
        this->coordinates()[1] += v.y();
        this->coordinates()[2] += v.z();
        return *this;
    }

    //! Move this point with the vector v
    template<class C>
    ThisType operator+ (const Vector<C, 3>& v) const
    {
        auto tmp(*this);
        tmp += v;
        return tmp;
    }

    //! Move this point with the vector v
    template<class C>
    ThisType& operator-= (const Vector<C, 3>& v)
    {
        this->coordinates()[0] -= v.x();
        this->coordinates()[1] -= v.y();
        this->coordinates()[2] -= v.z();
        return *this;
    }

    //! Move this point with the vector v
    template<class C>
    ThisType operator- (const Vector<C, 3>& v) const
    {
        auto tmp(*this);
        tmp -= v;
        return tmp;
    }
};

/*!
 * \brief Writes out a point in 1d space.
 * \relates Point
 *
 * \param s std::ostream to write to
 * \param p Point to write
 *
 * \returns the output stream (s)
 */
template<class CT>
std::ostream& operator<< (std::ostream& s, const Point<CT, 1>& p)
{
    s << p.x();
    return s;
}

/*!
 * \brief Writes out a point in 2d space.
 * \relates Point
 *
 * \param s std::ostream to write to
 * \param p Point to write
 *
 * \returns the output stream (s)
 */
template<class CT>
std::ostream& operator<< (std::ostream& s, const Point<CT, 2>& p)
{
    s << p.x() << " " << p.y();
    return s;
}

/*!
 * \brief Writes out a point in 3d space.
 * \relates Point
 *
 * \param s std::ostream to write to
 * \param p Point to write
 *
 * \returns the output stream (s)
 */
template<class CT>
std::ostream& operator<< (std::ostream& s, const Point<CT, 3>& p)
{
    s << p.x() << " " << p.y() << " " << p.z();
    return s;
}

} // end namespace Frackit

#endif // FRACKIT_GEOMETRY_POINT_HH
