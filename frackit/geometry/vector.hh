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
 * \brief Classes that implement vectors in n-dimensional space.
 */
#ifndef FRACKIT_GEOMETRY_VECTOR_HH
#define FRACKIT_GEOMETRY_VECTOR_HH

#include <cassert>
#include <cmath>
#include <array>
#include <iostream>
#include <algorithm>
#include <initializer_list>
#include <type_traits>

#include "point.hh"
#include "direction.hh"

namespace Frackit {

/*!
 * \brief Class that implements vectors in
 *        a coordinate space of dimension wd
 * \tparam CT The type used for coordinates
 * \tparam wd The dimension of the coordinate space.
 * \note Specializations for different space dimensions
 *       are provided below.
 */
template<class CT, int wd>
class Vector;

// Forward declarations of other geometry classes
template<class CT, int wd> class Point;
template<class CT, int wd> class Direction;

// namespace with implementation details
namespace Impl {

    /*!
     * \brief Computes the scalar product between two 1d vectors.
     * \tparam CT1 The type used for coordinates of vector 1
     * \tparam CT2 The type used for the coordinatex of vector 2
     *
     * \param v1 the first vector
     * \param v2 the second vector
     * \result the scalar product, cast into coordinate type of vector 1
     */
    template<class CT1, class CT2>
    CT1 scalarProduct(const Vector<CT1, 1>& v1, const Vector<CT2, 1>& v2)
    { return v1.x()*v2.x(); }

    //! Specialization for 2d vectors
    template<class CT1, class CT2>
    CT1 scalarProduct(const Vector<CT1, 2>& v1, const Vector<CT2, 2>& v2)
    { return v1.x()*v2.x() + v1.y()*v2.y(); }

    //! Specialization for 3d vectors
    template<class CT1, class CT2>
    CT1 scalarProduct(const Vector<CT1, 3>& v1, const Vector<CT2, 3>& v2)
    { return v1.x()*v2.x() + v1.y()*v2.y() + v1.z()*v2.z(); }

    /*!
     * \brief Computes the sum of two 1d vectors.
     * \tparam CT1 The type used for coordinates of vector 1
     * \tparam CT2 The type used for the coordinatex of vector 2
     *
     * \param v1 the first vector
     * \param v2 the second vector
     * \result the sum of the vectors, i.e. v1+v2
     * \note We use CT1 as the underlying data type for the result
     */
    template<class CT1, class CT2>
    Vector<CT1, 1> sumVectors(const Vector<CT1, 1>& v1, const Vector<CT2, 1>& v2)
    { return {v1.x()+v2.x()}; }

    //! Specialization for 2d vectors
    template<class CT1, class CT2>
    Vector<CT1, 2> sumVectors(const Vector<CT1, 2>& v1, const Vector<CT2, 2>& v2)
    { return {v1.x()+v2.x(), v1.y()+v2.y()}; }

    //! Specialization for 3d vectors
    template<class CT1, class CT2>
    Vector<CT1, 3> sumVectors(const Vector<CT1, 3>& v1, const Vector<CT2, 3>& v2)
    { return {v1.x()+v2.x(), v1.y()+v2.y(), v1.z()+v2.z()}; }

    /*!
     * \brief Computes the result of subtraction of two vectors.
     * \tparam CT1 The type used for coordinates of vector 1
     * \tparam CT2 The type used for the coordinatex of vector 2
     *
     * \param v1 the first vector
     * \param v2 the second vector
     * \result the vector being the result of the subtraction v1-v2
     * \note We use CT1 as the underlying data type for the result
     */
    template<class CT1, class CT2>
    Vector<CT1, 1> subtractVectors(const Vector<CT1, 1>& v1, const Vector<CT2, 1>& v2)
    { return {v1.x()-v2.x()}; }

    //! Specialization for 2d vectors
    template<class CT1, class CT2>
    Vector<CT1, 2> subtractVectors(const Vector<CT1, 2>& v1, const Vector<CT2, 2>& v2)
    { return {v1.x()-v2.x(), v1.y()-v2.y()}; }

    //! Specialization for 3d vectors
    template<class CT1, class CT2>
    Vector<CT1, 3> subtractVectors(const Vector<CT1, 3>& v1, const Vector<CT2, 3>& v2)
    { return {v1.x()-v2.x(), v1.y()-v2.y(), v1.z()-v2.z()}; }

    /*!
     * \brief Base class of the dimension-specific Vector implementations.
     * \tparam Impl The actual vector class implementation
     * \tparam CT The type used for coordinates
     * \tparam wd The dimension of the space
     */
    template<class Impl, class CoordType, int wd>
    class VectorBase
    {

    public:
        //! export dimensionality
        static constexpr int myDimension() { return 1; };
        static constexpr int worldDimension() { return wd; };

        //! export type used for coordinates
        using ctype = CoordType;

        //! export underlying point type
        using Point = typename Frackit::Point<ctype, wd>;
        using Vector = typename Frackit::Vector<ctype, wd>;
        using Direction = typename Frackit::Direction<ctype, wd>;

        /*!
         * \brief Constructs a vector from a single scalar value
         *        assigning this value to each coordinate direction.
         * \param value The coordinate value
         */
        VectorBase(ctype value = 0.0)
        {
            std::fill(coordinates_.begin(), coordinates_.end(), value);
        }

        /*!
         * \brief Constructs a vector from an initializer list.
         * \param l The initializer list containing the coordinates
         */
        VectorBase(const std::initializer_list<ctype>& l)
        {
            const auto minSize = std::min(static_cast<std::size_t>(wd), l.size());
            std::copy_n( l.begin(), minSize, coordinates_.begin() );
        }

        //! Return the name of this geometry class
        static std::string name()
        { return "Vector"; }

        //! Return the squared length of the vector
        ctype squaredLength() const
        {
            ctype result = 0;
            for (unsigned int i = 0; i < wd; ++i)
                result += coordinates_[i]*coordinates_[i];
            return result;
        }

        //! Return the length of the vector
        ctype length() const
        {
            using std::sqrt;
            return sqrt(squaredLength());
        }

        //! Scalar product with another vector
        template<class CT>
        ctype operator*(const Frackit::Vector<CT, wd>& other) const
        { return scalarProduct(asImp_(), other); }

        //! Add vector to this one
        template<class CT>
        Impl operator+(const Frackit::Vector<CT, wd>& other) const
        { return sumVectors(asImp_(), other); }

        //! Subtract vector from this one
        template<class CT>
        Impl operator-(const Frackit::Vector<CT, wd>& other) const
        { return subtractVectors(asImp_(), other); }

        //! Product of this vector with a scalar
        template<class CT>
        Impl& operator*=(CT value)
        {
            for (auto& coord : coordinates_)
                coord *= value;
            return asImp_();
        }

        //! Division of this vector with a scalar
        template<class CT>
        Impl& operator/=(CT value)
        {
            for (auto& coord : coordinates_)
                coord /= value;
            return asImp_();
        }

    protected:
        //! Provide access to the coordinates of the vector
        ctype operator[] (unsigned int coordIdx) const
        { return coordinates_[coordIdx]; }

        //! Provide access to the coordinates of the vector
        std::array<ctype, wd>& data()
        { return coordinates_; }

        //! Returns the implementation (static polymorphism)
        Impl &asImp_() { return *static_cast<Impl*>(this); }
        //! \copydoc asImp_()
        const Impl &asImp_() const { return *static_cast<const Impl*>(this); }

    private:
        // storage of the coordinates
        std::array<ctype, wd> coordinates_;
    };
} // end namespace Impl

/*!
 * \brief Vector class implementation.
 * \tparam CT The type used for coordinates
 * \tparam wd The dimension of the space
 */
template<class CT, int wd>
class Vector;

/*!
 * \brief Vector class implementation in 1d space.
 */
template<class CT>
class Vector<CT, 1>
: public Impl::VectorBase< Vector<CT, 1>, CT, 1 >
{
    using ThisType = Vector<CT, 1>;
    using ParentType = Impl::VectorBase<ThisType, CT, 1>;

public:
    //! pull up base class' constructors
    using ParentType::ParentType;

    /*!
     * \brief Constructs a vector from two points such that the
     *        resulting vector points from the first to the second point.
     * \param a The first point
     * \param b The second point
     */
    template<class ctype1, class ctype2>
    Vector(const Point<ctype1, 1>& a, const Point<ctype2, 1>& b)
    : ParentType({ b.x() - a.x() })
    {}

    /*!
     * \brief Constructs a vector from a given direction.
     * \param dir The direction
     */
    template<class ctype>
    Vector(const Direction<ctype, 1>& dir)
    : ParentType( {dir.x()} )
    {}

    /*!
     * \brief Constructs a vector with given coordinates.
     * \param x The x-coordinate
     */
    template<class ctype, std::enable_if_t<std::is_arithmetic<ctype>::value, int> = 0>
    Vector(ctype x)
    : ParentType( {x} )
    {}

    //! Return the x-coordinate of the vector
    CT x() const { return (*this)[0]; }
};

/*!
 * \brief Vector class implementation in 2d space.
 */
template<class CT>
class Vector<CT, 2>
: public Impl::VectorBase< Vector<CT, 2>, CT, 2 >
{
    using ThisType = Vector<CT, 2>;
    using ParentType = Impl::VectorBase<ThisType, CT, 2>;

public:
    //! pull up base class' constructors
    using ParentType::ParentType;

    /*!
     * \brief Constructs a vector from two points such that the
     *        resulting vector points from the first to the second point.
     * \param a The first point
     * \param b The second point
     */
    template<class ctype1, class ctype2>
    Vector(const Point<ctype1, 2>& a, const Point<ctype2, 2>& b)
    : ParentType({b.x() - a.x(), b.y() - a.y()})
    {}

    /*!
     * \brief Constructs a vector from a given direction.
     * \param dir The direction
     */
    template<class ctype>
    Vector(const Direction<ctype, 2>& dir)
    : ParentType( {dir.x(), dir.y()} )
    {}

    /*!
     * \brief Constructs a vector with given coordinates.
     * \param x The x-coordinate
     * \param y The y-coordinate
     */
    template<class ctype, std::enable_if_t<std::is_arithmetic<ctype>::value, int> = 0>
    Vector(ctype x, ctype y)
    : ParentType( {x, y} )
    {}

    //! Return the x-coordinate of the vector
    CT x() const { return (*this)[0]; }

    //! Return the y-coordinate of the vector
    CT y() const { return (*this)[1]; }
};

/*!
 * \brief Vector class implementation for 3d space.
 */
template<class CT>
class Vector<CT, 3>
: public Impl::VectorBase< Vector<CT, 3>, CT, 3 >
{
    using ThisType = Vector<CT, 3>;
    using ParentType = Impl::VectorBase<ThisType, CT, 3>;

public:
    //! pull up base class' constructors
    using ParentType::ParentType;

    /*!
     * \brief Constructs a vector from two points such that the
     *        resulting vector points from the first to the second point.
     * \param a The first point
     * \param b The second point
     */
    template<class ctype1, class ctype2>
    Vector(const Point<ctype1, 3>& a, const Point<ctype2, 3>& b)
    : ParentType({b.x() - a.x(), b.y() - a.y(), b.z() - a.z()})
    {}

    /*!
     * \brief Constructs a vector from a given direction.
     * \param dir The direction
     */
    template<class ctype>
    Vector(const Direction<ctype, 3>& dir)
    : ParentType( {dir.x(), dir.y(), dir.z()} )
    {}

    /*!
     * \brief Constructs a vector with given coordinates.
     * \param x The x-coordinate
     * \param y The y-coordinate
     * \param z The z-coordinate
     */
    template<class ctype, std::enable_if_t<std::is_arithmetic<ctype>::value, int> = 0>
    Vector(ctype x, ctype y, ctype z)
    : ParentType( {x, y, z} )
    {}

    //! Return the x-coordinate of the vector
    CT x() const { return (*this)[0]; }

    //! Return the y-coordinate of the vector
    CT y() const { return (*this)[1]; }

    //! Return the z-coordinate of the vector
    CT z() const { return (*this)[2]; }
};

/*!
 * \brief \todo TODO Doc me.
 *  \relates Vector
 *
 *  \param v a vector
 *  \returns a vector orthogonal to v
 */
template<class ctype>
Vector<ctype, 3> makeOrthogonalVector(const Vector<ctype, 3>& v)
{
    if (v.x() == 0.0)      return Vector<ctype, 3>({1.0, 0.0, 0.0});
    else if (v.y() == 0.0) return Vector<ctype, 3>({0.0, 1.0, 0.0});
    else if (v.z() == 0.0) return Vector<ctype, 3>({0.0, 0.0, 1.0});
    else                  return Vector<ctype, 3>({v.y(), - v.x(), 0.0});
}

/*!
 * \brief Writes out a vector in 1d space.
 * \relates Vector
 *
 * \param s std::ostream to write to
 * \param v Vector to write
 *
 * \returns the output stream (s)
 */
template<class CT>
std::ostream& operator<< (std::ostream& s, const Vector<CT, 1>& v)
{
    s << v.x();
    return s;
}

/*!
 * \brief Writes out a vector in 2d space.
 * \relates Vector
 *
 * \param s std::ostream to write to
 * \param v Vector to write
 *
 * \returns the output stream (s)
 */
template<class CT>
std::ostream& operator<< (std::ostream& s, const Vector<CT, 2>& v)
{
    s << v.x() << " " << v.y();
    return s;
}

/*!
 * \brief Writes out a vector in 2d space.
 * \relates Vector
 *
 * \param s std::ostream to write to
 * \param v Vector to write
 *
 * \returns the output stream (s)
 */
template<class CT>
std::ostream& operator<< (std::ostream& s, const Vector<CT, 3>& v)
{
    s << v.x() << " " << v.y() << " " << v.z();
    return s;
}

} // end namespace Frackit

#endif // FRACKIT_GEOMETRY_VECTOR_HH
