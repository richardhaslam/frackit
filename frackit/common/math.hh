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
 * \ingroup Common
 * \brief Defines some mathematical functions.
 */
#ifndef FRACKIT_MATH_HH
#define FRACKIT_MATH_HH

#include <cmath>
#include <cassert>
#include <array>
#include <vector>

#include <frackit/precision/precision.hh>
#include <frackit/geometry/vector.hh>
#include <frackit/geometry/direction.hh>

namespace Frackit {

// forward declaration
template<class ctype, int wd> class Vector;
template<class ctype, int wd> class Direction;

/*!
 * \ingroup Common
 * \brief Returns true if the two given intervals overlap.
 * \tparam Scalar The type used for scalar values (floating point type)
 * \param interval1 The first interval
 * \param interval2 The second interval
 * \param eps tolerance used (for floating point types)
 * \note This query includes the boundary, i.e. if the upper limit
 *       of interval1 is equal to the lower limit of interval2,
 *       true is returned.
 * \todo TODO: Overload for integer types without eps (if needed)
 */
template<class Scalar>
bool intervalsOverlap(const std::array<Scalar, 2>& interval1,
                      const std::array<Scalar, 2>& interval2,
                      Scalar eps)
{
    assert(interval1[1] >= interval1[0]);
    assert(interval2[1] >= interval2[0]);

    if (interval1[1] > interval2[1] + eps)
        return interval1[0] < interval2[1] + eps;
    else if (interval1[0] < interval2[0] - eps)
        return interval1[1] > interval2[0] - eps;
    return true;
}

/*!
 * \ingroup Common
 * \brief Returns the vector resulting from the cross product
 *        of two 3-dimensional vectors.
 */
template<class ctype>
Vector<ctype, 3> crossProduct(const Vector<ctype, 3>& v1,
                              const Vector<ctype, 3>& v2)
{
    return { v1.y()*v2.z() - v1.z()*v2.y(),
             v1.z()*v2.x() - v1.x()*v2.z(),
             v1.x()*v2.y() - v1.y()*v2.x() };
}

/*!
 * \ingroup Common
 * \brief Returns the scalar resulting from the cross product
 *        of two 2-dimensional vectors.
 */
template <class ctype>
ctype crossProduct(const Vector<ctype, 2>& v1,
                   const Vector<ctype, 2>& v2)
{ return v1.x()*v2.y() - v1.y()*v2.x(); }

/*!
 * \ingroup Common
 * \brief Returns the scalar resulting from the box product
 *        of three 3-dimensional vectors.
 */
template<class ctype>
ctype boxProduct(const Vector<ctype, 3>& v1,
                 const Vector<ctype, 3>& v2,
                 const Vector<ctype, 3>& v3)
{ return v1*crossProduct(v2, v3); }

/*!
 * \ingroup Common
 * \brief Returns true if the given 3-dimensional
 *        vectors form a right-hand system.
 */
template<class ctype>
bool isRightHandSystem(const Vector<ctype, 3>& v1,
                       const Vector<ctype, 3>& v2,
                       const Vector<ctype, 3>& v3)
{
    assert(std::abs(boxProduct(v1, v2, v3)) > v1.length()*Precision<ctype>::confusion()
                                              *v2.length()*Precision<ctype>::confusion()
                                              *v3.length()*Precision<ctype>::confusion());
    return boxProduct(v1, v2, v3) > 0.0;
}

/*!
 * \ingroup Common
 * \brief Converts a given angle in radians to degrees.
 */
template<class ctype>
ctype toDegrees(const ctype radians)
{ return radians*180.0/M_PI; }

/*!
 * \ingroup Common
 * \brief Converts a given angle in degrees to radians.
 */
template<class ctype>
ctype toRadians(const ctype degrees)
{ return degrees*M_PI/180.0; }

/*!
 * \ingroup Common
 * \brief Rotate a vector around an axis by a given angle.
 * \param v The vector to be rotated
 * \param axis The rotation axis
 * \param angle The rotation angle
 * \note The rotation matrix is taken from:
 *       https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
 */
template<class ctype>
void rotate(Vector<ctype, 3>& v,
            const Direction<ctype, 3>& axis,
            const ctype angle)
{
    using std::sin;
    using std::cos;

    const auto cosPhi = cos(angle);
    const auto sinPhi = sin(angle);

    Vector<ctype, 3> Rx(cosPhi + axis.x()*axis.x()*(1.0 - cosPhi),
                        axis.x()*axis.y()*(1.0 - cosPhi) - axis.z()*sinPhi,
                        axis.x()*axis.z()*(1.0 - cosPhi) - axis.y()*sinPhi);
    Vector<ctype, 3> Ry(axis.y()*axis.x()*(1.0 - cosPhi) + axis.z()*sinPhi,
                        cosPhi + axis.y()*axis.y()*(1.0 - cosPhi),
                        axis.y()*axis.z()*(1.0 - cosPhi) - axis.x()*sinPhi);
    Vector<ctype, 3> Rz(axis.z()*axis.x()*(1.0 - cosPhi) - axis.y()*sinPhi,
                        axis.z()*axis.y()*(1.0 - cosPhi) + axis.x()*sinPhi,
                        cosPhi + axis.z()*axis.z()*(1.0 - cosPhi));
    v = Vector<ctype, 3>(Rx*v, Ry*v, Rz*v);
}

/*!
 * \ingroup Common
 * \brief Rotate several vectors around an axis by a given angle.
 * \param vectors The vectors to be rotated
 * \param axis The rotation axis
 * \param angle The rotation angle
 * \note This overload is more efficient than calling the overload
 *       for a single vector several times, as the rotation matrix
 *       is only constructed once!
 * \note The rotation matrix is taken from:
 *       https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
 */
template<class ctype>
void rotate(std::vector<Vector<ctype, 3>>& vectors,
            const Direction<ctype, 3>& axis,
            const ctype angle)
{
    using std::sin;
    using std::cos;

    const auto cosPhi = cos(angle);
    const auto sinPhi = sin(angle);

    Vector<ctype, 3> Rx(cosPhi + axis.x()*axis.x()*(1.0 - cosPhi),
                        axis.x()*axis.y()*(1.0 - cosPhi) - axis.z()*sinPhi,
                        axis.x()*axis.z()*(1.0 - cosPhi) - axis.y()*sinPhi);
    Vector<ctype, 3> Ry(axis.y()*axis.x()*(1.0 - cosPhi) + axis.z()*sinPhi,
                        cosPhi + axis.y()*axis.y()*(1.0 - cosPhi),
                        axis.y()*axis.z()*(1.0 - cosPhi) - axis.x()*sinPhi);
    Vector<ctype, 3> Rz(axis.z()*axis.x()*(1.0 - cosPhi) - axis.y()*sinPhi,
                        axis.z()*axis.y()*(1.0 - cosPhi) + axis.x()*sinPhi,
                        cosPhi + axis.z()*axis.z()*(1.0 - cosPhi));

    for (auto& v : vectors)
        v = Vector<ctype, 3>(Rx*v, Ry*v, Rz*v);
}

} // end namespace Frackit

#endif // FRACKIT_MATH_HH
