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
#ifndef FRACKIT_MATH_HH
#define FRACKIT_MATH_HH

#include <cmath>

#include <frackit/geometry/vector.hh>
#include <frackit/precision/precision.hh>

namespace Frackit {

/*!
 * \brief \todo TODO doc me.
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
 * \todo TODO doc me.
 */
template <class ctype>
ctype crossProduct(const Vector<ctype, 2>& v1,
                   const Vector<ctype, 2>& v2)
{ return v1.x()*v2.y() - v1.y()*v2.x(); }

/*!
 * \brief \todo TODO doc me.
 */
template<class ctype>
ctype boxProduct(const Vector<ctype, 3>& v1,
                 const Vector<ctype, 3>& v2,
                 const Vector<ctype, 3>& v3)
{ return v1*crossProduct(v2, v3); }

/*!
 * \brief \todo TODO doc me.
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
 * \brief \todo TODO doc me.
 */
template<class ctype>
ctype toDegrees(const ctype radians)
{ return radians*180.0/M_PI; }

/*!
 * \brief \todo TODO doc me.
 */
template<class ctype>
ctype toRadians(const ctype degrees)
{ return degrees*M_PI/180.0; }

/*!
 * \brief \todo TODO doc me. Cite Wiki entry.
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
 * \brief \todo TODO doc me. Cite Wiki entry.
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
