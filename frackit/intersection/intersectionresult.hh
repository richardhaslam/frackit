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
#ifndef FRACKIT_INTERSECTION_RESULT_HH
#define FRACKIT_INTERSECTION_RESULT_HH

#include <variant>

namespace Frackit {
namespace IntersectionResult {

//! Obtain the object of an intersection result
template<class Geometry>
Geometry assign(const Geometry& geometry)
{ return geometry; }

} // end namespace IntersectionResult

//! Function to obtain the geometrical object of an intersection
template<class... T>
auto getIntersectionGeometry(const std::variant<T...>& is)
{ return std::visit([&] (auto&& g) { return IntersectionResult::assign(g); }, is); }

} // end Frackit OpenFrack

#endif // FRACKIT_INTERSECTION_RESULT_HH
