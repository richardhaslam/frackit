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
 * \ingroup Common
 * \brief Defines classes to obtain type information at compile-time.
 */
#ifndef FRACKIT_COMMON_TYPE_TRAITS_HH
#define FRACKIT_COMMON_TYPE_TRAITS_HH

#include <type_traits>

namespace Frackit {

/*!
 * \ingroup Common
 * \brief Helper struct to detect if a type T
 *        is contained in a parameter pack Ts
 * \tparam T The type of which an ocurrence in the pack is to be checked
 * \tparam Ts The parameter pack
 */
template<class T, class... Ts>
struct Contains : std::disjunction<std::is_same<T, Ts>...> {};

} // end namespace Frackit

#endif // FRACKIT_COMMON_TYPE_TRAITS_HH
