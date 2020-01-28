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
 * \brief Determines the type that results from arithmetic
 *        operations on two different number types.
 */
#ifndef FRACKIT_PROMOTED_TYPE_HH
#define FRACKIT_PROMOTED_TYPE_HH

#include <utility>

namespace Frackit {

/*!
 * \ingroup Common
 * \brief Determines the type that results from arithmetic
 *        operations on two different number types.
 */
template<class T1, class T2>
struct PromotionTraits
{
    using PromotedType = decltype(std::declval<T1>()+std::declval<T2>());
};

/*!
 * \ingroup Common
 * \brief Specialization of the PromotionTraits
 *        for the case of two equal types.
 */
template<class T1>
struct PromotionTraits<T1,T1> { typedef T1 PromotedType; };

/*!
 * \ingroup Common
 * \brief Convenience alias to get the promoted type
 */
template<class T1, class T2>
using PromotedType = typename PromotionTraits<T1, T2>::PromotedType;

} // end namespace Frackit

#endif // FRACKIT_PROMOTED_TYPE_HH
