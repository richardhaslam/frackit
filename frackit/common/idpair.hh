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
 * \brief Simple class to define pairs of ids.
 */
#ifndef FRACKIT_ID_PAIR_HH
#define FRACKIT_ID_PAIR_HH

#include <type_traits>

#include "id.hh"

namespace Frackit {

/*!
 * \ingroup Common
 * \brief Simple class to define pairs of ids.
 */
class IdPair
{
public:
    //! Default constructor
    IdPair() = default;

    /*!
     * \brief Construction from indices.
     */
    explicit IdPair(std::size_t id1, std::size_t id2)
    : first_(id1)
    , second_(id2)
    {}

    /*!
     * \brief Construction from ids.
     */
    explicit IdPair(Id id1, Id id2)
    : first_(id1)
    , second_(id2)
    {}

    /*!
     * \brief Retrieve the first id.
     */
    const Id& first() const { return first_; }

    /*!
     * \brief Retrieve the second id.
     */
    const Id& second() const { return second_; }
private:
    Id first_;
    Id second_;
};

} // end namespace Frackit

#endif // FRACKIT_ID_PAIR_HH
