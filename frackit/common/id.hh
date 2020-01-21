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
 * \brief Simple wrapper class to define ids/indices.
 */
#ifndef FRACKIT_ID_HH
#define FRACKIT_ID_HH

namespace Frackit {

/*!
 * \brief Simple wrapper class to define ids/indices.
 *        This can be used wherever indices are passed
 *        to interfaces to avoid implicit conversion of
 *        function arguments.
 */
class Id
{
public:
    /*!
     * \brief Construction from an index.
     */
    explicit Id(std::size_t id)
    : id_(id)
    {}

    /*!
     * \brief Retrieve the index.
     */
    std::size_t get() const
    { return id_; }

private:
    std::size_t id_;
};

} // end namespace Frackit

#endif // FRACKIT_ID_HH
