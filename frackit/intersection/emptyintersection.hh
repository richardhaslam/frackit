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
 * \brief Class used to define empty intersection.
 */
#ifndef FRACKIT_EMPTY_INTERSECTION_HH
#define FRACKIT_EMPTY_INTERSECTION_HH

#include <string>

namespace Frackit {

/*!
 * \brief Class used to define empty intersection.
 * \tparam wd The dimension of the coordinate space.
 * \tparam CT The type used for coordinates
 */
template<int wd = 3, class CT = double>
struct EmptyIntersection
{
    //! Technically, this does not have a defined dimension
    static constexpr int myDimension() { return 0; }

    //! Return the dimension of the coordinate space
    static constexpr int worldDimension() { return wd; }

    //! Return the name of this geometry
    static std::string name() { return "EmptyIntersection"; }

    //! Export coordinate type for compatibility
    using ctype = CT;
};

} // end namespace Frackit

#endif // FRACKIT_EMPTY_INTERSECTION_HH
