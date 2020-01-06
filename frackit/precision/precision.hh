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
 * \brief Defines base epsilons to be used for floating point arithmetic.
 */
#ifndef FRACKIT_PRECISION_HH
#define FRACKIT_PRECISION_HH

// standard tolerances for floating point comparison
#include <Precision.hxx>

namespace Frackit {

//! Alias for the OpenCascade Precision class
using OCCPrecision = Precision;

/*!
 * \brief Defines base epsilons to be used for floating point arithmetic.
 */
template<class ctype>
struct Precision
{
    /*!
     * \brief Tolerance value to be used for equality queries.
     */
    static ctype confusion() { return OCCPrecision::Confusion(); }

    /*!
     * \brief Tolerance value to be used for angle comparisons.
     *        For example, to determine if two vectors are parallel.
     * \note OpenCascade defines an angular precision of 1e-12, but
     *       in unit tests we saw that this might be too restrictive.
     */
    static ctype angular() { return OCCPrecision::Angular()*1e2; }
};

} // end namespace Frackit

#endif // FRACKIT_PRECISION_HH
