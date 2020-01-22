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
 * \brief Interface for classes that randomly generate points.
 */
#ifndef FRACKIT_POINT_SAMPLER_HH
#define FRACKIT_POINT_SAMPLER_HH

namespace Frackit {

/*!
 * \brief Defines the interface for point samplers.
 * \tparam P The type used for points
 */
template< class P >
class PointSampler
{
public:
    //! Export the sampled point type
    using Point = P;

    //! every abstract base class has a virtual destructor
    virtual ~PointSampler () {}

    /*!
     * \brief Return a point.
     */
    virtual Point operator() () = 0;
};

} // end namespace Frackit

#endif // FRACKIT_POINT_SAMPLER_HH
