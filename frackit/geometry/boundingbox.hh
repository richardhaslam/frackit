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
 * \ingroup Geometry
 * \brief Class that implements axis-aligned boxes in 3d space.
 */
#ifndef FRACKIT_GEOMETRY_BOUNDING_BOX_HH
#define FRACKIT_GEOMETRY_BOUNDING_BOX_HH

#include "box.hh"

namespace Frackit {

/*!
 * \ingroup Geometry
 * \brief Class that implements axis-aligned boxes in 3d space.
 * \tparam CT The type used for coordinates
 * \note For now this is only an alias for the Box class. But, we introduce
 *       this distinction here to be able to uniquely define interfaces where
 *       axis-aligned boxes are required such that they can easily identified
 *       once the Box class is made more general (to also describe non-axis-aligned
 *       boxes or maybe general hexahedra).
 */
template<class CT>
using BoundingBox = Box<CT>;

} // end namespace Frackit

#endif // FRACKIT_GEOMETRY_BOUNDING_BOX_HH
