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
 * \ingroup Intersection
 * \brief Functions to check if two axis-aligned bounding boxes
 *        intersect.
 * \todo TODO: implement function to determine the overlap box
 */
#ifndef FRACKIT_BOUNDING_BOX_INTERSECTION_HH
#define FRACKIT_BOUNDING_BOX_INTERSECTION_HH

#include <frackit/common/math.hh>
#include <frackit/precision/defaultepsilon.hh>
#include <frackit/geometry/boundingbox.hh>

namespace Frackit {

/*!
 * \ingroup Intersection
 * \brief Check if two bounding boxes intersect.
 * \tparam CT The type used for coordinates
 * \param bbox1 The first bounding box
 * \param bbox2 The second bounding box
 * \param eps The tolerance value to be used
 */
template<class CT>
bool doIntersect(const BoundingBox<CT>& bbox1, const BoundingBox<CT>& bbox2, CT eps)
{
    if ( !intervalsOverlap({bbox1.xMin(), bbox1.xMax()},
                           {bbox2.xMin(), bbox2.xMax()},
                           eps) )
        return false;

    if ( !intervalsOverlap({bbox1.yMin(), bbox1.yMax()},
                           {bbox2.yMin(), bbox2.yMax()},
                           eps) )
        return false;

    return intervalsOverlap({bbox1.zMin(), bbox1.zMax()},
                            {bbox2.zMin(), bbox2.zMax()},
                            eps);
}

/*!
 * \ingroup Intersection
 * \brief Check if two bounding boxes intersect.
 * \tparam CT The type used for coordinates
 * \param bbox1 The first bounding box
 * \param bbox2 The second bounding box
 * \note This overload chooses a default tolerance value
 */
template<class CT>
bool doIntersect(const BoundingBox<CT>& bbox1, const BoundingBox<CT>& bbox2)
{ return doIntersect(bbox1, bbox2, 0.5*(defaultEpsilon(bbox1) + defaultEpsilon(bbox2))); }

} // end namespace Frackit

#endif // FRACKIT_BOUNDING_BOX_INTERSECTION_HH
