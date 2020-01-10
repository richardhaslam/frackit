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
 * \brief Contains functionality for computing the
 *        magnitude (length/area/volume) of the part
 *        of a geometry contained within another geometry.
 */
#ifndef FRACKIT_CONTAINED_MAGNITUDE_HH
#define FRACKIT_CONTAINED_MAGNITUDE_HH

#include <type_traits>

#include <frackit/precision/defaultepsilon.hh>
#include <frackit/occ/breputilities.hh>
#include "magnitude.hh"

namespace Frackit {

/*!
 * \brief Returns the magnitude of zero-dimensional geometries (points)
 * \note The magnitude of points is always zero independent on if they
 *       are contained in the geometry or not.
 */
template<class Geometry, class Domain, std::enable_if_t<Geometry::myDimension() == 0, int> = 0>
typename Geometry::ctype computeContainedMagnitude(const Geometry& geometry,
                                                   const Domain& domain)
{ return 0.0; }

/*!
 * \brief Returns the length of the part of a one-dimensional
 *        geometry that is contained in a domain geometry.
 */
template<class Geometry, class Domain, std::enable_if_t<Geometry::myDimension() == 1, int> = 0>
typename Geometry::ctype computeContainedMagnitude(const Geometry& geometry,
                                                   const Domain& domain)
{
    const auto geomShape = OCCUtilities::getShape(geometry);
    const auto domainShape = OCCUtilities::getShape(domain);
    const auto is = OCCUtilities::intersect(geomShape, domainShape, defaultEpsilon(domain));
    const auto isEdges = OCCUtilities::getEdges(is);

    if (isEdges.empty())
        return 0.0;

    typename Geometry::ctype size = 0.0;
    for (const auto& edge : isEdges)
        size += computeMagnitude(edge);
    return size;
}

/*!
 * \brief Returns the surface area of the part of a two-dimensional
 *        geometry that is contained in a domain geometry.
 */
template<class Geometry, class Domain, std::enable_if_t<Geometry::myDimension() == 2, int> = 0>
typename Geometry::ctype computeContainedMagnitude(const Geometry& geometry,
                                                   const Domain& domain)
{
    const auto geomShape = OCCUtilities::getShape(geometry);
    const auto domainShape = OCCUtilities::getShape(domain);
    const auto is = OCCUtilities::intersect(geomShape, domainShape, defaultEpsilon(domain));
    const auto isFaces = OCCUtilities::getFaces(is);

    if (isFaces.empty())
        return 0.0;

    typename Geometry::ctype size = 0.0;
    for (const auto& face : isFaces)
        size += computeMagnitude(face,
                                 defaultEpsilon(domain),
                                 OCCUtilities::point(geometry.center()));
    return size;
}

} // end namespace Frackit

#endif // FRACKIT_CONTAINED_MAGNITUDE_HH
