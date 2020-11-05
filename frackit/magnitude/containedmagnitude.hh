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
 * \ingroup Magnitude
 * \brief Contains functionality for computing the
 *        magnitude (length/area/volume) of the part
 *        of a geometry contained within another geometry.
 */
#ifndef FRACKIT_CONTAINED_MAGNITUDE_HH
#define FRACKIT_CONTAINED_MAGNITUDE_HH

#include <type_traits>

#include <frackit/occ/breputilities.hh>
#include <frackit/common/extractdimension.hh>
#include <frackit/precision/defaultepsilon.hh>

#include <frackit/geometry/ctype.hh>
#include <frackit/geometry/geometry.hh>
#include <frackit/geometryutilities/applyongeometry.hh>

#include "magnitude.hh"

namespace Frackit {

/*!
 * \ingroup Magnitude
 * \brief Returns the magnitude of zero-dimensional geometries (points)
 * \note The magnitude of points is always zero independent on if they
 *       are contained in the geometry or not.
 */
template<class Geometry, class Domain,
         std::enable_if_t<DimensionalityTraits<Geometry>::geometryDimension() == 0, int> = 0>
typename CoordinateTypeTraits<Geometry>::type
computeContainedMagnitude(const Geometry& geometry,
                          const Domain& domain)
{ return 0.0; }

/*!
 * \ingroup Magnitude
 * \brief Returns the length of the part of a one-dimensional
 *        geometry that is contained in a domain geometry.
 */
template<class Geometry, class Domain,
         std::enable_if_t<DimensionalityTraits<Geometry>::geometryDimension() == 1, int> = 0>
typename CoordinateTypeTraits<Geometry>::type
computeContainedMagnitude(const Geometry& geometry,
                          const Domain& domain)
{
    const auto geomShape = OCCUtilities::getShape(geometry);
    const auto domainShape = OCCUtilities::getShape(domain);
    const auto is = OCCUtilities::intersect(geomShape, domainShape, defaultEpsilon(domain));
    const auto isEdges = OCCUtilities::getEdges(is);

    if (isEdges.empty())
        return 0.0;

    typename CoordinateTypeTraits<Geometry>::type size = 0.0;
    for (const auto& edge : isEdges)
        size += computeMagnitude(edge);
    return size;
}

/*!
 * \ingroup Magnitude
 * \brief Returns the surface area of the part of a two-dimensional
 *        geometry that is contained in a domain geometry.
 */
template<class Geometry, class Domain,
         std::enable_if_t<DimensionalityTraits<Geometry>::geometryDimension() == 2, int> = 0>
typename CoordinateTypeTraits<Geometry>::type
computeContainedMagnitude(const Geometry& geometry,
                          const Domain& domain)
{
    const auto geomShape = OCCUtilities::getShape(geometry);
    const auto domainShape = OCCUtilities::getShape(domain);
    const auto is = OCCUtilities::intersect(geomShape, domainShape, defaultEpsilon(domain));
    const auto isFaces = OCCUtilities::getFaces(is);

    if (isFaces.empty())
        return 0.0;

    typename CoordinateTypeTraits<Geometry>::type size = 0.0;
    for (const auto& face : isFaces)
    {
        // TODO: Compute center point for TopoDS_Face?
        if constexpr (std::is_same_v<Geometry, TopoDS_Face>)
            size += computeMagnitude(face);
        else
            size += computeMagnitude(face,
                                     defaultEpsilon(domain),
                                     OCCUtilities::point(geometry.center()));
    }

    return size;
}

/*!
 * \ingroup Magnitude
 * \brief Returns the magnitude of that part of a
 *        geometry that is contained in a domamin.
 */
template<class Domain>
typename CoordinateTypeTraits<Domain>::type
computeContainedMagnitude(std::shared_ptr<Geometry> geometry,
                          const Domain& domain)
{
    using ctype = typename CoordinateTypeTraits<Domain>::type;

    // lambda to evaluate the contained magnitude
    auto doComputation = [&] (const auto& geometry) -> ctype
    { return computeContainedMagnitude(geometry, domain); };

    return applyOnGeometry(doComputation, geometry);
}

} // end namespace Frackit

#endif // FRACKIT_CONTAINED_MAGNITUDE_HH
