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
 * \ingroup Distance
 * \brief Contains functionality for computing
 *        the distance between geometries.
 */
#ifndef FRACKIT_DISTANCE_HH
#define FRACKIT_DISTANCE_HH

#include <cmath>
#include <limits>
#include <vector>
#include <variant>
#include <stdexcept>
#include <type_traits>

#include <TopoDS_Shape.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <Extrema_ExtAlgo.hxx>
#include <Extrema_ExtFlag.hxx>

#include <frackit/precision/precision.hh>
#include <frackit/common/promotedtype.hh>
#include <frackit/occ/breputilities.hh>

#include <frackit/geometry/point.hh>
#include <frackit/geometry/vector.hh>
#include <frackit/geometry/line.hh>
#include <frackit/geometry/segment.hh>

namespace Frackit {
namespace Impl {
    // Convenience alias to extract the promoted type
    // of the types used for coordinates of two geometries.
    // This makes the code below more readable in some places
    template<class Geom1, class Geom2>
    using PCT = PromotedType<typename Geom1::ctype,
                             typename Geom2::ctype>;
}

/*!
 * \ingroup Distance
 * \brief Computes the distance between two TopoDS_Shape objects.
 * \param shape1 The first shape
 * \param shape2 The second shape
 * \param deflection The epsilon used in the BrepExtrema command
 * \param extFlag The flag passed to the BrepExtrema command (MIN/MAX/MINMAX)
 * \param extAlgo The algorithm passed to the BrepExtrema command (TREE/GRAD)
 */
template<class ctype = double>
ctype computeDistance(const TopoDS_Shape& shape1,
                      const TopoDS_Shape& shape2,
                      ctype deflection = Precision<ctype>::confusion(),
                      Extrema_ExtFlag extFlag = Extrema_ExtFlag_MINMAX,
                      Extrema_ExtAlgo extAlgo = Extrema_ExtAlgo_Grad)
{
    BRepExtrema_DistShapeShape algo(shape1, shape2, deflection, extFlag, extAlgo);
    if (!algo.IsDone())
        throw std::runtime_error("Could not compute BRep distance");
    return algo.Value();
}

/*!
 * \ingroup Distance
 * \brief In the general case we compute the distance based on
 *        the basis of the BRep of the geometries. Overloads for
 *        geometries for which the distance can be computed more
 *        easily are provided below.
 * \param geo1 The first geometry
 * \param geo2 The second geometry
 * \param deflection The epsilon used in the BrepExtrema command
 * \param extFlag The flag passed to the BrepExtrema command (MIN/MAX/MINMAX)
 * \param extAlgo The algorithm passed to the BrepExtrema command (TREE/GRAD)
 * \note This overload fails if the provided geometry types do not export
 *       the type "ctype" as the promoted type, i.e. the call to Impl::PCT<,>
 *       requires that. Overloads for one of the geometries being a class of
 *       the TopoDS package are provided below.
 */
template<class Geom1, class Geom2>
Impl::PCT<Geom1, Geom2>
computeDistance(const Geom1& geo1,
                const Geom2& geo2,
                Impl::PCT<Geom1, Geom2> deflection = Precision<Impl::PCT<Geom1, Geom2>>::confusion(),
                Extrema_ExtFlag extFlag = Extrema_ExtFlag_MINMAX,
                Extrema_ExtAlgo extAlgo = Extrema_ExtAlgo_Grad)
{
    return computeDistance(OCCUtilities::getShape(geo1),
                           OCCUtilities::getShape(geo2),
                           deflection,
                           extFlag,
                           extAlgo);
}

/*!
 * \ingroup Distance
 * \brief Computes the distance between a TopoDS_Shape and an internal geometry type.
 * \param shape The shape of one geometry
 * \param geo The second geometry
 * \param deflection The epsilon used in the BrepExtrema command
 * \param extFlag The flag passed to the BrepExtrema command (MIN/MAX/MINMAX)
 * \param extAlgo The algorithm passed to the BrepExtrema command (TREE/GRAD)
 */
template<class Geom>
typename Geom::ctype
computeDistance(const TopoDS_Shape& shape,
                const Geom& geo,
                typename Geom::ctype deflection = Precision<typename Geom::ctype>::confusion(),
                Extrema_ExtFlag extFlag = Extrema_ExtFlag_MINMAX,
                Extrema_ExtAlgo extAlgo = Extrema_ExtAlgo_Grad)
{ return computeDistance(shape, OCCUtilities::getShape(geo), deflection, extFlag, extAlgo); }

/*!
 * \ingroup Distance
 * \brief Computes the distance between an internal geometry type and a TopoDS_Shape.
 * \param geo The second geometry
 * \param shape The shape of one geometry
 * \param deflection The epsilon used in the BrepExtrema command
 * \param extFlag The flag passed to the BrepExtrema command (MIN/MAX/MINMAX)
 * \param extAlgo The algorithm passed to the BrepExtrema command (TREE/GRAD)
 */
template<class Geom>
typename Geom::ctype
computeDistance(const Geom& geo,
                const TopoDS_Shape& shape,
                typename Geom::ctype deflection = Precision<typename Geom::ctype>::confusion(),
                Extrema_ExtFlag extFlag = Extrema_ExtFlag_MINMAX,
                Extrema_ExtAlgo extAlgo = Extrema_ExtAlgo_Grad)
{ return computeDistance(OCCUtilities::getShape(geo), shape, deflection, extFlag, extAlgo); }

/*!
 * \ingroup Distance
 * \brief Returns the euclidian distance between two points.
 * \param p1 The first point
 * \param p2 The second point
 */
template<class ctype1, class ctype2, int worldDim>
PromotedType<ctype1, ctype2> computeDistance(const Point<ctype1, worldDim>& p1,
                                             const Point<ctype2, worldDim>& p2)
{ return Vector<PromotedType<ctype1, ctype2>, worldDim>(p1, p2).length(); }

/*!
 * \ingroup Distance
 * \brief Returns the euclidian distance between a point and a line.
 * \param p The point
 * \param line The line
 */
template<class ctype1, class ctype2, int worldDim>
PromotedType<ctype1, ctype2> computeDistance(const Point<ctype1, worldDim>& p,
                                             const Line<ctype2, worldDim>& line)
{ return Vector<PromotedType<ctype1, ctype2>, worldDim>(p, line.projection(p)).length(); }

/*!
 * \ingroup Distance
 * \brief Returns the euclidian distance between a line and a point.
 * \param line The line
 * \param p The point
 */
template<class ctype1, class ctype2, int worldDim>
PromotedType<ctype1, ctype2> computeDistance(const Line<ctype1, worldDim>& line,
                                             const Point<ctype2, worldDim>& p)
{ return computeDistance(p, line); }

/*!
 * \ingroup Distance
 * \brief Returns the euclidian distance between a point and a segment.
 * \param p The point
 * \param seg The segment
 */
template<class ctype1, class ctype2, int worldDim>
PromotedType<ctype1, ctype2> computeDistance(const Point<ctype1, worldDim>& p,
                                             const Segment<ctype2, worldDim>& seg)
{
    using ctype = PromotedType<ctype1, ctype2>;

    const auto proj = seg.supportingLine().projection(p);
    if (seg.contains(proj))
        return Vector<ctype, worldDim>(p, proj).length();

    const auto d1 = Vector<ctype, worldDim>(p, seg.source()).length();
    const auto d2 = Vector<ctype, worldDim>(p, seg.target()).length();

    using std::min;
    return min(d1, d2);
}

/*!
 * \ingroup Distance
 * \brief Returns the euclidian distance between a segment and a point.
 * \param seg The segment
 * \param p The point
 */
template<class ctype1, class ctype2, int worldDim>
PromotedType<ctype1, ctype2> computeDistance(const Segment<ctype1, worldDim>& seg,
                                             const Point<ctype2, worldDim>& p)
{ return computeDistance(p, seg); }

/*!
 * \ingroup Distance
 * \brief Returns the minimum distance between a geometry and
 *        std::variant, for instance, resulting from intersection
 *        of two geometries.
 * \param geo The first geometry
 * \param is The intersection geometry
 */
template<class Geo, class... Geos>
auto computeDistance(const Geo& geo1,
                     const std::variant<Geos...>& is)
{ return std::visit([&] (const auto& isGeo) { return computeDistance(geo1, isGeo); }, is); }

/*!
 * \ingroup Distance
 * \brief Returns the minimum distance between a geometry and
 *        vector of geometries.
 * \param geo1 The first geometry
 * \param geos2 The vector of geometries
 * \return the minimum of the distances with the individual
 *         geometries provided in geos2.
 */
template<class Geo1, class Geo2>
auto computeDistance(const Geo1& geo1,
                     const std::vector<Geo2>& geos2)
{
    using std::min;
    using ctype = std::decay_t<decltype(computeDistance(geo1, geos2[0]))>;

    ctype minDist = std::numeric_limits<ctype>::max();
    for (const auto& geo2 : geos2)
        minDist = min(minDist, computeDistance(geo1, geo2));

    return minDist;
}

} // end namespace Frackit

#endif // FRACKIT_DISTANCE_HH
