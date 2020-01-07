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
 * \brief Contains functionality for computing
 *        the distance between geometries.
 */
#ifndef FRACKIT_DISTANCE_HH
#define FRACKIT_DISTANCE_HH

#include <cmath>

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
    using PCT = PromotedType<typename Geom1::ctype, typename Geom2::ctype>;
}

/*!
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
        throw std::runtime_error(std::string("Could not compute BRep distance"));

    return algo.Value();
}

/*!
 * \brief In the general case we compute the distance based on
 *        the basis of the BRep of the geometries. Overloads for
 *        geometries for which the distance can be computed more
 *        easily are provided below.
 * \param geo1 The first geometry
 * \param geo2 The second geometry
 * \param deflection The epsilon used in the BrepExtrema command
 * \param extFlag The flag passed to the BrepExtrema command (MIN/MAX/MINMAX)
 * \param extAlgo The algorithm passed to the BrepExtrema command (TREE/GRAD)
 */
template<class Geom1, class Geom2>
Impl::PCT<Geom1, Geom2>
computeDistance(const Geom1& geo1,
                const Geom2& geo2,
                Impl::PCT<Geom1, Geom2> deflection = Precision<Impl::PCT<Geom1, Geom2>>::confusion(),
                Extrema_ExtFlag extFlag = Extrema_ExtFlag_MINMAX,
                Extrema_ExtAlgo extAlgo = Extrema_ExtAlgo_Grad)
{ return computeDistance(OCCUtilities::getShape(geo1), OCCUtilities::getShape(geo2)); }

/*!
 * \brief Returns the euclidian distance between two points.
 */
template<class ctype1, class ctype2, int worldDim>
PromotedType<ctype1, ctype2> computeDistance(const Point<ctype1, worldDim>& p1,
                                             const Point<ctype2, worldDim>& p2)
{ return Vector<PromotedType<ctype1, ctype2>, worldDim>(p1, p2).length(); }

/*!
 * \brief Returns the euclidian distance between a point and a line.
 */
template<class ctype1, class ctype2, int worldDim>
PromotedType<ctype1, ctype2> computeDistance(const Point<ctype1, worldDim>& p,
                                             const Line<ctype2, worldDim>& line)
{ return Vector<PromotedType<ctype1, ctype2>, worldDim>(p, line.projection(p)).length(); }

/*!
 * \brief Returns the euclidian distance between a line and a point.
 */
template<class ctype1, class ctype2, int worldDim>
PromotedType<ctype1, ctype2> computeDistance(const Line<ctype1, worldDim>& line,
                                             const Point<ctype2, worldDim>& p)
{ return computeDistance(p, line); }

/*!
 * \brief Returns the euclidian distance between a point and a segmennt.
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
 * \brief Returns the euclidian distance between a segment and a point.
 */
template<class ctype1, class ctype2, int worldDim>
PromotedType<ctype1, ctype2> computeDistance(const Segment<ctype1, worldDim>& seg,
                                             const Point<ctype2, worldDim>& p)
{ return computeDistance(p, seg); }

} // end namespace Frackit

#endif // FRACKIT_DISTANCE_HH
