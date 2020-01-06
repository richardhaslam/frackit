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
 * \brief \todo TODO doc me.
 */
#ifndef FRACKIT_ELLIPSE_ARC_HH
#define FRACKIT_ELLIPSE_ARC_HH

#include <cmath>
#include <stdexcept>

#include <gp_Vec.hxx>
#include <gp_Dir.hxx>
#include <gp_Elips.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Ellipse.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <Standard_Handle.hxx>

#include "precision.hh"
#include "ellipse.hh"
#include "vector.hh"

namespace Frackit {

/*!
 * \brief \todo TODO doc me.
 */
template<class CT, int worldDim>
class EllipseArc;

/*!
 * \brief \todo TODO doc me.
 */
template<class CT>
class EllipseArc<CT, /*worldDim=*/3>
: public Ellipse<CT, /*worldDim=*/3>
{
    using ParentType = Frackit::Ellipse<CT, /*worldDim=*/3>;

public:
    //! export dimensionality
    static constexpr int myDimension() { return 1; }
    static constexpr int worldDimension() { return 3; }

    //! export type used for coordinates
    using typename ParentType::ctype;

    //! export underlying geometry types
    using typename ParentType::Point;
    using Ellipse = ParentType;

    /*!
     * \brief \todo TODO doc me.
     * \note if source and target points are the same
     *       (within the precision), we create a full
     *       ellipse instead of zero-length arc.
     */
    EllipseArc(const Ellipse& ellipse,
               const Point& source,
               const Point& target)
    : ParentType(ellipse)
    , source_(source)
    , target_(target)
    {
        const auto eps = ellipse.minorAxisLength()*Precision<ctype>::confusion();
        assert(ParentType::contains(source, eps));
        assert(ParentType::contains(target, eps));

        if (source.isEqual(target, eps))
        {
            isFullEllipse_ = true;
            sourceAngle_ = this->getAngle(source);
            if (sourceAngle_ < Precision<ctype>::angular()) targetAngle_ = 2.0*M_PI;
            else targetAngle_ = sourceAngle_;
        }
        else
        {
            isFullEllipse_ = false;
            sourceAngle_ = this->getAngle(source);
            targetAngle_ = this->getAngle(target);
        }
    }

    //! \todo TODO doc me.
    static std::string name() { return "EllipseArc"; }

    //! \todo TODO doc me.
    const Point& source() const { return source_; }
    //! \todo TODO doc me.
    const Point& target() const { return target_; }

    //! Returns true if the arc describes a full ellipse
    bool isFullEllipse() const { return isFullEllipse_; }

    //! \todo TODO doc me.
    ctype length() const
    {
        if (isFullEllipse())
             return M_PI*(this->majorAxisLength() + this->minorAxisLength());

        // create Geom_Ellipse and trim to describe arc
        gp_Dir normal(this->normal().x(), this->normal().y(), this->normal().z());
        gp_Dir majorAx(this->majorAxis().x(), this->majorAxis().y(), this->majorAxis().z());
        gp_Pnt center(this->center().x(), this->center().y(), this->center().z());
        gp_Elips gpEllipse(gp_Ax2(center, normal, majorAx),
                           this->majorAxisLength(),
                           this->minorAxisLength());
        Handle(Geom_Curve) ellipseHandle = new Geom_Ellipse(gpEllipse);
        Handle(Geom_Curve) arcHandle = targetAngle_ < sourceAngle_
                                       ? new Geom_TrimmedCurve(ellipseHandle, sourceAngle_, targetAngle_)
                                       : new Geom_TrimmedCurve(ellipseHandle, sourceAngle_, targetAngle_ + 2.0*M_PI);

        const auto uMin = arcHandle->FirstParameter();
        const auto uMax = arcHandle->LastParameter();
        GeomAdaptor_Curve adaptorCurve(arcHandle, uMin, uMax);
        return GCPnts_AbscissaPoint::Length(adaptorCurve, uMin, uMax);
    }

    //! Return the ellipse that supports this arc
    Ellipse supportingEllipse() const
    {
        return Ellipse(this->center(),
                       this->majorAxis(),
                       this->minorAxis(),
                       this->majorAxisLength(),
                       this->minorAxisLength());
    }

    //! Returns true if a point is on the arc
    //! \todo note about choice of eps
    bool contains(const Point& p, ctype eps, bool checkIfOnEllipse = true) const
    {
        if (checkIfOnEllipse)
            if (!ParentType::contains(p, eps))
                return false;

        const auto angle = this->getAngle(p, checkIfOnEllipse);
        if (sourceAngle_ > targetAngle_)
            return angle >= sourceAngle_ || angle <= targetAngle_;
        else
            return angle >= sourceAngle_ && angle <= targetAngle_;
    }

    //! Returns true if a point is on the arc
    //! \todo note about choice of eps
    bool contains(const Point& p, bool checkIfOnEllipse = true) const
    {
        using std::min;
        auto eps = min(this->majorAxisLength(), this->minorAxisLength());
        eps *= Precision<ctype>::confusion();

        return contains(p, eps, checkIfOnEllipse);
    }

    //! Returns the point on the arc for the given parameter
    //! \note It has to be 0.0 <= param <= 1.0, where 0.0
    //!       corresponds to the source and 1.0 to the target.
    Point getPoint(ctype param) const
    {
        assert(param >= 0.0 && param <= 1.0);

        ctype deltaAngle;
        if (isFullEllipse_) deltaAngle = 2.0*M_PI;
        if (sourceAngle_ < targetAngle_) deltaAngle = targetAngle_ - sourceAngle_;
        else deltaAngle = 2.0*M_PI - sourceAngle_ + targetAngle_;

        return ParentType::getPointFromAngle(sourceAngle_ + deltaAngle*param);
    }

private:
    Point source_;
    Point target_;

    ctype sourceAngle_;
    ctype targetAngle_;
    bool isFullEllipse_;
};

} // end namespace Frackit

#endif // FRACKIT_ELLIPSE_ARC_HH
