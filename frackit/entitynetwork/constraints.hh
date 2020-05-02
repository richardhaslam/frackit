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
 * \ingroup Constraints
 * \brief Class that enforces geometric
 *        constraints on networks of entities.
 */
#ifndef FRACKIT_ENTITYNETWORK_CONSTRAINTS_HH
#define FRACKIT_ENTITYNETWORK_CONSTRAINTS_HH

#include <stdexcept>
#include <vector>
#include <type_traits>

#include <frackit/distance/distance.hh>
#include <frackit/magnitude/magnitude.hh>

#include <frackit/geometry/geometry.hh>
#include <frackit/geometryutilities/applyongeometry.hh>

#include <frackit/intersection/intersect.hh>
#include <frackit/intersection/emptyintersection.hh>

// default engines to be used for intersection angles, etc.
#include <frackit/intersection/intersectionangle.hh>

#include "constraints_impl/admissibledimension.hh"
#include "constraints_impl/admissiblemagnitude.hh"
#include "constraints_impl/distancetoboundary.hh"

namespace Frackit {

/*!
 * \brief Forward declaration of the constraints class.
 */
template<class ST = double,
         class AE = IntersectionAngle<ST> >
class EntityNetworkConstraints;

/*!
 * \ingroup Constraints
 * \brief Convenience function to construct an instance of
 *        the constraints class using the default template arguments.
 */
template<class ST = double>
EntityNetworkConstraints<ST> makeDefaultConstraints()
{ return EntityNetworkConstraints<ST>(IntersectionAngle<ST>()); }

/*!
 * \ingroup Constraints
 * \brief Class which defines and checks constraints on
 *        the geometric relationships between entities of
 *        a network. As constraints one can define:
 *        - minimum distance between entities
 *        - minimum intersection angle
 *        - minimum intersection magnitude
 *        - minimum distance between the intersection and entity boundaries
 * \note For the evaluation of the distance between the boundary of entity
 *       intersections and the entity boundaries, only those boundary parts
 *       of the intersection geometry are considered that to not intersect
 *       themselves with the boundary of an entity, as this would return
 *       zero distance. However, this case is admissible, only small features
 *       related to the non-intersecting boundaries are of interest here.
 * \tparam ST The type used for the scalar constraint values
 * \tparam AE The engine used for angle computations between intersecting geometries.
 *            This engine is required  to return the angle between geometries when
 *            the operator() is called with the two geometries and the intersection
 *            result as arguments.
 */
template<class ST, class AE>
class EntityNetworkConstraints
{
public:
    //! Export type used for constraints
    using Scalar = ST;

    //! Exort the engine used for angle computations
    using AngleComputationEngine = AE;

    /*!
     * \brief Default constructor.
     * \note This is only available if the engines are default constructible
     */
    EntityNetworkConstraints()
    : EntityNetworkConstraints(AngleComputationEngine())
    {
        static_assert(std::is_default_constructible<AE>::value,
                      "Angle computation engine not default constructible. "
                      "Use constructor taking the engine as argument instead");
    }

    /*!
     * \brief The constructor.
     * \param angleEngine An instance of the engine used for angle computations.
     */
    EntityNetworkConstraints(const AngleComputationEngine& angleEngine)
    : angleEngine_(angleEngine)
    , minDistance_(),            minIsAngle_()
    , minIsMagnitude_(),         minIsDistance_()
    , useMinDistance_(false),    useMinIsAngle_(false)
    , useMinIsMagnitude_(false), useMinIsDistance_(false)
    , intersectionEps_(),        useIntersectionEps_(false)
    {
        // per default, we do not allow equi-dimensional intersections
        allowEquiDimIS_ = false;
    }

    //! Set the constraint for the minimum distance between entities
    void setMinDistance(Scalar minDistance)
    {
        minDistance_ = minDistance;
        useMinDistance_ = true;
    }

    //! Set the constraint for the minimum intersection angle
    void setMinIntersectingAngle(Scalar minIsAngle)
    {
        minIsAngle_ = minIsAngle;
        useMinIsAngle_ = true;
    }

    //! Set the constraint for the minimum intersection magnitude
    void setMinIntersectionMagnitude(Scalar minIsMagnitude)
    {
        minIsMagnitude_ = minIsMagnitude;
        useMinIsMagnitude_ = true;
    }

    //! Set constraint for the min distance of an intersection to entity boundaries
    void setMinIntersectionDistance(Scalar minIsDistance)
    {
        minIsDistance_ = minIsDistance;
        useMinIsDistance_ = true;
    }

    //! Deactivate the distance constraint
    void neglectMinDistance() { useMinDistance_ = false; }
    //! Deactivate the intersecting angle constraint
    void neglectMinIntersectionAngle() { useMinIsAngle_ = false; }
    //! Deactivate the intersection magnitude constraint
    void neglectMinIntersectionMagnitude() { useMinIsMagnitude_ = false; }
    //! Deactivate the intersecting distance constraint
    void neglectMinIntersectionDistance() { useMinIsDistance_ = false; }

    //! Set an epsilon value to be used in the intersection computation
    void setIntersectionEpsilon(Scalar eps)
    {
        intersectionEps_ = eps;
        useIntersectionEps_ = true;
    }

    //! Reactivate the use of the default epsilon
    void setDefaultIntersectionEpsilon()
    { useIntersectionEps_ = false; }

    /*!
     * \brief Define if intersections that have the same dimension
     *        as the entites are admissible. For example, two-dimensional
     *        surface intersections as intersections between disks.
     * \param value the flag to be set (true/false)
     */
    void allowEquiDimensionalIntersections(bool value)
    { allowEquiDimIS_ = value; }

    /*!
     * \brief Check if a pair of geometries fulfills the constraints
     * \param geo1 The first geometry
     * \param geo2 The second geometry
     * \returns True if all constraints are fulfilled, false otherwise
     */
    template<class Geo1, class Geo2>
    bool evaluate(const Geo1& geo1, const Geo2& geo2) const
    {
        // if the provided geometries are pointers (we support shared_ptr here)
        // to the abstract base class, we have to first cast them into the derived types
        static constexpr bool isBasePtr1 = std::is_same_v<Geo1, std::shared_ptr<Geometry>>;
        static constexpr bool isBasePtr2 = std::is_same_v<Geo2, std::shared_ptr<Geometry>>;

        if constexpr (isBasePtr1 && !isBasePtr2)
        {
            auto evalCast1 = [&] (const auto& actualGeo1) { return evaluate_(actualGeo1, geo2); };
            return applyOnGeometry(evalCast1, geo1);
        }
        else if constexpr (!isBasePtr1 && isBasePtr2)
        {
            auto evalCast2 = [&] (const auto& actualGeo2) { return evaluate_(geo1, actualGeo2); };
            return applyOnGeometry(evalCast2, geo2);
        }
        else if constexpr (isBasePtr1 && isBasePtr2)
        {
            auto evalCast1 = [&] (const auto& actualGeo1)
            {
                auto evalCast2 = [&] (const auto& actualGeo2) { return evaluate_(actualGeo1, actualGeo2); };
                return applyOnGeometry(evalCast2, geo2);
            };
            return applyOnGeometry(evalCast1, geo1);
        }
        // both geometries are actual geometry types
        else
            return evaluate_(geo1, geo2);
    }

    /*!
     * \brief Check if an entity fulfills the constraints for
     *        all entities of the provided entity set.
     * \tparam Geo1 The geometry type of the entities in the set
     * \tparam Geo2 The geometry type of the entity to be checked
     * \param entitySet An entity set
     * \param entity The geometry of the entity to be checked
     * \returns True if all constraints are fulfilled, false otherwise
     */
    template<class Geo1, class Geo2>
    bool evaluate(const std::vector<Geo1>& entitySet, const Geo2& entity) const
    {
        return std::all_of(entitySet.begin(),
                           entitySet.end(),
                           [&] (const auto& e) { return evaluate(e, entity); });
    }

    /*!
     * \brief Check if an entity fulfills the constraints for
     *        all entities of the provided entity set.
     * \tparam Geo1 The geometry type of the entity to be checked
     * \tparam Geo2 The geometry type of the entities in the set
     * \param entity The geometry of the entity to be checked
     * \param entitySet An entity set
     * \returns True if all constraints are fulfilled, false otherwise
     */
    template<class Geo1, class Geo2>
    bool evaluate(const Geo1& entity, const std::vector<Geo2>& entitySet) const
    { return evaluate(entitySet, entity); }

private:
    /*!
     * \brief Check if a pair of geometries fulfills the constraints
     * \param geo1 The first geometry
     * \param geo2 The second geometry
     * \returns True if all constraints are fulfilled, false otherwise
     */
    template<class Geo1, class Geo2>
    bool evaluate_(const Geo1& geo1, const Geo2& geo2) const
    {
        static constexpr int dim = DimensionalityTraits<Geo1>::geometryDimension();
        static_assert(dim == DimensionalityTraits<Geo2>::geometryDimension(),
                      "We expect entities to be of same dimension");

        const bool checkIs = useMinIsMagnitude_ || useMinIsAngle_ || useMinIsDistance_;
        if (!checkIs && !useMinDistance_)
            return true;

        const auto isection = !useIntersectionEps_ ? intersect(geo1, geo2)
                                                   : intersect(geo1, geo2, intersectionEps_);

        if ( !isEmptyIntersection(isection) )
        {
            // check  if dimensionality constraint is violated
            if (!allowEquiDimIS_ && !ConstraintImpl::isAdmissibleDimension(isection, dim-1))
                return false;

            // magnitude constraint
            if ( useMinIsMagnitude_ && !ConstraintImpl::isAdmissibleMagnitude(isection, minIsMagnitude_) )
                return false;

            // angle constraint
            if ( useMinIsAngle_ && angleEngine_(geo1, geo2, isection) < minIsAngle_ )
                return false;

            // constraint on distance of intersection to geometry boundaries
            if ( !useMinIsDistance_ )
                return true;

            using namespace ConstraintImpl;
            if (!isAdmissibleDistanceToBoundary(isection, geo1, minIsDistance_))
                return false;
            return isAdmissibleDistanceToBoundary(isection, geo2, minIsDistance_);
        }

        // no intersection - check distance constraint
        else if (useMinDistance_)
            return computeDistance(geo1, geo2) >= minDistance_;

        return true;
    }

    AngleComputationEngine angleEngine_; //! Algorithms to compute angles between intersecting geometries

    Scalar minDistance_;    //!< Minimum distance allowed between two entities
    Scalar minIsAngle_;     //!< Minimum angle in which two entities are allowed to intersect
    Scalar minIsMagnitude_; //!< Minimum magnitude of the intersection geometry
    Scalar minIsDistance_;  //!< Minimum distance between the intersection and the two entity boundaries

    //! flags which constraints are active
    bool useMinDistance_;
    bool useMinIsAngle_;
    bool useMinIsMagnitude_;
    bool useMinIsDistance_;
    bool allowEquiDimIS_;

    Scalar intersectionEps_;  //! Tolerance value to be used for intersections
    bool useIntersectionEps_; //! Stores wether or not a user-defined epsilon value was set
};

} // end namespace Frackit

#endif // FRACKIT_ENTITYNETWORK_CONSTRAINTS_HH
