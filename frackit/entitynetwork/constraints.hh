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
 * \brief Class that enforces geometric
 *        constraints on networks of entities.
 */
#ifndef FRACKIT_ENTITYNETWORK_CONSTRAINTS_HH
#define FRACKIT_ENTITYNETWORK_CONSTRAINTS_HH

#include <stdexcept>
#include <variant>
#include <vector>

#include <frackit/distance/distance.hh>
#include <frackit/distance/distancetoboundary.hh>

#include <frackit/magnitude/magnitude.hh>
#include <frackit/intersection/intersect.hh>
#include <frackit/intersection/intersectionpredicates.hh>

namespace Frackit {

/*!
 * \brief Class which defines and checks constraints on
 *        the geometric relationships between entities of
 *        a network. As constraints one can define:
 *        - minimum distance
 *        - minimum intersection angle
 *        - minimum intersection magnitude
 *        - minimum distance between the intersection and an entity boundary
 * \tparam ST The type used for the scalar constraint values
 */
template<class ST = double>
class EntityNetworkConstraints
{

public:
    //! Export type used for constraints
    using Scalar = ST;

    /*!
     * \brief Default constructor, deactivates all constraints.
     */
    EntityNetworkConstraints()
    : minDistance_(),            minIsAngle_()
    , minIsMagnitude_(),         minIsDistance_()
    , useMinDistance_(false),    useMinIsAngle_(false)
    , useMinIsMagnitude_(false), useMinIsDistance_(false)
    , intersectionEps_(),        useIntersectionEps_(false)
    {}

    /*!
     * \brief The constructor defining all constraints
     * \param minDistance Minimum distance allowed between two (non-intersecting) entities
     * \param minIsAngle Minimum angle in which two entities are allowed to intersect
     * \param minIsMagnitude Minimum magnitude of intersection allowed
     * \param minIsDistance Minimum distance of an intersection to intersecting entity boundaries
     * \note The default epsilon is set to be used for intersection computations
     */
    EntityNetworkConstraints(Scalar minDistance,
                             Scalar minIsAngle,
                             Scalar minIsMagnitude,
                             Scalar minIsDistance)
    : minDistance_(minDistance),       minIsAngle_(minIsAngle)
    , minIsMagnitude_(minIsMagnitude), minIsDistance_(minIsDistance)
    , useMinDistance_(false),          useMinIsAngle_(false)
    , useMinIsMagnitude_(false),       useMinIsDistance_(false)
    , intersectionEps_(),              useIntersectionEps_(false)
    {}

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
     * \brief Check if a pair of geometries fulfills the constraints
     * \param geo1 The first geometry
     * \param geo2 The second geometry
     * \returns True if all constraints are fulfilled, false otherwise
     */
    template<class Geo1, class Geo2>
    bool evaluate(const Geo1& geo1, const Geo2& geo2)
    {
        const bool checkIs = useMinIsMagnitude_ || useMinIsAngle_ || useMinIsDistance_;
        if (!checkIs && !useMinDistance_)
            return true;

        const auto is = !useIntersectionEps_ ? intersect(geo1, geo2)
                                             : intersect(geo1, geo2, intersectionEps_);

        if ( !IntersectionPredicates::isEmpty(is) )
        {
            // magnitude constraint
            if ( useMinIsMagnitude_ && !isAdmissibleMagnitude_(is) )
                return false;

            // angle constraint
            if ( useMinIsAngle_ && IntersectionPredicates::angle(geo1, geo2) < minIsAngle_ )
                return false;

            // constraint on distance of intersection to geometry boundaries
            if ( !useMinIsDistance_ )
                return true;

            bool isAdmissible;
            std::visit([&] (auto&& is) { isAdmissible = isAdmissibleDistanceToBoundary_(is, geo1); }, is);
            if (!isAdmissible)
                return false;

            return std::visit([&] (auto&& is) { return isAdmissibleDistanceToBoundary_(is, geo2); }, is);
        }

        // no intersection - check distance constraint
        else if (useMinDistance_)
            return computeDistance(geo1, geo2) >= minDistance_;

        return true;
    }

    /*!
     * \brief Check if an entity fulfills the constraints for
     *        all entities of the provided entity set.
     * \tparam Geo1 The geometry type of the entities in the set
     * \tparam Geo2 The geometry type of the entity to be checked
     * \param entitySet An entity network
     * \param entity The geometry of the entity to be checked
     * \returns True if all constraints are fulfilled, false otherwise
     */
    template<class Geo1, class Geo2>
    bool evaluate(const std::vector<Geo1>& entitySet, const Geo2& entity)
    {
        return std::all_of(entitySet.begin(),
                           entitySet.end(),
                           [&] (const auto& e) { return evaluate(e, entity); });
    }

private:
    Scalar minDistance_;    //!< Minimum distance allowed between two entities
    Scalar minIsAngle_;     //!< Minimum angle in which two entities are allowed to intersect
    Scalar minIsMagnitude_; //!< Minimum magnitude of the intersection geometry
    Scalar minIsDistance_;  //!< Minimum distance between the intersection and the two entity boundaries

    //! flags which constraints are active
    bool useMinDistance_;
    bool useMinIsAngle_;
    bool useMinIsMagnitude_;
    bool useMinIsDistance_;

    Scalar intersectionEps_;  //! Tolerance value to be used for intersections
    bool useIntersectionEps_; //! Stores wether or not a user-defined epsilon value was set

    //! \todo TODO Doc me.
    template<class T>
    bool isAdmissibleMagnitude_(const std::vector<T>& intersection)
    {
        return std::all_of(intersection.begin(),
                           intersection.end(),
                           [&] (const auto& is) { return isAdmissibleMagnitude_(is); });
    }

    //! \todo TODO Doc me.
    template<class... T>
    bool isAdmissibleMagnitude_(const std::variant<T...>& intersection)
    { return std::visit([&] (auto&& is) { return isAdmissibleMagnitude_(is); }, intersection); }

    //! \todo TODO Doc me.
    template<class Geo>
    bool isAdmissibleMagnitude_(const Geo& isGeometry) const
    {
        // zero dimensional (point) intersections fulfill any magnitude constraint
        if (Geo::myDimension() == 0) return true;
        // empty intersections fulfill any magnitude constraint
        if (IsEmptyIntersection<Geo>::value) return true;
        // any other intersection geometry
        return computeMagnitude(isGeometry) >= minIsMagnitude_;
    }

    //! \todo TODO Doc me.
    template<class T, class Geo>
    bool isAdmissibleDistanceToBoundary_(const std::vector<T>& intersection, const Geo& geo)
    {
        return std::all_of(intersection.begin(),
                           intersection.end(),
                           [&] (const auto& is) { return isAdmissibleDistanceToBoundary_(is, geo); });
    }

    //! \todo TODO Doc me.
    template<class... T, class Geo>
    bool isAdmissibleDistanceToBoundary_(const std::variant<T...>& intersection, const Geo& geo)
    { return std::visit([&] (auto&& is) { return isAdmissibleDistanceToBoundary_(is, geo); }, intersection); }

    //! \todo TODO Doc me.
    template<class IsGeom, class Geo>
    bool isAdmissibleDistanceToBoundary_(const IsGeom& isGeometry, const Geo& geometry) const
    {
        // empty intersections fulfill the distance constraint
        if (IsEmptyIntersection<Geo>::value)
            return true;

        return computeDistanceToBoundary(isGeometry, geometry) >= minIsDistance_;
    }
};

} // end namespace Frackit

#endif // FRACKIT_ENTITYNETWORK_CONSTRAINTS_HH
