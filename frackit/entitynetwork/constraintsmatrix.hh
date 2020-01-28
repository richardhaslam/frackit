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
 * \brief Class that allows defining a matrix of constraints.
 *        This is useful in contexts where multiple entity sets
 *        are defined, and where different constraints should be
 *        fulfilled for different pairs of entity sets.
 */
#ifndef FRACKIT_ENTITYNETWORK_CONSTRAINTS_MATRIX_HH
#define FRACKIT_ENTITYNETWORK_CONSTRAINTS_MATRIX_HH

#include <stdexcept>
#include <utility>
#include <vector>
#include <unordered_map>
#include <initializer_list>

#include <frackit/common/id.hh>
#include <frackit/common/idpair.hh>

namespace Frackit {

/*!
* \ingroup Constraints
 * \brief Class that allows defining a matrix of constraints.
 *        For each constraint that is defined, a list of id pairs
 *        is passed that defines between which ids (representing
 *        the ids of the corresponding entity sets) the constraints
 *        should hold.
 */
template<class C>
class EntityNetworkConstraintsMatrix
{
    using IndexSet = std::vector<std::size_t>;
    using IndexPair = std::pair<std::size_t, std::size_t>;
    using IndexPairVector = std::vector<IndexPair>;

public:
    //! export underlying constraints type
    using Constraints = C;

    /*!
     * \brief Add a constraint that should be fulfilled
     *        between the two entity sets defined in the
     *        provided id pair.
     */
    void addConstraints(const Constraints& c, const IdPair& idPair)
    {
        addConstraint_(c);
        updateIndexMaps_(idPair);
    }

    /*!
     * \brief Add a constraint that should be fulfilled
     *        between all pairs of entity sets defined
     *        in the provided list of id pairs.
     */
    void addConstraints(const Constraints& c,
                        const std::initializer_list<IdPair>& idPairs)
    {
        addConstraint_(c);
        for (const auto& idPair : idPairs)
            updateIndexMaps_(idPair);
    }

    /*!
     * \brief Evaluate if an entity fulfills all constraints.
     * \param entitySets The class containing all entity sets.
     * \param entity The entity to be checked.
     * \param id The id of the entity set to which the entity belongs.
     */
    template<class EntitySets, class Entity>
    bool evaluate(const EntitySets& entitySets,
                  const Entity& entity,
                  const Id& id)
    {
        auto it = constraintsIndexMap_.find(id.get());
        if (it == constraintsIndexMap_.end())
            throw std::runtime_error("evaluate: no constraints set for given id");

        for (const auto& idxPair : it->second)
        {
            const auto otherId = Id(idxPair.first);
            const auto& constraints = constraints_[idxPair.second];

            // lambda for the evaluation of the constraints
            bool isAdmissible;
            auto evalConstraintOnSet = [&] (const auto& entitySet) -> void
            { isAdmissible = constraints.evaluate(entitySet, entity); };

            // applyOnSet() returns false if lambda was not applied
            // This is the case e.g. when the set is empty. Here,
            // return false only if applyOnSet() was successful.
            if (entitySets.applyOnSet(otherId, evalConstraintOnSet))
                if (!isAdmissible)
                    return false;
        }

        return true;
    }

private:
    /*!
     * \brief Adds a constraints instance to the container.
     */
    void addConstraint_(const Constraints& c)
    {
        constraints_.push_back(c);
    }

    /*!
     * \brief Updates the index maps after a call to addConstraint_().
     */
    void updateIndexMaps_(const IdPair& idPair)
    {
        const auto firstIdx = idPair.first().get();
        const auto secondIdx = idPair.second().get();

        auto it = relations_.find(idPair.first().get());
        if (it != relations_.end())
            if (std::count(it->second.begin(),
                           it->second.end(),
                           idPair.second().get()))
                throw std::runtime_error("ConstraintsMatrix: Relation was already defined!");

        const auto cIdx = constraints_.size() - 1;
        constraintsIndexMap_[firstIdx].push_back(std::make_pair(secondIdx, cIdx));
        relations_[firstIdx].push_back(secondIdx);
    }

    //! stores all constraints
    std::vector<Constraints> constraints_;
    //! maps to each binary relation the index of the constraints
    std::unordered_map<std::size_t, IndexPairVector> constraintsIndexMap_;
    //! stores to each id with which other ids a relation was defined
    std::unordered_map<std::size_t, IndexSet> relations_;
};

} // end namespace Frackit

#endif // FRACKIT_ENTITYNETWORK_CONSTRAINTS_HH
