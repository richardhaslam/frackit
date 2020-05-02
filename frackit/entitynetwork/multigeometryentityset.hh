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
 * \ingroup EntityNetwork
 * \brief Class in which sets of entities of different
 *        geometry types can be stored. Each entity set
 *        receives a unique identifier, which is deduced
 *        automatically while entities are added.
 */
#ifndef FRACKIT_MULTI_GEOMETRY_ENTITY_SET_HH
#define FRACKIT_MULTI_GEOMETRY_ENTITY_SET_HH

#include <tuple>
#include <vector>
#include <string>
#include <type_traits>
#include <unordered_map>

#include <frackit/geometryutilities/name.hh>
#include <frackit/geometryutilities/assign.hh>
#include <frackit/common/typetraits.hh>
#include <frackit/common/id.hh>

namespace Frackit {

/*!
 * \ingroup EntityNetwork
 * \brief Class in which sets of entities of
 *        different geometry types can be stored.
 */
template<class... Geometries>
class MultiGeometryEntitySet
{
    template<std::size_t i>
    using GeometryType = typename std::tuple_element<i, std::tuple<Geometries...>>::type;

    template<class G>
    using GeometryVectorIdPair = std::pair< std::vector<G>, Id >;

    template<class G>
    using GeometryVectorIdPairVector = std::vector< GeometryVectorIdPair<G> >;

public:

    /*!
     * \brief Adds an entity to the set with the given id.
     * \param entity The entity to be added.
     * \param id The id of the entity set.
     */
    template<class Entity>
    void addEntity(const Entity& entity, const Id& id)
    {
        static_assert(Contains<Entity, Geometries...>::value,
                      "The provided entity is of a type that was not "
                      "defined for this MultiGeometryEntitySet instance");

        // get the vector in which we store entity sets of this type
        auto& entityIdVectors = std::get< GeometryVectorIdPairVector<Entity> >(entitySets_);

        // there exists an entity set for this id already
        auto it = entitySetGeometries_.find(id.get());
        if (it != entitySetGeometries_.end())
        {
            if (it->second != geometryName(entity))
            {
                std::string msg = "MultiGeometryEntitySet::addEntity(): ";
                msg += "The type of the provided entity does not match the ";
                msg += "geometry type of the entity set with the given id!";
                throw std::runtime_error(msg);
            }

            // find the corresponding set and add entity
            for (auto& pair : entityIdVectors)
                if (pair.second == id)
                { pair.first.push_back(entity); return; }

            throw std::runtime_error("MultiGeometryEntitySet::addEntity(): Could not find entity set!");
        }

        // define a new entity set
        entityIdVectors.push_back( GeometryVectorIdPair<Entity>({}, id.get()) );
        entityIdVectors.back().first.push_back(entity);
        entitySetGeometries_[id.get()] = geometryName(entity);
    }

    /*!
     * \brief Adds an entity to the set with the given id.
     * \param entity Pointer to the abstract geometry base class
     * \param id The id of the entity set.
     */
    void addEntity(std::shared_ptr<Geometry> entity, const Id& id)
    {
        // lambda to parse pointer into an entity and add to corresponding set
        bool wasAdded = false;
        auto parse = [&] (auto& geomVectorIdPairVector)
        {
            if (wasAdded) return;

            using GeomVecIdPairVector = std::decay_t<decltype(geomVectorIdPairVector)>;
            using PairType = typename GeomVecIdPairVector::value_type;
            using GeomVectorType = typename PairType::first_type;
            using GeomType = typename GeomVectorType::value_type;

            GeomType geometry;
            if (geometryName(geometry) == entity->name())
            {
                assign(entity, geometry);
                addEntity(geometry, id);
                wasAdded = true;
            }
        };

        std::apply([&](auto& ...x){ (..., parse(x)); }, entitySets_);

        if (!wasAdded)
            throw std::runtime_error("Could not add entity from abstract base class");
    }

    /*!
     * \brief Applies the provided function to the
     *        entity set with the provided id.
     * \param id The id of the entity set on which the function should be applied.
     * \param applyFunc The function which is to be applied.
     * \return true, if the function was applied succesfully, false otherwise.
     */
    template<class ApplyFunc>
    bool applyOnSet(const Id& id, ApplyFunc&& applyFunc) const
    {
        // find the entity set for this id
        auto it = entitySetGeometries_.find(id.get());
        if (it == entitySetGeometries_.end())
            return false;

        // lambda to be executed on each vector of entity sets
        bool wasApplied = false;
        auto findAndApply = [&] (const auto& geomVectorIdPairVector)
        {
            if (wasApplied) return;

            for (const auto& vectorIdPair : geomVectorIdPairVector)
                if (vectorIdPair.second == id)
                { applyFunc(vectorIdPair.first); wasApplied = true; break; }
        };

        // find the id and apply function
        std::apply([&](auto& ...x){(..., findAndApply(x));}, entitySets_);

        // return true if the function was applied succesfully
        return wasApplied;
    }

    /*!
     * \brief Returns the number of entities contained
     *        in the set with the provided id.
     */
    std::size_t numEntities(const Id& id) const
    {
        std::size_t size = 0;
        auto getSize = [&] (const auto& entitySet)
        { size = entitySet.size(); };

        if (applyOnSet(id, getSize))
            return size;

        return 0;
    }

    /*!
     * \brief Returns the overall number of entities.
     */
    std::size_t numEntities() const
    {
        std::size_t numEnt = 0;
        for (const auto& mapEntry : entitySetGeometries_)
            numEnt += numEntities( Id(mapEntry.first) );
        return numEnt;
    }

    /*!
     * \brief Export all entity sets into the provided
     *        entity network builder.
     * \param builder An instance of a builder class.
     * \param subDomainId The id of the sub-domain for which
     *                    the entities should be defined.
     */
    template<class NetworkBuilder>
    void exportEntitySets(NetworkBuilder& builder, const Id& subDomainId) const
    {
        // lambda to add sets of a specific geometry type
        auto addSets = [&] (const auto& geomVectorIdPairVector)
        {
            for (const auto& vectorIdPair : geomVectorIdPairVector)
                builder.addSubDomainEntities(vectorIdPair.first, subDomainId);
        };

        // find the id and apply function
        std::apply([&](auto& ...x){(..., addSets(x));}, entitySets_);
    }

    /*!
     * \brief Export all entity sets into the provided
     *        entity network builder.
     * \param builder An instance of a builder class.
     * \note This overload works with builders that create
     *       uncontained network, i.e. only the network is built
     *       independent of the embedding domain.
     */
    template<class NetworkBuilder>
    void exportEntitySets(NetworkBuilder& builder) const
    {
        // lambda to add sets of a specific geometry type
        auto addSets = [&] (const auto& geomVectorIdPairVector)
        {
            for (const auto& vectorIdPair : geomVectorIdPairVector)
                builder.addEntities(vectorIdPair.first);
        };

        // find the id and apply function
        std::apply([&](auto& ...x){(..., addSets(x));}, entitySets_);
    }

    /*!
     * \brief Export various entity sets into the
     *        provided entity network builder.
     * \param idList The list containing the ids of the entity sets to be exported.
     * \param builder An instance of a builder class.
     * \param subDomainId The id of the sub-domain for which
     *                    the entities should be defined.
     */
    template<class NetworkBuilder>
    void exportEntitySets(const std::initializer_list<Id>& idList,
                          NetworkBuilder& builder,
                          const Id& subDomainId) const
    {
        std::vector<Id> ids(idList);

        // lambda to add sets of a specific geometry type
        auto addSets = [&] (const auto& geomVectorIdPairVector)
        {
            for (const auto& vectorIdPair : geomVectorIdPairVector)
                if (std::count(ids.begin(), ids.end(), vectorIdPair.second))
                    builder.addSubDomainEntities(vectorIdPair.first, subDomainId);
        };

        // find the id and apply function
        std::apply([&](auto& ...x){(..., addSets(x));}, entitySets_);
    }

private:
    // We store vectors of entities for each geometry types.
    std::tuple< GeometryVectorIdPairVector<Geometries>... > entitySets_;

    // Memorize to each added entity set the geometry type used for it
    std::unordered_map<std::size_t, std::string> entitySetGeometries_;
};

} // end namespace Frackit

#endif // FRACKIT_MULTI_GEOMETRY_ENTITY_SET_HH
