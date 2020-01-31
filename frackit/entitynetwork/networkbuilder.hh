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
 * \brief Contains builder classes for entity networks.
 */
#ifndef FRACKIT_ENTITY_NETWORK_BUILDER_HH
#define FRACKIT_ENTITY_NETWORK_BUILDER_HH

#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <memory>

#include <TopTools_DataMapOfShapeInteger.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopoDS_Shape.hxx>

#include <frackit/common/id.hh>
#include <frackit/common/extractdimension.hh>

#include <frackit/precision/defaultepsilon.hh>
#include <frackit/precision/precision.hh>
#include <frackit/occ/breputilities.hh>

#include "entitynetwork.hh"
#include "containedentitynetwork.hh"

namespace Frackit {

// Forward declaration of the builder class for entity networks
template<class ctype = double>
class EntityNetworkBuilder;

// Forward declaration of the builder class for contained entity networks
template<class ctype = double>
class ContainedEntityNetworkBuilder;

/*!
 * \ingroup EntityNetwork
 * \brief Base class for builders of entity networks.
 *        Stores data related to the entities and sub-domains.
 */
template<class ctype = double>
class EntityNetworkBuilderBase
{
public:

    /*!
     * \brief Define a tolerance value to be used
     *        for boolean operations on the entity shapes.
     * \param eps Tolerance value to be used
     */
    void setEpsilon(ctype eps)
    {
        eps_ = eps;
        epsilonIsSet_ = true;
    }

    /*!
     * \brief Computes and sets a tolerance value using the extents
     *        of ths entities inserted into the network so far.
     */
    void setDefaultEpsilon()
    {
        // if no entities have been set, set default precision
        if (subDomains_.empty() && networks_.empty())
            eps_ = Precision<ctype>::confusion();
        else
        {
            // Use the minimum of all default epsilons of the sub-domains or entities
            eps_ = std::numeric_limits<ctype>::max();
            using std::min;

            if (!subDomains_.empty())
                for (const auto& sdPair : subDomains_)
                    eps_ = min(eps_, defaultEpsilon(sdPair.second));
            else if (!networks_.empty())
                for (const auto& netPair : networks_)
                    for (const auto& entity : netPair.second)
                        eps_ = min(eps_, defaultEpsilon(entity));
        }

        epsilonIsSet_ = true;
    }

    /*!
     * \brief Returns true if the network has been
     *        built, i.e. if a call to build() has ocurred.
     */
    bool isBuilt() const
    { return built_; }

    /*!
     * \brief Defines a sub-domain which potentially
     *        contains an embedded entity network.
     *        Networks defined for this sub-domain
     *        will not be confined by the sub-domain's
     *        boundary.
     * \param domain The sub-domain geometry
     * \param subDomainId The identifier (index) to be used for this sub-domain
     */
    template<class Domain>
    void addSubDomain(const Domain& domain, Id subDomainId)
    {
        const auto domainDim = getDimension(domain);
        if (!networks_.empty() && domainDim <= entityDimension_)
            throw std::runtime_error("Sub-domain dimension must be greater than entity dimension");
        if (!subDomains_.empty() && domainDim != domainDimension_)
            throw std::runtime_error("Sub-domain dimension does not match to previously added");
        if (subDomains_.find(subDomainId.get()) != subDomains_.end())
            throw std::runtime_error("Sub-domain id already taken!");

        domainDimension_ = domainDim;
        subDomains_[subDomainId.get()] = OCCUtilities::getShape(domain);
        subDomainBoundsNetwork_[subDomainId.get()] = false;
    }

    /*!
     * \brief Defines a sub-domain which potentially
     *        contains an embedded entity network.
     *        Using this interface, the sub-domain is
     *        defined as confining, that is, the
     *        network embedded in it will be confined
     *        by the sub-domain's boundary.
     * \param domain The sub-domain geometry
     * \param subDomainId The identifier (index) to be used for this sub-domain
     */
    template<class Domain>
    void addConfiningSubDomain(const Domain& domain, Id subDomainId)
    {
        addSubDomain(domain, subDomainId);
        subDomainBoundsNetwork_[subDomainId.get()] = true;
    }

    /*!
     * \brief Adds an entity to the network embedded in
     *        the sub-domain with index subDomainIndex.
     * \param entity The entity geometry
     * \param subDomainId The identifier (index) to be used for this sub-domain
     */
    template<class Entity>
    void addSubDomainEntity(const Entity& entity, Id subDomainId)
    {
        const auto entityDim = getDimension(entity);
        if (!networks_.empty() && getDimension(entity) != entityDimension_)
            throw std::runtime_error("Entity dimension does not match to previously added");
        if (!subDomains_.empty() && getDimension(entity) >= domainDimension_)
            throw std::runtime_error("Entity dimension must be smaller than domain dimension");

        entityDimension_ = entityDim;
        networks_[subDomainId.get()].Append(OCCUtilities::getShape(entity));
    }

    /*!
     * \brief Adds a set of entities to the network
     *        embedded in sub-domain with index subDomainIndex.
     * \param entities The set of entities
     * \param subDomainId The identifier (index) to be used for this sub-domain
     */
    template<class EntityNetwork>
    void addSubDomainEntities(const EntityNetwork& entities, Id subDomainId)
    {
        for (const auto& entity : entities)
            addSubDomainEntity(entity, subDomainId);
    }

    /*!
     * \brief Deletes all data.
     */
    void clear()
    {
        clearBuild_();
        subDomains_.clear();
        subDomainBoundsNetwork_.clear();
        networks_.clear();
    }

protected:

    /*!
     * \brief Confines the networks to the
     *        (sub-)domains they are embedded in.
     */
    void confineNetworks_()
    {
        for (const auto& indexNetworkPair : networks_)
        {
            const auto subDomainIndex = indexNetworkPair.first;
            const auto& network = indexNetworkPair.second;

            // initialize index map with network entities
            std::size_t idx = 1;
            auto& fragmentList = networkFragmentLists_[subDomainIndex];
            auto& fragmentMap = networkFragmentMaps_[subDomainIndex];
            for (const auto& entityShape : network)
            {
                fragmentList.Append(entityShape);
                fragmentMap.Bind(entityShape, idx++);
            }

            // confine network only if sub-domain was set
            if (subDomains_.find(subDomainIndex) == subDomains_.end())
                continue;

            // confine network with sub-domain or entire domain
            TopTools_ListOfShape domainShapeList;
            if (subDomainBoundsNetwork_.at(subDomainIndex))
                domainShapeList.Append(subDomains_.at(subDomainIndex));
            else
                domainShapeList.Append(getSubDomainUnion_());

            BRepAlgoAPI_Common intersection;
            intersection.SetArguments(network);
            intersection.SetTools(domainShapeList);
            intersection.SetRunParallel(false);
            intersection.SetFuzzyValue(getEpsilon_());
            intersection.Build();
            if(!intersection.IsDone())
                throw std::runtime_error("Could not confine network");

            // keep track of original entity indices
            TopTools_ListOfShape newList;
            TopTools_DataMapOfShapeInteger newIndexMap;
            for (const auto& entityShape : network)
            {
                int entityIdx;
                if (!fragmentMap.Find(entityShape, entityIdx))
                    throw std::runtime_error("Could not find entity index");

                // map all fragments to this new index
                const bool isDeleted = intersection.IsDeleted(entityShape);
                const auto& generated = intersection.Generated(entityShape);
                const auto& modified = intersection.Modified(entityShape);
                if (!isDeleted && generated.IsEmpty() && modified.IsEmpty())
                {
                    newIndexMap.Bind(entityShape, entityIdx);
                    newList.Append(entityShape);
                }
                for (const auto& fragment : generated)
                {
                    if (!newIndexMap.IsBound(fragment)) newIndexMap.Bind(fragment, entityIdx);
                    if (!newList.Contains(fragment)) newList.Append(fragment);
                }
                for (const auto& fragment : modified)
                {
                    if (!newIndexMap.IsBound(fragment)) newIndexMap.Bind(fragment, entityIdx);
                    if (!newList.Contains(fragment)) newList.Append(fragment);
                }
            }

            // update the fragment map
            fragmentList = std::move(newList);
            fragmentMap = std::move(newIndexMap);
        }
    }

    /*!
     * \brief Computes the complete fragmentation
     *        of the network and/or the sub-domains.
     * \param fragmentNetworkOnly If set to true, the domains will not be
     *                            fragmented. This can be used if one is only
     *                            interested in the entity network and not in the domains.
     */
    void fragment_(bool fragmentNetworkOnly)
    {
        // List with all entities of the network
        TopTools_ListOfShape allEntities;
        TopTools_DataMapOfShapeInteger entityToSubDomainIndex;
        for (const auto& fragmentListPair : networkFragmentLists_)
        {
            for (const auto& fragment : fragmentListPair.second)
            {
                allEntities.Append(fragment);
                entityToSubDomainIndex.Bind(fragment, fragmentListPair.first);
            }
        }

        // All sub-domains (omit this if only network is considered)
        TopTools_ListOfShape allDomains;
        TopTools_DataMapOfShapeInteger domainToSubDomainIndex;
        if (!fragmentNetworkOnly)
        {
            for (const auto& subDomainDataPair : subDomains_)
            {
                allDomains.Append(subDomainDataPair.second);
                domainToSubDomainIndex.Bind(subDomainDataPair.second,
                                            subDomainDataPair.first);
            }
        }

        // List of both entity and sub-domain shapes
        TopTools_ListOfShape allShapes;
        for (const auto& entityShape : allEntities)
            allShapes.Append(entityShape);
        for (const auto& domainShape : allDomains)
            allShapes.Append(domainShape);

        // fragment all shapes
        BRepAlgoAPI_BuilderAlgo fragmentAlgo;
        fragmentAlgo.SetRunParallel(false);
        fragmentAlgo.SetArguments(allShapes);
        fragmentAlgo.SetFuzzyValue(getEpsilon_());
        fragmentAlgo.Build();
        if(!fragmentAlgo.IsDone())
            throw std::runtime_error("Could not perform fragmentation");

        // map the resulting entities back to the original entity index
        std::unordered_map<std::size_t, TopTools_ListOfShape> newFragmentLists;
        std::unordered_map<std::size_t, TopTools_DataMapOfShapeInteger> newFragmentMaps;
        for (const auto& shape : allEntities)
        {
            // get subdomain index
            int subDomainIdx;
            if (!entityToSubDomainIndex.Find(shape, subDomainIdx))
                throw std::runtime_error("Could not find sub-domain index");

            // get original entity index
            int entityIdx;
            const auto& map = networkFragmentMaps_.at(subDomainIdx);
            if (!map.Find(shape, entityIdx))
                throw std::runtime_error("Could not find original entity index");

            // map all fragments to this new index
            const bool isDeleted = fragmentAlgo.IsDeleted(shape);
            const auto& generated = fragmentAlgo.Generated(shape);
            const auto& modified = fragmentAlgo.Modified(shape);
            auto& newList = newFragmentLists[subDomainIdx];
            auto& newMap = newFragmentMaps[subDomainIdx];
            if (!isDeleted && generated.IsEmpty() && modified.IsEmpty())
            {
                newMap.Bind(shape, entityIdx);
                newList.Append(shape);
            }
            for (const auto& fragment : generated)
            {
                if (!newMap.IsBound(fragment)) newMap.Bind(fragment, entityIdx);
                if (!newList.Contains(fragment)) newList.Append(fragment);
            }
            for (const auto& fragment : modified)
            {
                if (!newMap.IsBound(fragment)) newMap.Bind(fragment, entityIdx);
                if (!newList.Contains(fragment)) newList.Append(fragment);
            }
        }

        // overwrite fragment maps with the newly obtained
        for (auto& dataPair : newFragmentLists)
            networkFragmentLists_[dataPair.first] = std::move(dataPair.second);
        for (auto& dataPair : newFragmentMaps)
            networkFragmentMaps_[dataPair.first] = std::move(dataPair.second);

        // map the resulting domain entities back to the original entity index
        for (const auto& domainShape : allDomains)
        {
            // get subdomain index
            int subDomainIdx;
            if (!domainToSubDomainIndex.Find(domainShape, subDomainIdx))
                throw std::runtime_error("Could not find sub-domain index");

            // store the fragments of this sub-domain
            const bool isDeleted = fragmentAlgo.IsDeleted(domainShape);
            const auto& generated = fragmentAlgo.Generated(domainShape);
            const auto& modified = fragmentAlgo.Modified(domainShape);

            auto& sdFragmentList = subDomainFragmentLists_[subDomainIdx];
            if (!isDeleted && generated.IsEmpty() && modified.IsEmpty())
                sdFragmentList.Append(domainShape);
            for (const auto& fragment : generated)
                if (!sdFragmentList.Contains(fragment))
                    sdFragmentList.Append(fragment);
            for (const auto& fragment : modified)
                if (!sdFragmentList.Contains(fragment))
                    sdFragmentList.Append(fragment);
        }
    }

    /*!
     * \brief Returns the tolerance value.
     */
    ctype getEpsilon_()
    {
        if (!epsilonIsSet_)
            setDefaultEpsilon();
        return eps_;
    }

    /*!
     * \brief Returns the shape of the union of all sub-domains.
     */
    const TopoDS_Shape& getSubDomainUnion_()
    {

        if (subDomainUnion_.IsNull())
        {
            if (subDomains_.size() == 0)
                throw std::runtime_error("No sub-domains defined yet");

            if (subDomains_.size() == 1)
                subDomainUnion_ = subDomains_.begin()->second;
            else
            {
                // compute union of all subdomains
                std::vector<TopoDS_Shape> subDomains;
                for (const auto& sdPair  : subDomains_)
                    subDomains.push_back(sdPair.second);
                subDomainUnion_ = OCCUtilities::fuse(subDomains, getEpsilon_());
            }
        }

        return subDomainUnion_;
    }

    /*!
     * \brief Initializes the fragment lists prior to a build.
     */
    void initializeBuild_()
    {
        // make sure the fragment maps exist for all subdomains
        // even if there are no networks defined for some of them
        for (const auto& sdDataPair : subDomains_)
        {
            networkFragmentLists_[sdDataPair.first] = TopTools_ListOfShape();
            networkFragmentMaps_[sdDataPair.first] = TopTools_DataMapOfShapeInteger();
        }
    }

    /*!
     * \brief Finishes a building process.
     */
    void finishBuild_()
    {
        built_ = true;
    }

    /*!
     * \brief Clears data generated during the build process
     *        such that a new call to build() can be made.
     */
    void clearBuild_()
    {
        built_ = false;
        subDomainFragmentLists_.clear();
        networkFragmentLists_.clear();
        networkFragmentMaps_.clear();
        subDomainUnion_.Nullify();
    }

    // User-defined sub-domains and entity networks
    std::unordered_map<std::size_t, TopoDS_Shape> subDomains_;
    std::unordered_map<std::size_t, bool> subDomainBoundsNetwork_;
    std::unordered_map<std::size_t, TopTools_ListOfShape> networks_;

    // dimensionalities
    int domainDimension_;
    int entityDimension_;

    // The shape of the union of the user-defined sub-domains
    TopoDS_Shape subDomainUnion_;

    // Epsilon value to be used for fragmentation
    ctype eps_;

    // state variables
    bool epsilonIsSet_ = false;
    bool built_ = false;

    // store the fragments of each subdomain after fragmentation
    std::unordered_map<std::size_t, TopTools_ListOfShape> subDomainFragmentLists_;

    // fragment maps, map each fragment to the index of the original entity
    std::unordered_map<std::size_t, TopTools_ListOfShape> networkFragmentLists_;
    std::unordered_map<std::size_t, TopTools_DataMapOfShapeInteger> networkFragmentMaps_;
};

/*!
 * \ingroup EntityNetwork
 * \brief Builder class for entity networks,
 *        which neglect the sub-domains in which
 *        the network is embedded in and only carries
 *        information about the network itself.
 */
template<class ctype>
class EntityNetworkBuilder
: public EntityNetworkBuilderBase<ctype>
{

public:
    //! Export network type to be built
    using EntityNetwork = Frackit::EntityNetwork;

    /*!
     * \brief Adds an entity to the network.
     *        This function can be used when no sub-domain
     *        specifications are made. Internally, the entity
     *        is added to the sub-domain with index 0.
     */
    template<class Entity>
    void addEntity(const Entity& entity)
    {
        const auto entityDim = getDimension(entity);
        if (!this->networks_.empty() && getDimension(entity) != this->entityDimension_)
            throw std::runtime_error("Entity dimension does not match to previously added");
        if (!this->subDomains_.empty() && getDimension(entity) >= this->domainDimension_)
            throw std::runtime_error("Entity dimension must be smaller than domain dimension");

        this->entityDimension_ = entityDim;
        this->networks_[0].Append(OCCUtilities::getShape(entity));
    }

    /*!
     * \brief Adds a set of entities to the network.
     */
    template<class EntityNetwork>
    void addEntities(const EntityNetwork& entities)
    {
        for (const auto& entity : entities)
            addEntity(entity);
    }

    /*!
     * \brief Builds the network.
     */
    const EntityNetwork& build()
    {
        // prepare build
        this->clearBuild_();
        this->initializeBuild_();

        // first confine each network to its domain
        this->confineNetworks_();

        // Then fragment the network entities
        this->fragment_(/*fragmentNetworkOnly*/true);

        // merge all fragments into one list
        TopTools_ListOfShape fragmentList;
        for (auto& sdFragmentDataPair : this->networkFragmentLists_)
            fragmentList.Append(sdFragmentDataPair.second);

        // merge all index maps into one
        std::size_t entityIndexOffset = 0;
        TopTools_DataMapOfShapeInteger fragmentIndexMap;
        for (auto& dataPair : this->networkFragmentMaps_)
        {
            const auto& indexMap = dataPair.second;
            for (TopTools_DataMapIteratorOfDataMapOfShapeInteger it(indexMap); it.More(); it.Next())
                fragmentIndexMap.Bind(it.Key(), it.Value() + entityIndexOffset);
            entityIndexOffset += indexMap.Extent();
        }

        // create the network
        network_ = std::make_unique<EntityNetwork>(this->entityDimension_,
                                                   std::move(fragmentList),
                                                   std::move(fragmentIndexMap));

        // set flag that network fragmentation was carried out
        this->finishBuild_();

        return *network_;
    }

    /*!
     * \brief Returns the built network.
     */
    const EntityNetwork& network() const
    {
        if (!this->isBuilt())
            throw std::runtime_error("Network has not been built yet!");
        return *network_;
    }

private:
    std::unique_ptr<EntityNetwork> network_;
};

/*!
 * \ingroup EntityNetwork
 * \brief Builder class for entity networks contained
 *        in (possibly multiple) sub-domains.
 */
template<class ctype>
class ContainedEntityNetworkBuilder
: public EntityNetworkBuilderBase<ctype>
{
public:

    //! Export network type to be built
    using EntityNetwork = Frackit::ContainedEntityNetwork;

    /*!
     * \brief Build the network.
     */
    const EntityNetwork& build()
    {
        // maybe print a warning
        if (std::any_of(this->networks_.begin(),
                        this->networks_.end(),
                        [&] (const auto& pair)
                        { return this->subDomains_.find(pair.first) == this->subDomains_.end(); }))
            std::cout << "\t --> build(): there are networks with no sub-domain defined" << std::endl;

        // prepare build
        this->clearBuild_();
        this->initializeBuild_();

        // first confine each network to its domain
        this->confineNetworks_();

        // Then fragment the complete set of shapes
        this->fragment_(/*fragmentNetworkOnly*/false);

        // create the network
        network_ = std::make_unique<EntityNetwork>(this->entityDimension_,
                                                   this->domainDimension_,
                                                   std::move(this->subDomainFragmentLists_),
                                                   std::move(this->networkFragmentLists_),
                                                   std::move(this->networkFragmentMaps_));

        // set flag that network fragmentation was carried out
        this->finishBuild_();

        return *network_;
    }

    /*!
     * \brief Return the built network.
     */
    const EntityNetwork& network() const
    {
        if (!this->isBuilt())
            throw std::runtime_error("Network has not been built yet!");
        return *network_;
    }

private:
    std::unique_ptr<EntityNetwork> network_;
};

} // end namespace Frackit

#endif // FRACKIT_ENTITY_NETWORK_BUILDER_HH
