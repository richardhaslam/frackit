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
#ifndef FRACKIT_PYTHON_ENTITYNETWORK_BUILDER_WRAPPERS_HH
#define FRACKIT_PYTHON_ENTITYNETWORK_BUILDER_WRAPPERS_HH

#include <frackit/python/geometry/brepwrapper.hh>
#include <frackit/entitynetwork/networkbuilder.hh>
#include <frackit/common/id.hh>

namespace Frackit::Python {

namespace Detail {

    // Wrappers around the builder classes to support brep wrappers
    // we use the convenience function getUnwrappedShape() to support all types
    template<class ctype>
    class EntityNetworkBuilderWrapper
    : public EntityNetworkBuilder<ctype>
    {
        using ParentType = EntityNetworkBuilder<ctype>;

    public:
        template<class Entity>
        void addSubDomainEntity(const Entity& entity, Id subDomainId)
        { ParentType::addSubDomainEntity(getUnwrappedShape(entity), subDomainId); }

        template<class Entity>
        void addEntity(const Entity& entity)
        { addSubDomainEntity(entity, Id(0)); }

        template<class Domain>
        void addConfiningSubDomain(const Domain& domain, Id subDomainId)
        { ParentType::addConfiningSubDomain(getUnwrappedShape(domain), subDomainId); }

        template<class Domain>
        void addSubDomain(const Domain& domain, Id subDomainId)
        { ParentType::addSubDomain(getUnwrappedShape(domain), subDomainId); }
    };

    template<class ctype>
    class ContainedEntityNetworkBuilderWrapper
    : public ContainedEntityNetworkBuilder<ctype>
    {
        using ParentType = ContainedEntityNetworkBuilder<ctype>;

    public:
        template<class Entity>
        void addSubDomainEntity(const Entity& entity, Id subDomainId)
        { ParentType::addSubDomainEntity(getUnwrappedShape(entity), subDomainId); }

        template<class Domain>
        void addConfiningSubDomain(const Domain& domain, Id subDomainId)
        { ParentType::addConfiningSubDomain(getUnwrappedShape(domain), subDomainId); }

        template<class Domain>
        void addSubDomain(const Domain& domain, Id subDomainId)
        { ParentType::addSubDomain(getUnwrappedShape(domain), subDomainId); }
    };

} // end namespace Detail
} // end namespace Frackit::Python

#endif
