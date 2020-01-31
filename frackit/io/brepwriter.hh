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
 * \ingroup IO
 * \brief Class that writes entity networks into the .brep file format.
 */
#ifndef FRACKIT_BREP_WRITER_HH
#define FRACKIT_BREP_WRITER_HH

#include <string>
#include <stdexcept>
#include <unordered_map>

#include <BRepTools.hxx>
#include <BRep_Builder.hxx>

#include <TopTools_IndexedMapOfShape.hxx>
#include <TopTools_DataMapOfShapeInteger.hxx>
#include <TopTools_ListOfShape.hxx>
#include <Standard_TypeDef.hxx>

#include <TopExp_Explorer.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS.hxx>

#include <frackit/occ/breputilities.hh>
#include <frackit/entitynetwork/entitynetwork.hh>
#include <frackit/entitynetwork/containedentitynetwork.hh>

namespace Frackit {

/*!
 * \ingroup IO
 * \brief Class that writes entity networks into the .brep file format.
 *        This creates a single TopoDS_Compound shape in which each sub-shape
 *        is uniquely defined.
 */
class BRepWriter
{
public:

    /*!
     * \brief Construction from an entity network.
     */
    BRepWriter(const EntityNetwork& network)
    {
        if (network.entityDimension() != 2)
            throw std::runtime_error("BRepWriter only implemented for 2d networks so far");

        makeEntityMap_(network);
        findEntityIntersectionShapes_();
        makeSubShapeMaps_(network);
        makeCompound_();
    }

    /*!
     * \brief Construction from a contained entity network.
     */
    BRepWriter(const ContainedEntityNetwork& network)
    {
        if (network.entityDimension() != 2)
            throw std::runtime_error("BRepWriter only implemented for 2d networks so far");
        if (network.domainDimension() != network.entityDimension() + 1)
            throw std::runtime_error("BRepWriter expects entityDim = domainDim - 1");

        makeEntityMap_(network);
        findEntityIntersectionShapes_();
        makeSubShapeMaps_(network);
        makeCompound_();
    }

    /*!
     * \brief Write the entity network to disk.
     */
    void write(const std::string& fileName)
    {
        const std::string brepFile = fileName + ".brep";
        BRepTools::Write(compound_, brepFile.c_str());
    }

protected:
    /*!
     * \brief Returns the map that maps a primary entity index to the list
     *        of fragment indices that were inserted to the compound for
     *        that primary entity index. The indices refer to the indices
     *        within the compound.
     */
    const std::unordered_map<std::size_t, std::vector<std::size_t>>& entityToFragmentsMap() const
    { return entityToFragmentsMap_; }

    /*!
     * \brief Returns the map that maps a primary entity intersectionindex to the
     *        list of intersection fragments that were inserted to the compound for
     *        that intersection. The indices refer to the indices within the compound.
     */
    const std::unordered_map<std::size_t, std::vector<std::size_t>>& entityIntersectionsToFragmentsMap() const
    { return entityIntersectionsToFragmentsMap_; }

    /*!
     * \brief Returns the map that maps to each primary sub-domain index the list
     *        of domain fragment indices that were inserted to the compound for
     *        that sub-domain. The indices refer to the indices within the compound.
     */
    const std::unordered_map<std::size_t, std::vector<std::size_t>>& domainToFragmentsMap() const
    { return domainToFragmentsMap_; }

private:
    /*!
     * \brief Copies all entities of a network into a single map.
     * \note This overload is for contained entity networks.
     */
    void makeEntityMap_(const ContainedEntityNetwork& network)
    {
        std::size_t subNetworkOffset = 0;
        for (const auto& id : network.subDomainIds())
        {
            // the map maps to the entity index within the sub-domain network
            // Therefore, we add the offset in order to ensure that we don't add
            // entities from another network to the same primary entity index
            const auto& map = network.subDomainEntityFragmentsIndexMap(id);
            for (TopTools_DataMapIteratorOfDataMapOfShapeInteger it(map); it.More(); it.Next())
            {
                allEntities_.Append(it.Key());
                allEntitiesIndexMap_.Bind(it.Key(), it.Value() + subNetworkOffset);
            }

            // the next sub-domain
            subNetworkOffset += map.Extent();
        }
    }

    /*!
     * \brief Copies all entities of a network into a single map.
     * \note This overload is for uncontained entity networks.
     */
    void makeEntityMap_(const EntityNetwork& network)
    {
        const auto& map = network.entityFragmentsIndexMap();
        for (TopTools_DataMapIteratorOfDataMapOfShapeInteger it(map); it.More(); it.Next())
        {
            allEntities_.Append(it.Key());
            allEntitiesIndexMap_.Bind(it.Key(), it.Value());
        }
    }

    /*!
     * \brief Determins the shapes that describe intersections of entities.
     */
    void findEntityIntersectionShapes_()
    {
        // TODO: Start second iterator from ++it, somehow. A first implementation
        //       of that led to some intersection edges not being meshed at the end
        for (TopTools_ListIteratorOfListOfShape it(allEntities_); it.More(); it.Next())
        {
            for (TopTools_ListIteratorOfListOfShape it2(allEntities_); it2.More(); it2.Next())
            {
                if (it.Value().IsSame(it2.Value()))
                    continue;

                const auto edges = OCCUtilities::getOverlapEdges(TopoDS::Face(it.Value()),
                                                                 TopoDS::Face(it2.Value()));

                std::vector<std::size_t> handled;
                for (unsigned int eIdx = 0; eIdx < edges.size(); ++eIdx)
                    for (const auto& isList : entityIntersections_)
                        if (isList.Contains(edges[eIdx]))
                            handled.push_back(eIdx);

                if (handled.size() != edges.size())
                {
                    entityIntersections_.resize(entityIntersections_.size()+1);
                    for (unsigned int eIdx = 0; eIdx < edges.size(); ++eIdx)
                        if (!std::count(handled.begin(), handled.end(), eIdx))
                            entityIntersections_.back().Append(edges[eIdx]);
                }
            }
        }
    }

    /*!
     * \brief Fill the sub-shape maps with all shapes that describe the network.
     * \note This overload is for contained entity networks.
     */
    void makeSubShapeMaps_(const ContainedEntityNetwork& network)
    {
        for (auto id : network.subDomainIds())
            for (const auto& fragment : network.subDomainFragments(id))
                addSolids_(fragment, id.get());

        for (auto id : network.subDomainIds())
            for (const auto& fragment : network.subDomainEntityFragments(id))
                addFaces_(fragment);
    }

    /*!
     * \brief Fill the sub-shape maps with all shapes that describe the network.
     * \note This overload is for uncontained entity networks.
     */
    void makeSubShapeMaps_(const EntityNetwork& network)
    {
        for (const auto& fragment : network.entityFragments())
            addFaces_(fragment);
    }

    /*!
     * \brief Construct the single compound describing the network.
     */
    void makeCompound_()
    {
        BRep_Builder b;
        b.MakeCompound(compound_);

        // use Standard_Integer (as OCC uses this type) instead of
        // std::size_t here to avoid compiler warning
        for(Standard_Integer i = 1; i <= vmap_.Extent(); i++) b.Add(compound_, vmap_(i));
        for(Standard_Integer i = 1; i <= emap_.Extent(); i++) b.Add(compound_, emap_(i));
        for(Standard_Integer i = 1; i <= wmap_.Extent(); i++) b.Add(compound_, wmap_(i));
        for(Standard_Integer i = 1; i <= fmap_.Extent(); i++) b.Add(compound_, fmap_(i));
        for(Standard_Integer i = 1; i <= shmap_.Extent(); i++) b.Add(compound_, shmap_(i));
        for(Standard_Integer i = 1; i <= somap_.Extent(); i++) b.Add(compound_, somap_(i));
    }

    /*!
     * \brief Add solids to the sub-shape maps
     * \param shape The shape of which the solids are to be extracted
     * \param subDomainIdx The index of the sub-domain these solids are embedded in
     */
    void addSolids_(const TopoDS_Shape& shape, std::size_t subDomainIdx)
    {
        for(TopExp_Explorer solidExp(shape, TopAbs_SOLID); solidExp.More(); solidExp.Next())
        {
            TopoDS_Solid solid = TopoDS::Solid(solidExp.Current());
            if(somap_.FindIndex(solid) < 1) // solid not yet in map
            {
                const auto curEntityIdx = somap_.Add(solid);
                domainToFragmentsMap_[subDomainIdx].push_back(curEntityIdx);

                for(TopExp_Explorer shellExp(solid, TopAbs_SHELL); shellExp.More(); shellExp.Next())
                {
                    TopoDS_Shell shell = TopoDS::Shell(shellExp.Current());
                    if(shmap_.FindIndex(shell) < 1)
                    {
                        shmap_.Add(shell);
                        addFaces_(shellExp.Current());
                    }
                }
            }
        }
    }

    /*!
     * \brief Add faces to the sub-shape maps
     * \param shape The shape of which the faces are to be extracted
     * \todo TODO: For support of 1d networks (in 2d space) we would have
     *             to add an overload taking the sub-domain index which is used
     *             in that case and which would not search for the entity index.
     */
    void addFaces_(const TopoDS_Shape& shape)
    {
        for(TopExp_Explorer faceExp(shape, TopAbs_FACE); faceExp.More(); faceExp.Next())
        {
            TopoDS_Face face = TopoDS::Face(faceExp.Current());
            const auto faceIdx = fmap_.FindIndex(face);

            if (faceIdx < 1)
            {
                const auto curEntityIdx = fmap_.Add(face);

                int entityIdx;
                if (allEntitiesIndexMap_.Find(faceExp.Current(), entityIdx))
                    entityToFragmentsMap_[entityIdx].push_back(curEntityIdx);

                for (TopExp_Explorer wireExp(face.Oriented(TopAbs_FORWARD), TopAbs_WIRE); wireExp.More(); wireExp.Next())
                {
                    TopoDS_Wire wire = TopoDS::Wire(wireExp.Current());
                    if(wmap_.FindIndex(wire) < 1)
                    {
                        wmap_.Add(wire);
                        addEdges_(wireExp.Current());
                    }
                }
            }
        }
    }

    /*!
     * \brief Add edges to the sub-shape maps
     * \param shape The shape of which the edges are to be extracted
     * \todo TODO: For support of 1d networks (in 2d space) we would have
     *             to distinguish here and search for the entity index instead
     *             of the intersection index.
     */
    void addEdges_(const TopoDS_Shape& shape)
    {
        for (TopExp_Explorer edgeExp(shape, TopAbs_EDGE); edgeExp.More(); edgeExp.Next())
        {
            TopoDS_Edge edge = TopoDS::Edge(edgeExp.Current());
            const auto edgeIdx = emap_.FindIndex(edge);

            if(edgeIdx < 1)
            {
                const auto curEntityIdx = emap_.Add(edge);

                // start from index 1 in the map
                for (std::size_t isIdx = 0; isIdx < entityIntersections_.size(); ++isIdx)
                    if (entityIntersections_[isIdx].Contains(edge))
                        entityIntersectionsToFragmentsMap_[isIdx+1].push_back(curEntityIdx);

                addVertices_(edgeExp.Current());
            }
        }
    }

    /*!
     * \brief Add vertices to the sub-shape maps
     * \param shape The shape of which the vertices are to be extracted
     * \todo TODO: For support of 1d networks (in 2d space) we would have
     *            to  check if the vertices describe entity intersections.
     */
    void addVertices_(const TopoDS_Shape& shape)
    {
        for (TopExp_Explorer vertExp(shape, TopAbs_VERTEX); vertExp.More(); vertExp.Next())
        {
            TopoDS_Vertex vertex = TopoDS::Vertex(vertExp.Current());
            if(vmap_.FindIndex(vertex) < 1)
                vmap_.Add(vertex);
        }
    }

    // containers to store the primary entities
    TopTools_ListOfShape allEntities_;
    TopTools_DataMapOfShapeInteger allEntitiesIndexMap_;
    std::vector<TopTools_ListOfShape> entityIntersections_;

    // sub-shape maps
    TopTools_IndexedMapOfShape vmap_, emap_, wmap_, fmap_, shmap_, somap_;

    // single compound describing the overall network
    TopoDS_Compound compound_;

    // index maps from primary to fragment indices
    std::unordered_map<std::size_t, std::vector<std::size_t>> entityToFragmentsMap_;
    std::unordered_map<std::size_t, std::vector<std::size_t>> entityIntersectionsToFragmentsMap_;
    std::unordered_map<std::size_t, std::vector<std::size_t>> domainToFragmentsMap_;
};

} // end namespace Frackit

#endif // FRACKIT_BREP_WRITER_HH
