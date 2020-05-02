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
 * \ingroup Sampling
 * \brief Sampler which can sample multiple geometry types.
 *        For each geometry type, various sampler classes can be defined,
 *        each of which is associated with a unique id. The sampling strategy
 *        defines the id of the sampler for which to obtain the next geometry.
 */
#ifndef FRACKIT_MULTI_GEOMETRY_SAMPLER_HH
#define FRACKIT_MULTI_GEOMETRY_SAMPLER_HH

#include <tuple>
#include <vector>
#include <string>
#include <memory>
#include <type_traits>

#include <frackit/geometry/geometry.hh>
#include <frackit/common/typetraits.hh>
#include <frackit/common/id.hh>

#include "geometrysampler.hh"
#include "samplingstrategy.hh"
#include "sequentialsamplingstrategy.hh"

namespace Frackit {

/*!
 * \ingroup Sampling
 * \brief Sampler which can sample multiple geometry types.
 * \note Per default, we use a sequential sampling strategy here.
 */
template<class... Geometries>
class MultiGeometrySampler
{
    template<std::size_t i>
    using GeometryType = typename std::tuple_element<i, std::tuple<Geometries...>>::type;

    template<class G>
    using SamplerPtrIdPair = std::pair< std::shared_ptr< GeometrySampler<G> >, Id >;

    template<class G>
    using SamplerPtrIdPairVector = std::vector< SamplerPtrIdPair<G> >;

    template<class G>
    std::shared_ptr<Geometry> sample(std::shared_ptr<GeometrySampler<G>> samplerPtr)
    { return std::make_shared<G>( (*samplerPtr)() ); }

public:

    /*!
     * \brief Default constructor. Per default, we sample sequentially.
     */
    MultiGeometrySampler()
    : samplingStrategy_(std::make_shared<SequentialSamplingStrategy>())
    {}

    /*!
     * \brief Constructor taking a sampling strategy.
     */
    template<class Strategy>
    MultiGeometrySampler(const Strategy& strategy)
    {
        setSamplingStrategy(strategy);
    }

    /*!
     * \brief Set a sampling strategy to be used.
     */
    template<class Strategy>
    void setSamplingStrategy(const Strategy& strategy)
    {
        static_assert(std::is_base_of_v<SamplingStrategy, Strategy>,
                      "The provided sampling strategy does not inherit from the strategy interface");
        samplingStrategy_ = std::make_shared<Strategy>(strategy);
    }

    /*!
     * \brief Defines a new geometry sampler.
     * \param sampler The sampler class
     * \param id The id with which this sampler is to be associated
     */
    template<class Sampler>
    void addGeometrySampler(const Sampler& sampler, const Id& id)
    {
        using SampledGeometry = typename Sampler::Geometry;
        static_assert(Contains<SampledGeometry, Geometries...>::value,
                      "The provided sampler does not sample any of the "
                      "geometries that were defined for this multi sampler instance");

        // check if id was already inserted
        if (samplingStrategy_->hasId(id))
            throw std::runtime_error("addGeometrySampler(): Id already taken!");

        // get the vector in which we store samplers of this type
        auto& samplerPtrPairVector = std::get< SamplerPtrIdPairVector<SampledGeometry> >(samplers_);

        // create new shared pointer on a copy of the given sampler object
        using CurrentPtr = std::shared_ptr<GeometrySampler<SampledGeometry>>;
        CurrentPtr currentPtr = std::make_shared<Sampler>(sampler);

        // add the provided sampler to it
        samplerPtrPairVector.push_back( std::make_pair(currentPtr, id) );

        // pass the id to the strategy
        samplingStrategy_->addId(id);
    }

    /*!
     * \brief Sample a geometry.
     * \param id An instance of Id, in which the id of
     *           the sampler is stored that has been used
     *           to generate the returned geometry.
     */
    std::shared_ptr<Geometry> operator() (Id& id)
    {
        if (!samplingStrategy_)
            throw std::runtime_error("MultiGeometrySampler::operator(): no sampling strategy set");

        id = samplingStrategy_->getNextId();

        // shared pointer to the geometry to be sampled
        std::shared_ptr<Geometry> result(nullptr);

        // lambda to be executed on each tuple entry
        auto doSampling = [&] (auto& samplerPtrIdPairVec)
        {
            // skip this if the sampler was found already
            if (result) return;

            // otherwise search this sampler vector for the id
            for (auto& pair : samplerPtrIdPairVec)
                if (pair.second == id)
                { assert(!result); result = sample(pair.first); break; }
        };

        // find the id and sample geometry
        std::apply([&](auto& ...x){(..., doSampling(x));}, samplers_);

        // error handling
        if (!result)
        {
            std::string msg = "MultiGeometrySampler::operator(): ";
            msg += "no sampler has been added for this id.\n";
            msg += "In case you use a custom sampling strategy, ";
            msg += "please verify that the strategy you provided ";
            msg += "contains the same ids as the MultiGeometrySampler\n";
            throw std::runtime_error(msg);
        }

        return result;
    }

private:
    // We store a vector of pointers to samplers for each geometry type.
    // Thus, for each geometry type, an arbitrary number of samplers can be added.
    // Moreover, to each sampler a unique ID is assigned which the user defines upon add.
    std::tuple< SamplerPtrIdPairVector<Geometries>... > samplers_;

    // Store a pointer on the sampling strategy used
    std::shared_ptr<SamplingStrategy> samplingStrategy_;
};

} // end namespace Frackit

#endif // FRACKIT_MULTI_GEOMETRY_SAMPLER_HH
