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
 * \brief Class to randomly generate convex polygons in n-dimensional space.
 */
#ifndef FRACKIT_POLYGON_SAMPLER_HH
#define FRACKIT_POLYGON_SAMPLER_HH

#include <array>
#include <cmath>
#include <memory>
#include <limits>
#include <random>
#include <algorithm>
#include <stdexcept>
#include <type_traits>

#include <frackit/common/math.hh>
#include <frackit/precision/precision.hh>

#include <frackit/geometry/point.hh>
#include <frackit/geometry/vector.hh>
#include <frackit/geometry/direction.hh>
#include <frackit/geometry/polygon.hh>

#include "pointsampler.hh"
#include "geometrysampler.hh"

namespace Frackit {

//! Forward declaration of the default traits class.
template<class ctype, int worldDim>
struct DefaultPolygonSamplerTraits;

/*!
 * \ingroup Sampling
 * \brief Default traits class for sampling polygons in 3d space.
 */
template<class ctype>
struct DefaultPolygonSamplerTraits<ctype, 3>
{
    // Orientation of the polygon plane within the x-y plane
    using StrikeAngleDistribution = std::normal_distribution<ctype>;

    // Dip angle of the polygon plane measured against the x-y plane
    using DipAngleDistribution = std::normal_distribution<ctype>;

    // Distributions for the sizes (lengths) of the polygons in the two directions
    using StrikeLengthDistribution = std::uniform_real_distribution<ctype>;
    using DipLengthDistribution = std::uniform_real_distribution<ctype>;

    // distribution for the number of corners of the polygons
    using NumCornersDistribution = std::uniform_int_distribution<int>;
};

//! Forward declaration of the dimension-specific sampler class for polygons
template< int worldDim,
          class ctype = double,
          class T = DefaultPolygonSamplerTraits<ctype, worldDim> >
class PolygonSampler;

/*!
 * \ingroup Sampling
 * \brief Sampler for convex polygons in 3d space.
 *
 *        The algorithm is a modification of the algorithm proposed in:
 *        Pavel Valtr. “Probability that n random points are in convex position.”
 *        Discrete & Computational Geometry 13.1 (1995): 637-643.
 *        http://dx.doi.org/10.1007/BF02574070
 *
 *        Instead of random combination of sampled x- and y-values to form the
 *        edges of the polygon, we combine the randomly sampled values such that
 *        the lengths of the polygon edges are maximized.
 */
template< class ctype, class T >
class PolygonSampler<3, ctype, T>
: public GeometrySampler< Polygon<ctype, 3> >
{
    using Point = Frackit::Point<ctype, 3>;
    using Vector = Frackit::Vector<ctype, 3>;
    using Direction = Frackit::Direction<ctype, 3>;

    // distribution used to sample coordinates in the unit interval
    using IntervalDistribution = std::uniform_real_distribution<ctype>;

    //! underlying distributions
    using StrikeAngleDistribution = typename T::StrikeAngleDistribution;
    using DipAngleDistribution = typename T::DipAngleDistribution;
    using StrikeLengthDistribution = typename T::StrikeLengthDistribution;
    using DipLengthDistribution = typename T::DipLengthDistribution;
    using NumCornersDistribution = typename T::NumCornersDistribution;

public:
    //! export underlying geometry type
    using Polygon = Frackit::Polygon<ctype, 3>;
    using Geometry = Polygon;

    //! export traits class
    using Traits = T;

    /*!
     * \brief The constructor.
     * \param pointSampler Point sampler object to sample center points
     * \param strikeAngleDistro The distribution used for the strike angle
     * \param dipAngleDistro The distribution used for the dip angle
     * \param strikeLengthDistro The distribution used for the size in strike direction
     * \param dipLengthDistro The distribution used for the size in dip direction
     * \param numCornersDistro Distribution to sample the number of corners from
     * \note The number of corners must be >= 3 to define a polygon. Thus, users
     *       must choose distributions in numCornersDistro that are able to produce
     *       integer values above or equal to 3. Otherwise, one is trapped in the
     *       loop within the () operator. After 100 tries, an error is thrown.
     */
    template<class PointSamplerImpl>
    PolygonSampler(const PointSamplerImpl& pointSampler,
                   const StrikeAngleDistribution& strikeAngleDistro,
                   const DipAngleDistribution& dipAngleDistro,
                   const StrikeLengthDistribution& strikeLengthDistro,
                   const DipLengthDistribution& dipLengthDistro,
                   const NumCornersDistribution& numCornersDistro)
    : pointSampler_(std::make_shared<PointSamplerImpl>(pointSampler))
    , generator_(std::random_device{}())
    , p_strike_angle_(strikeAngleDistro)
    , p_dip_angle_(dipAngleDistro)
    , p_strike_length_(strikeLengthDistro)
    , p_dip_length_(dipLengthDistro)
    , p_num_corners_(numCornersDistro)
    , p_interval_(0.0, 1.0)
    {
        static_assert(std::is_base_of<PointSampler<Point>, PointSamplerImpl>::value,
                      "The provided point sampler does not inherit from the point sampler interface");
    }

    /*!
     * \brief Generate a random polygon.
     */
    Polygon operator() () override
    { return Polygon(sampleCorners()); }

    /*!
     * \brief Randomly generate a point cloud that represents
     *        the corner points of a convex polygon.
     */
    std::vector<Point> sampleCorners()
    {
        int count = 0;
        int numCorners = 0;
        while (numCorners < 3 && count < 100)
        { numCorners = p_num_corners_(generator_); count++; };

        if (count == 100)
        {
            std::string msg = "Did not obtain a value >= 3 for the number of corners after 100 tries. ";
            msg += "Please revisit your choice of distribution for that parameter.";
            throw std::runtime_error(msg);
        }

        // randomly sample within interval [0, 1]
        std::vector<ctype> xValues; xValues.reserve(numCorners);
        std::vector<ctype> yValues; yValues.reserve(numCorners);
        for (int i = 0; i < numCorners; ++i) xValues.push_back(p_interval_(generator_));
        for (int i = 0; i < numCorners; ++i) yValues.push_back(p_interval_(generator_));

        auto it = std::min_element(xValues.begin(), xValues.end()); auto xMin = *it; xValues.erase(it);
        it = std::min_element(yValues.begin(), yValues.end());      auto yMin = *it; yValues.erase(it);
        it = std::max_element(xValues.begin(), xValues.end());      auto xMax = *it; xValues.erase(it);
        it = std::max_element(yValues.begin(), yValues.end());      auto yMax = *it; yValues.erase(it);

        // make sure the interval is exploited
        if (xMin > 0.1) { xMin = 0.0; } if (yMin > 0.1) { yMin = 0.0; }
        if (xMax < 0.9) { xMax = 1.0; } if (yMax < 0.9) { yMax = 1.0; }

        // split values into two sets
        std::array<std::vector<ctype>, 2> xVals;
        std::array<std::vector<ctype>, 2> yVals;
        for (auto& v : xVals) { v.reserve(numCorners); v.push_back(xMin); }
        for (auto& v : yVals) { v.reserve(numCorners); v.push_back(yMin); }

        for (int i = 0; i < numCorners-2; ++i)
        {
            xVals[i%2].push_back(xValues[i]);
            yVals[i%2].push_back(yValues[i]);
        }

        // sort the arrays and add max values
        for (auto& v : xVals) { std::sort(v.begin(), v.end()); v.push_back(xMax); }
        for (auto& v : yVals) { std::sort(v.begin(), v.end()); v.push_back(yMax); }

        // transform values into deltas
        std::adjacent_difference(xVals[0].begin(), xVals[0].end(), xVals[0].begin());
        std::adjacent_difference(xVals[1].begin(), xVals[1].end(), xVals[1].begin());
        std::adjacent_difference(yVals[0].begin(), yVals[0].end(), yVals[0].begin());
        std::adjacent_difference(yVals[1].begin(), yVals[1].end(), yVals[1].begin());

        // flip sign of second arrays
        std::for_each(xVals[1].begin(), xVals[1].end(), [] (auto& v) { v *= -1.0; });
        std::for_each(yVals[1].begin(), yVals[1].end(), [] (auto& v) { v *= -1.0; });

        // store results in delta arrays (skip first entry, is no delta)
        std::vector<ctype> deltaX; deltaX.reserve(numCorners);
        std::vector<ctype> deltaY; deltaY.reserve(numCorners);
        for (int i : std::array<int, 2>{{0, 1}}) std::move(xVals[i].begin()+1, xVals[i].end(), std::back_inserter(deltaX));
        for (int i : std::array<int, 2>{{0, 1}}) std::move(yVals[i].begin()+1, yVals[i].end(), std::back_inserter(deltaY));
        assert(static_cast<int>(deltaX.size()) == numCorners);
        assert(static_cast<int>(deltaY.size()) == numCorners);

        // compute the maximum deltaX/deltaY values
        std::sort(deltaX.begin(), deltaX.end());
        auto maxDeltaX = deltaX.back();
        auto maxDeltaY = *std::max_element(deltaY.begin(), deltaY.end());
        const auto maxDx2 = maxDeltaX*maxDeltaX + maxDeltaY*maxDeltaY;

        // construct edge vectors by selecting pairs of deltaX
        // and deltaY values whose length is closest to maxDx2
        using std::abs;
        std::vector<std::pair<ctype, ctype>> deltas;
        deltas.reserve(numCorners);
        for (auto dx : deltaX)
        {
            const auto tmp = maxDx2 - dx*dx;
            auto diff = [tmp] (auto& dy) { return abs(tmp - dy*dy); };
            auto comp = [&diff] (auto& dy1, auto& dy2) { return diff(dy1) < diff(dy2); };

            auto it = std::min_element(deltaY.begin(), deltaY.end(), comp);
            deltas.emplace_back(std::make_pair(dx, *it));
            deltaY.erase(it);
        }

        // scale vectors with actual dimensions (make sure they are > 0.0)
        auto strikeLength = p_strike_length_(generator_);
        auto dipLength = p_dip_length_(generator_);

        int countStrike = 0;
        while (strikeLength < 0.0 && countStrike < 100)
        { strikeLength = p_strike_length_(generator_); countStrike++; }
        if (countStrike == 100)
            throw std::runtime_error("Could not sample positive strike length in 100 tries");

        int countDip = 0;
        while (dipLength < 0.0 && countDip < 100)
        { dipLength = p_dip_length_(generator_); countDip++; }
        if (countStrike == 100)
            throw std::runtime_error("Could not sample positive dip length in 100 tries");

        std::vector<Vector> edgeVectors;
        edgeVectors.reserve(numCorners);
        for (auto [dx, dy] : deltas)
            edgeVectors.emplace_back(Vector{dx*dipLength, dy*strikeLength, 0.0});

        auto getAngle = [] (const auto& v)
        {
            if (abs(v.y()) < 1e-6) return v.x() > 0 ? 0 : M_PI;
            if (abs(v.x()) < 1e-6) return v.y() > 0 ? M_PI/2.0 : 3.0*M_PI/2.0;

            using std::atan;
            const auto angle = atan(v.y()/v.x());
            if (v.x() > 0.0) return angle < 0.0 ? 2.0*M_PI + angle : angle;
            else return M_PI + angle;
        };

        // sort edges by angle
        std::sort(edgeVectors.begin(), edgeVectors.end(),
                  [&] (const auto& v1, const auto& v2) { return getAngle(v1) < getAngle(v2); });

        // rotate them by strike and dip angle
        using std::sin; using std::cos;
        const auto strikeAngle = p_strike_angle_(generator_);
        const auto strikeDir = Direction(Vector(-sin(strikeAngle), cos(strikeAngle), 0.0));
        rotate(edgeVectors, Direction(Vector(0.0, 0.0, 1.0)), strikeAngle);
        rotate(edgeVectors, strikeDir, p_dip_angle_(generator_));

        // combine them to a polygon (leave out last vertex - is same as first)
        using std::min; using std::max;
        std::vector<Point> vertices;
        vertices.reserve(numCorners);
        vertices.emplace_back();

        // first vertex inserted above is (0.0, 0.0, 0.0)
        xMin = 0.0; xMax = 0.0; yMin = 0.0; yMax = 0.0;
        ctype zMin = 0.0; ctype zMax = 0.0;
        for (unsigned int i = 0; i < edgeVectors.size()-1; ++i)
        {
            vertices.emplace_back(vertices.back() + edgeVectors[i]);
            xMin = min(xMin, vertices.back().x()); xMax = max(xMax, vertices.back().x());
            yMin = min(yMin, vertices.back().y()); yMax = max(yMax, vertices.back().y());
            zMin = min(zMin, vertices.back().z()); zMax = max(zMax, vertices.back().z());
        }

        // check if last and first vertex are the same
        assert((vertices.back() + edgeVectors.back()).isEqual(vertices.front()));

        // move to sampled center
        const auto c = (*pointSampler_)();
        const Point bboxCenter({xMin + 0.5*(xMax-xMin),
                                yMin + 0.5*(yMax-yMin),
                                zMin + 0.5*(zMax-zMin)});
        const Vector dispVec(bboxCenter, c);
        for (auto& v : vertices)
            v += dispVec;

        return vertices;
    }

 private:
    std::shared_ptr<PointSampler<Point>> pointSampler_; //!< pointer to the sampler for polygon position
    std::default_random_engine generator_;              //!< Random number generator

    StrikeAngleDistribution p_strike_angle_;   //!< Distribution used for the strike angle
    DipAngleDistribution p_dip_angle_;         //!< Distribution used for the dip angle
    StrikeLengthDistribution p_strike_length_; //!< Distribution used for the size in strike direction
    DipLengthDistribution p_dip_length_;       //!< Distribution used for the size in dip direction
    NumCornersDistribution p_num_corners_;     //!< Distribution used for the number of polygon corners
    IntervalDistribution p_interval_;          //!< Distribution used for coordinate sampling in interval [0, 1]
};

} // end namespace Frackit

#endif // FRACKIT_POLYGON_SAMPLER_HH
