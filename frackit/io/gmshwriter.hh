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
 * \ingroup IO
 * \brief Class that writes entity networks into Gmsh .geo file format.
 *        Effectively, this will write a .brep file and a .geo file in
 *        which the .brep file is included. Then, physical definitions
 *        are added to all entity and sub-domain fragments.
 */
#ifndef FRACKIT_GMSH_WRITER_HH
#define FRACKIT_GMSH_WRITER_HH

#include <string>
#include <fstream>

#include "brepwriter.hh"

namespace Frackit {

/*!
 * \ingroup IO
 * \brief Class that writes entity networks into Gmsh .geo file format.
 *        Effectively, this will write a .brep file and a .geo file in
 *        which the .brep file is included. Then, physical definitions
 *        are added to all entity and sub-domain fragments.
 */
class GmshWriter : public BRepWriter
{

public:

    /*!
     * \brief Constructor from an entity network.
     */
    template<class EntityNetwork>
    GmshWriter(const EntityNetwork& network)
    : BRepWriter(network)
    , entityDimension_(network.entityDimension())
    {}

    /*!
     * \brief Write the .geo file.
     * \param fileName The body of the file name to be used.
     * \param meshSizeAtEntities Characteristic mesh size at entities.
     * \param meshSizeAtBoundary Characteristic mesh size at domain boundaries.
     * \note Two files are created: _fileName_.brep and _fileName_.geo
     */
    template<class ctype>
    void write(const std::string& fileName,
               ctype meshSizeAtEntities,
               ctype meshSizeAtBoundary)
    {
        BRepWriter::write(fileName);

        std::ofstream geoFile(fileName + ".geo");
        geoFile << "Merge \"";
        geoFile << fileName + ".brep";
        geoFile << "\";\n";

        if (!this->domainToFragmentsMap().empty())
        {
            geoFile << "\n// Physical domain definitions\n";
            auto it = this->domainToFragmentsMap().begin();
            for (; it != this->domainToFragmentsMap().end(); ++it)
            {
                geoFile << "Physical ";
                if (entityDimension_ == 2) geoFile << "Volume(";
                else if (entityDimension_ == 1) geoFile << "Surface(";
                geoFile << it->first;
                geoFile << ") = {";
                geoFile << it->second[0];
                for (unsigned int j = 1; j < it->second.size(); ++j)
                    geoFile << ", " << it->second[j];
                geoFile << "};\n";

                // mesh size specification
                geoFile << "Characteristic Length{ PointsOf{";
                if (entityDimension_ == 2) geoFile << "Volume{";
                else if (entityDimension_ == 1) geoFile << "Surface{";
                geoFile << it->second[0];
                for (unsigned int j = 1; j < it->second.size(); ++j)
                    geoFile << ", " << it->second[j];
                geoFile << "};} } = ";
                geoFile << meshSizeAtBoundary << ";\n";
            }
        }

        if (!this->entityToFragmentsMap().empty())
        {
            geoFile << "\n// Physical entity definitions\n";
            auto it = this->entityToFragmentsMap().begin();
            for (; it != this->entityToFragmentsMap().end(); ++it)
            {
                geoFile << "Physical ";
                if (entityDimension_ == 2) geoFile << "Surface(";
                else if (entityDimension_ == 1) geoFile << "Line(";
                geoFile << it->first;
                geoFile << ") = {";
                geoFile << it->second[0];
                for (unsigned int j = 1; j < it->second.size(); ++j)
                    geoFile << ", " << it->second[j];
                geoFile << "};\n";

                // mesh size specification
                geoFile << "Characteristic Length{ PointsOf{";
                if (entityDimension_ == 2) geoFile << "Surface{";
                else if (entityDimension_ == 1) geoFile << "Line{";
                geoFile << it->second[0];
                for (unsigned int j = 1; j < it->second.size(); ++j)
                    geoFile << ", " << it->second[j];
                geoFile << "};} } = ";
                geoFile << meshSizeAtEntities << ";\n";
            }
        }

        if (!this->entityIntersectionsToFragmentsMap().empty())
        {
            geoFile << "\n// Physical entity intersections definitions\n";
            auto it = this->entityIntersectionsToFragmentsMap().begin();
            for (; it != this->entityIntersectionsToFragmentsMap().end(); ++it)
            {
                geoFile << "Physical ";
                if (entityDimension_ == 2) geoFile << "Line(";
                else if (entityDimension_ == 1) geoFile << "Point(";
                geoFile << it->first;
                geoFile << ") = {";
                geoFile << it->second[0];
                for (unsigned int j = 1; j < it->second.size(); ++j)
                    geoFile << ", " << it->second[j];
                geoFile << "};\n";
            }
        }
    }

    /*!
     * \brief Write the .geo file.
     * \param fileName The body of the file name to be used.
     * \param meshSizeAtEntities Characteristic mesh size at entities.
     * \note This uses the provided mesh size on both entities and boundaries.
     * \note Two files are created: _fileName_.brep and _fileName_.geo.
     */
    template<class ctype>
    void write(const std::string& fileName, ctype meshSizeAtEntities)
    { write(fileName, meshSizeAtEntities, meshSizeAtEntities); }

private:
    int entityDimension_;
};

} // end namespace Frackit

#endif // FRACKIT_GMSH_WRITER_HH
