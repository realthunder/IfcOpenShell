/********************************************************************************
 *                                                                              *
 * This file is part of IfcOpenShell.                                           *
 *                                                                              *
 * IfcOpenShell is free software: you can redistribute it and/or modify         *
 * it under the terms of the Lesser GNU General Public License as published by  *
 * the Free Software Foundation, either version 3.0 of the License, or          *
 * (at your option) any later version.                                          *
 *                                                                              *
 * IfcOpenShell is distributed in the hope that it will be useful,              *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of               *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                 *
 * Lesser GNU General Public License for more details.                          *
 *                                                                              *
 * You should have received a copy of the Lesser GNU General Public License     *
 * along with this program. If not, see <http://www.gnu.org/licenses/>.         *
 *                                                                              *
 ********************************************************************************/

// This file was generated with the assistance of an AI coding tool.

#include "mapping.h"

#include <set>
#include <vector>

#define mapping POSTFIX_SCHEMA(mapping)
using namespace ifcopenshell::geom;

// Which other walls a wall runs into, and at which of its own ends.
//
// This is what tells the kernel a wall's material layers have somewhere to
// fold to. Only a connection made at one of *our* ends can fold: a wall
// meeting us ATPATH runs into our side, which divides nothing of ours, and a
// connection whose other end is ATPATH is that same case seen from the other
// wall -- it butts into us and its own layers stop dead, so there is no pair
// of layer stacks to reconcile.
//
// The relationship is IfcRelConnectsPathElements, and it reaches a wall
// through two inverse attributes: ConnectedTo holds the ones where we are the
// relating element, ConnectedFrom the ones where we are the related element.
// Which of the two we came through is what says which connection type is ours,
// so the elements need not be compared.

namespace {
    ifcopenshell::geom::connection_type to_connection_type(typename IfcSchema::IfcConnectionTypeEnum::Value v) {
        switch (v) {
        case IfcSchema::IfcConnectionTypeEnum::IfcConnectionType_ATPATH: return ifcopenshell::geom::ATPATH;
        case IfcSchema::IfcConnectionTypeEnum::IfcConnectionType_ATSTART: return ifcopenshell::geom::ATSTART;
        case IfcSchema::IfcConnectionTypeEnum::IfcConnectionType_ATEND: return ifcopenshell::geom::ATEND;
        default: return ifcopenshell::geom::NOTDEFINED;
        }
    }
}

bool mapping::get_wall_neighbours(const express::base& p, std::vector<endpoint_connection>& neighbours) {
    auto element = p.as<IfcSchema::IfcElement>();
    if (!element) {
        return false;
    }

    // A file may state the same junction twice, and the fold must be made once
    // -- built twice it is cut by its own result the second time.
    std::set<express::base> seen;

    auto gather = [&](const std::vector<IfcSchema::IfcRelConnectsElements>& rels, bool we_are_relating) {
        for (auto& rel : rels) {
            auto path = rel.as<IfcSchema::IfcRelConnectsPathElements>();
            if (!path) {
                continue;
            }

            const auto own = to_connection_type(we_are_relating
                ? path.RelatingConnectionType()
                : path.RelatedConnectionType());
            const auto other = to_connection_type(we_are_relating
                ? path.RelatedConnectionType()
                : path.RelatingConnectionType());

            if (own != ifcopenshell::geom::ATSTART && own != ifcopenshell::geom::ATEND) {
                continue;
            }
            if (other == ifcopenshell::geom::ATPATH) {
                continue;
            }

            IfcSchema::IfcElement neighbour = we_are_relating
                ? path.RelatedElement()
                : path.RelatingElement();

            if (!neighbour || !neighbour.declaration().is(IfcSchema::IfcWall::Class())) {
                continue;
            }

            express::base base = neighbour;
            if (base == p || !seen.insert(base).second) {
                continue;
            }

            neighbours.push_back(std::make_tuple(own, other, base));
        }
    };

    gather(element.ConnectedTo(), true);
    gather(element.ConnectedFrom(), false);

    return !neighbours.empty();
}
