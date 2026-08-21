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

#include "mapping.h"

#include "../../ifcparse/si_prefix.h"

#include <optional>
#include <vector>

#define mapping POSTFIX_SCHEMA(mapping)
using namespace ifcopenshell::geom;

// Units belong to a context, not to a file, and a file is allowed to hold more
// than one context: from IFC4 on an IfcProjectLibrary is a context of its own,
// and merged or federated exports put several IfcProject instances in a file
// even though IfcSingleProjectInstance forbids it. Rather than pick one set of
// units for everything, each context's geometry is mapped in that context's
// own units, and anything whose context declares none keeps the units last
// selected.

namespace {
    // RepresentationContexts is mandatory in IFC2x3 and optional from IFC4 on,
    // so the attribute comes back as a list in one and an optional list in the
    // other. Only the overload that matches is instantiated.
    template <typename T>
    const std::vector<T>* as_list(const std::vector<T>& list) {
        return &list;
    }

    template <typename T>
    const std::vector<T>* as_list(const std::optional<std::vector<T>>& list) {
        return list ? &*list : nullptr;
    }
}

std::pair<bool, bool> mapping::read_units_(const IfcSchema::IfcUnitAssignment& assignment, double& length, double& angle, std::string* length_name) {
    bool length_seen = false, angle_seen = false;

    if (!assignment) {
        return { length_seen, angle_seen };
    }

    try {
        auto units = assignment.Units();
        if (units.empty()) {
            logger_.warning("GEO", 310, "No unit information found");
            return { length_seen, angle_seen };
        }
        for (auto& base : units) {
            auto named_unit = base.as<IfcSchema::IfcNamedUnit>();
            if (!named_unit) {
                continue;
            }
            const bool is_length = named_unit.UnitType() == IfcSchema::IfcUnitEnum::IfcUnit_LENGTHUNIT;
            const bool is_angle = named_unit.UnitType() == IfcSchema::IfcUnitEnum::IfcUnit_PLANEANGLEUNIT;
            if (!is_length && !is_angle) {
                continue;
            }
            const double magnitude = ifcopenshell::get_SI_equivalent<IfcSchema>(named_unit);
            if (magnitude == 0.) {
                continue;
            }
            if (is_length) {
                length = magnitude;
                length_seen = true;
                if (length_name) {
                    std::string name;
                    if (auto conversion_based = named_unit.as<IfcSchema::IfcConversionBasedUnit>()) {
                        name = conversion_based.Name();
                    } else if (auto si_unit = named_unit.as<IfcSchema::IfcSIUnit>()) {
                        if (si_unit.Prefix()) {
                            name = IfcSchema::IfcSIPrefix::ToString(*si_unit.Prefix());
                        }
                        name += IfcSchema::IfcSIUnitName::ToString(si_unit.Name());
                    }
                    *length_name = name;
                }
            } else {
                angle = magnitude;
                angle_seen = true;
            }
        }
    } catch (const ifcopenshell::exception& e) {
        logger_.message(ifcopenshell::logger::LOG_ERROR, "GEO", 311, std::string("Failed to determine unit information '") + e.what() + "'");
    }

    return { length_seen, angle_seen };
}

void mapping::build_context_units_() {
#ifdef SCHEMA_HAS_IfcContext
    auto contexts = file_->instances_by_type<IfcSchema::IfcContext>();
#else
    auto contexts = file_->instances_by_type<IfcSchema::IfcProject>();
#endif

    if (contexts.size() < 2) {
        return;
    }

    for (auto& context : contexts) {
        auto assignment = context.UnitsInContext();
        if (!assignment) {
            // No units of its own, so it keeps whatever was last selected.
            continue;
        }

        double length = length_unit_, angle = angle_unit_;
        read_units_(assignment, length, angle, nullptr);

        if (length != length_unit_ || angle != angle_unit_) {
            units_vary_ = true;
        }

        auto representation_contexts = context.RepresentationContexts();
        if (auto* list = as_list(representation_contexts)) {
            for (auto& representation_context : *list) {
                units_by_context_[representation_context.identity()] = { length, angle };
            }
        }
    }

    if (!units_vary_) {
        // Every context says the same thing, so there is nothing to switch
        // between and nothing below this point needs to do any work.
        units_by_context_.clear();
        return;
    }

    // A mapped item carries the units it was built in, so it cannot be handed
    // to a context that measures differently. Geometry is no longer reused
    // between representations, which costs time on exactly the files that
    // would otherwise come out wrong.
    use_caching() = false;

    logger_.warning("GEO", 339, "Contexts of this file do not declare the same units; each is mapped in its own and geometry is not reused between them");
}

void mapping::select_units_(const express::base& context) {
    if (units_by_context_.empty() || !context) {
        return;
    }

    auto it = units_by_context_.find(context.identity());
    if (it != units_by_context_.end()) {
        length_unit_ = it->second.first;
        angle_unit_ = it->second.second;
        return;
    }

    // A sub-context is measured in the units of the context it hangs off.
    if (auto sub_context = context.as<IfcSchema::IfcGeometricRepresentationSubContext>()) {
        select_units_(sub_context.ParentContext());
    }

    // Anything else keeps the units last selected, which is what a context
    // declaring none is entitled to.
}
