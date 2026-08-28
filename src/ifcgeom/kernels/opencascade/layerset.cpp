#include "layerset.h"
#include "opencascade_conversion_result.h"

#include "base_utils.h"
#include "boolean_utils.h"

#include "../../../ifcparse/logger.h"

#include <BRep_Tool.hxx>

#include <TopoDS.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopExp_Explorer.hxx>

#include <Standard_Macro.hxx>
#include <TopoDS_Shape.hxx>
#include <NCollection_List.hxx>

#include <string>

#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <gp_Pln.hxx>
#include <Geom_Plane.hxx>
#include <cstring>

#include <BRepBuilderAPI_MakeShell.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepAlgoAPI_Splitter.hxx>
#include <BRepPrimAPI_MakeHalfSpace.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepOffsetAPI_Sewing.hxx>
#include <BOPAlgo_PaveFiller.hxx>
#include <Standard_Version.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <ShapeFix_Shape.hxx>
#include <NCollection_IncAllocator.hxx>

using namespace ifcopenshell::geom;

namespace {

	double face_area(const TopoDS_Shape& f) {
		GProp_GProps props;
		BRepGProp::SurfaceProperties(f, props);
		return props.Mass();
	}

	void subshapes(const TopoDS_Shape& in, std::list<TopoDS_Shape>& out) {
		TopoDS_Iterator sit(in);
		for (; sit.More(); sit.Next()) {
			out.push_back(sit.Value());
		}
	}

	// Which of the two pieces a single boundary divides a body into carries the
	// first layer: the one nearer the first boundary. The many-boundary branch
	// places pieces by index and gives style 0 to the piece that touches only
	// boundary 1 -- the piece on the side boundary 0 is on -- and this is the
	// same rule measured rather than assumed, because the dividing surface's
	// normal does not point the same way for a plane lifted from a slab's
	// extrusion as for a surface offset from a wall's axis. Assuming either
	// sense swaps the two materials of every two-layer product of the other
	// kind, which is invisible wherever both layers are the same thickness.
	bool nearer_first_boundary(const TopoDS_Shape& a, const TopoDS_Shape& b, const opencascade::handle<Geom_Surface>& first) {
		if (first.IsNull()) {
			return true;
		}
		GProp_GProps ga, gb;
		BRepGProp::VolumeProperties(a, ga);
		BRepGProp::VolumeProperties(b, gb);
		GeomAPI_ProjectPointOnSurf pa(ga.CentreOfMass(), first);
		GeomAPI_ProjectPointOnSurf pb(gb.CentreOfMass(), first);
		if (!pa.IsDone() || !pb.IsDone() || pa.NbPoints() == 0 || pb.NbPoints() == 0) {
			return true;
		}
		return pa.LowerDistance() <= pb.LowerDistance();
	}

	// A shell cannot be divided into layers. Splitting a shell by a face only
	// cuts the shell's own faces along the intersection curve -- it never adds
	// a face lying on the cutting surface -- so every piece comes back with no
	// face on any layer boundary and none of them can be placed. Close the body
	// into a solid first, which is also what the two-layer path does by way of
	// its half-space boolean.
	TopoDS_Shape close_for_splitting(const TopoDS_Shape& shape, double tol) {
		TopoDS_Shape sld = ifcopenshell::geom::util::ensure_fit_for_subtraction(shape, tol);
		if (TopExp_Explorer(sld, TopAbs_SOLID).More() || !TopExp_Explorer(sld, TopAbs_SHELL).More()) {
			return sld;
		}
		TopoDS_Shape solid;
		if (ifcopenshell::geom::util::create_solid_from_compound(sld, solid, tol) && TopExp_Explorer(solid, TopAbs_SOLID).More()) {
			return solid;
		}
		return sld;
	}

#if OCC_VERSION_HEX >= 0x70200
    bool split(const TopoDS_Shape& input, const NCollection_List<TopoDS_Shape>& operands, double eps, std::vector<std::vector<TopoDS_Shape>>& slices) {
		if (operands.Extent() < 2) {
			// Needs to have at least two cutting surfaces for the ordering based on surface containment to work.
			return false;
		}

		BRepAlgoAPI_Splitter split;
        NCollection_List<TopoDS_Shape> input_list;
		input_list.Append(input);
		split.SetArguments(input_list);
		split.SetTools(operands);
		split.SetNonDestructive(true);
		split.SetFuzzyValue(eps);
		split.Build();

		if (!split.IsDone()) {
			return false;
		} else {

			std::map<Geom_Surface*, int> surfaces;

			// NB 1, since first surface has been excluded
			int i = 1;
			for (NCollection_List<TopoDS_Shape>::Iterator it(operands); it.More(); it.Next(), ++i) {
				TopExp_Explorer exp(it.Value(), TopAbs_FACE);
				for (; exp.More(); exp.Next()) {
					surfaces.insert(std::make_pair(BRep_Tool::Surface(TopoDS::Face(exp.Current())).get(), i));
				}
			}

			auto result_shape = split.Shape();
			std::list<TopoDS_Shape> subs;
			subshapes(result_shape, subs);

			// Sometimes there is more nesting of compounds, so when we find a single compound we again try to explode it into a list.
			if (subs.size() == 1 && (subs.front().ShapeType() == TopAbs_COMPSOLID || subs.front().ShapeType() == TopAbs_COMPOUND)) {
				auto s = subs.front();
				subs.clear();
				subshapes(s, subs);
			}

			// One bucket per layer, not per piece. A layer is not always a
			// single piece: an opening, or anything else that divides the body
			// across its whole thickness, leaves each layer in several pieces,
			// and they all belong to the same material. Sizing this by the
			// piece count instead made every such body fail -- the second piece
			// of a layer found the bucket taken and the whole product was
			// handed back unsliced.
			slices.resize(operands.Extent() + 1);

			// Count the ways a piece failed to place, so that a refusal says
			// which of them happened rather than only that one did.
			int unmatched = 0, spanning = 0, empty_layers = 0;

			for (auto& s : subs) {

				// Iterate over the faces of solid to find correspondence to original
				// splitting surfaces. For the outmost slices, there will be a single
				// corresponding surface, because the outmost surfaces that align with
				// the body geometry have not been added as operands. For intermediate
				// slices, two surface indices should be find that should be next to
				// each other in the array of input surfaces.

				TopExp_Explorer exp(s, TopAbs_FACE);
				int min = std::numeric_limits<int>::max();
				int max = std::numeric_limits<int>::min();
				for (; exp.More(); exp.Next()) {
					auto ssrf = BRep_Tool::Surface(TopoDS::Face(exp.Current()));
					auto it = surfaces.find(ssrf.get());
					if (it != surfaces.end()) {
						if (it->second < min) {
							min = it->second;

						}
						if (it->second > max) {
							max = it->second;
						}
					}
				}

				int idx = std::numeric_limits<int>::max();
				if (min != std::numeric_limits<int>::max()) {
					if (min == 1 && max == 1) {
						idx = 0;
					} else if (min + 1 == max || min == max) {
						idx = min;
					}
				}

				if (idx >= 0 && idx < (int)slices.size()) {
					slices[idx].push_back(s);
					continue;
				}

				// A piece could not be placed. Pieces are placed by which
				// cutting surface their faces sit on, compared by surface
				// object, which holds as long as the splitter hands back the
				// object it was given -- it does not when it has had to rebuild
				// the face. Ordering the pieces by where they sit instead was
				// tried and taken back out: it admits pieces that are not
				// layers, on bodies whose own thickness does not match the
				// layer set they carry. What guards that case now is the
				// thickness check the kernel runs over the finished slices.
				if (min == std::numeric_limits<int>::max()) {
					// No face of this piece sits on any cutting surface.
					++unmatched;
				} else {
					// Its faces sit on boundaries that are not neighbours, so
					// it is not one layer thick.
					++spanning;
				}
			}

			for (auto& bucket : slices) {
				if (bucket.empty()) {
					++empty_layers;
				}
			}

			if (unmatched || spanning || empty_layers) {
				// Where each cutting surface sat relative to the body: the
				// signed distance from its plane to the body's furthest points
				// either side. A boundary that does not straddle zero never
				// had anything to cut.
				std::string reach;
				for (NCollection_List<TopoDS_Shape>::Iterator oit(operands); oit.More(); oit.Next()) {
					TopExp_Explorer fe(oit.Value(), TopAbs_FACE);
					if (!fe.More()) {
						continue;
					}
					auto srf = BRep_Tool::Surface(TopoDS::Face(fe.Current()));
					GeomAdaptor_Surface ga(srf);
					if (ga.GetType() != GeomAbs_Plane) {
						reach += " (not a plane)";
						continue;
					}
					gp_Pln pln = ga.Plane();
					double dmin = std::numeric_limits<double>::max();
					double dmax = -std::numeric_limits<double>::max();
					for (TopExp_Explorer ve(input, TopAbs_VERTEX); ve.More(); ve.Next()) {
						double d = pln.SignedDistance(BRep_Tool::Pnt(TopoDS::Vertex(ve.Current())));
						dmin = std::min(dmin, d);
						dmax = std::max(dmax, d);
					}
					reach += " [" + std::to_string(dmin) + "," + std::to_string(dmax) + "]";
				}

				std::string kinds;
				for (auto& s : subs) {
					int nf = 0;
					for (TopExp_Explorer e(s, TopAbs_FACE); e.More(); e.Next(), ++nf) {}
					kinds += " " + std::to_string((int)s.ShapeType()) + "/" + std::to_string(nf) + "f";
				}
				ifcopenshell::logger::root().error("GEO", 171,
					"Unable to map layer geometry to material index: " +
					std::to_string(subs.size()) + " pieces for " +
					std::to_string(operands.Extent() + 1) + " layers, of which " +
					std::to_string(unmatched) + " on no boundary and " +
					std::to_string(spanning) + " not between neighbouring boundaries, leaving " +
					std::to_string(empty_layers) + " layers with no piece" +
					" [input type " + std::to_string((int)input.ShapeType()) +
					", " + std::to_string(surfaces.size()) + " boundary faces; pieces type/faces:" + kinds +
					"; body either side of each boundary:" + reach + "]");
				return false;
			}
		}

		return true;
	}
#else
	bool split(const TopoDS_Shape& input, const TopTools_ListOfShape& operands, double, std::vector<std::vector<TopoDS_Shape>>& slices) {
		TopTools_ListIteratorOfListOfShape it(operands);
		TopoDS_Shape i = input;
		for (; it.More(); it.Next()) {
			const TopoDS_Shape& s = it.Value();
			TopoDS_Shape a, b;

			Handle(Geom_Surface) surf;
			if (s.ShapeType() == TopAbs_FACE) {
				surf = BRep_Tool::Surface(TopoDS::Face(s));
			}

			if ((s.ShapeType() == TopAbs_FACE && ifcopenshell::geom::util::split_solid_by_surface(i, surf, a, b)) ||
				(s.ShapeType() == TopAbs_SHELL && ifcopenshell::geom::util::split_solid_by_shell(i, s, a, b))) {
				slices.push_back({ b });
				i = a;
			} else {
				return false;
			}
		}
		slices.push_back({ i });
		return true;
	}
#endif
}


bool ifcopenshell::geom::util::apply_folded_layerset(const std::vector<conversion_result>& items, const std::vector< std::vector<opencascade::handle<Geom_Surface>>>& surfaces, const std::vector<ifcopenshell::geom::taxonomy::style::ptr>& styles, std::vector<conversion_result>& result, double tol, double widen) {
	Bnd_Box bb;
	TopoDS_Shape input;
	flatten_shape_list(items, input, false, false, tol);

	typedef std::vector< std::vector<opencascade::handle<Geom_Surface>> > folded_surfaces_t;
	typedef std::vector< std::pair< TopoDS_Face, std::pair<gp_Pnt, gp_Pnt> > > faces_with_mass_t;

	NCollection_List<TopoDS_Shape> shells;

	// The first and last boundaries are the wall's own faces: they order the
	// slices but they do not divide anything, and they are never folded. This
	// is the same contract apply_layerset works to, so that the reference
	// boundary below is there to be read -- without it a two layer wall got an
	// arbitrary order and came out with its materials swapped whenever it was
	// folded.
	for (size_t idx = 1; idx + 1 < surfaces.size(); ++idx) {
		folded_surfaces_t::const_iterator it = surfaces.begin() + idx;
		if (it->empty()) {
			continue;
		} else if (it->size() == 1) {
			const opencascade::handle<Geom_Surface>& surface = (*it)[0];
			double u1, v1, u2, v2;
			if (!project(surface, input, u1, v1, u2, v2, widen)) {
				continue;
			}
			shells.Append(BRepBuilderAPI_MakeShell(surface, u1, v1, u2, v2).Shell());
		} else {
			faces_with_mass_t solids;
			// What each surface is and how each cut below changed it, kept so
			// that a failed sew can say which faces never met. A cut that
			// leaves a face whole still counts one face, so the na/nb checks
			// alone cannot tell a real trim from a miss.
			std::vector<std::string> kinds;
			for (folded_surfaces_t::value_type::const_iterator jt = it->begin(); jt != it->end(); ++jt) {
				// The fold surfaces continue this boundary inside a
				// neighbour, so they lie beyond this wall's own body -- by up
				// to the axis gap the junction guard allows plus the
				// neighbour's stack. The default reach of a tenth of a unit
				// left every deeper fold plane just out of range: the face
				// stopped short of it, both half space cuts were no-ops, and
				// the sew below reported honestly that the pieces were apart.
				const opencascade::handle<Geom_Surface>& surface = *jt;
				double u1, v1, u2, v2;
				if (!project(surface, input, u1, v1, u2, v2, widen)) {
					continue;
				}
				TopoDS_Face face = BRepBuilderAPI_MakeFace(surface, u1, u2, v1, v2, 1.e-7).Face();
				gp_Pnt p, p1, p2; gp_Vec vu, vv, n;
				surface->D1((u1 + u2) / 2., (v1 + v2) / 2., p, vu, vv);
				n = vu ^ vv;
				p1 = p.Translated(n);
				p2 = p.Translated(-n);
				solids.push_back(std::make_pair(face, std::make_pair(p1, p2)));

				std::string k(surface->DynamicType()->Name());
				const size_t us = k.rfind('_');
				if (us != std::string::npos) {
					k = k.substr(us + 1);
				}
				char buf[640];
				snprintf(buf, sizeof(buf), "%s u=%.3f..%.3f v=%.3f..%.3f a=%.4f", k.c_str(), u1, u2, v1, v2, face_area(face));
				if (auto gpln = opencascade::handle<Geom_Plane>::DownCast(surface)) {
					const gp_Pnt& o = gpln->Position().Location();
					const gp_Dir& d = gpln->Position().Direction();
					const gp_Dir& xd = gpln->Position().XDirection();
					char buf2[192];
					snprintf(buf2, sizeof(buf2), " o=(%.3f,%.3f,%.3f) n=(%.3f,%.3f,%.3f) x=(%.3f,%.3f,%.3f)",
						o.X(), o.Y(), o.Z(), d.X(), d.Y(), d.Z(), xd.X(), xd.Y(), xd.Z());
					strncat(buf, buf2, sizeof(buf) - strlen(buf) - 1);
				}
				{
					Bnd_Box fb;
					BRepBndLib::Add(face, fb);
					double x1, y1, z1, x2, y2, z2;
					fb.Get(x1, y1, z1, x2, y2, z2);
					char buf3[160];
					snprintf(buf3, sizeof(buf3), " box=(%.3f,%.3f,%.3f)..(%.3f,%.3f,%.3f)", x1, y1, z1, x2, y2, z2);
					strncat(buf, buf3, sizeof(buf) - strlen(buf) - 1);
				}
				kinds.push_back(buf);
			}


			if (solids.empty()) {
				continue;
			}

			int uncut = 0;
			std::string detail;

			faces_with_mass_t::iterator jt = solids.begin();
			TopoDS_Face& A = jt->first;
			TopoDS_Shape An = BRepPrimAPI_MakeHalfSpace(A, jt->second.second).Solid();
			for (++jt; jt != solids.end(); ++jt) {
				TopoDS_Face& B = jt->first;
				TopoDS_Shape Bn = BRepPrimAPI_MakeHalfSpace(B, jt->second.second).Solid();

				const double a_before = face_area(A), b_before = face_area(B);

				TopoDS_Shape a = BRepAlgoAPI_Cut(A, Bn);
				const int na = util::count(a, TopAbs_FACE);
				if (na == 1) {
					A = TopoDS::Face(TopExp_Explorer(a, TopAbs_FACE).Current());
				}

				TopoDS_Shape b = BRepAlgoAPI_Cut(B, An);
				const int nb = util::count(b, TopAbs_FACE);
				if (nb == 1) {
					B = TopoDS::Face(TopExp_Explorer(b, TopAbs_FACE).Current());
				}

				{
					const size_t bi = jt - solids.begin();
					char buf[96];
					snprintf(buf, sizeof(buf), "; cut 0-%d: %.4f->%.4f / %.4f->%.4f",
						(int)bi, a_before, face_area(A), b_before, face_area(B));
					kinds[0] += buf;
				}

				// Either face left whole is a fold that did not turn: the two
				// then cross instead of meeting, and the sew below cannot make
				// a shell of them.
				if (na != 1 || nb != 1) {
					uncut += 1;
					detail = "boundary cut into " + std::to_string(na) + " and fold into " + std::to_string(nb);
				}
			}

			// The kernel's own precision rather than the sewer's 1e-6 default,
			// since the faces arrive from boolean cuts and agree only to the
			// tolerance the rest of the conversion works to. Measured on King
			// this changes nothing either way -- where the sew fails there the
			// two faces are genuinely apart, not just out of tolerance -- so
			// it is the right default, not a fix for that.
			BRepOffsetAPI_Sewing builder(tol);
			for (faces_with_mass_t::const_iterator kt = solids.begin(); kt != solids.end(); ++kt) {
				builder.Add(kt->first);
			}

			builder.Perform();
			TopoDS_Shape s = builder.SewedShape();
			if (s.ShapeType() == TopAbs_SHELL) {
				shells.Append(TopoDS::Shell(s));
			} else {
				// Which faces never meet: the smallest distance between each
				// pair of the faces the sewer was given. A pair a long way
				// apart was left untouched by a cut that was meant to trim it.
				std::string gaps;
				for (size_t gi = 0; gi < solids.size(); ++gi) {
					for (size_t gj = gi + 1; gj < solids.size(); ++gj) {
						BRepExtrema_DistShapeShape dss(solids[gi].first, solids[gj].first);
						char buf[64];
						if (dss.IsDone()) {
							snprintf(buf, sizeof(buf), " %d-%d=%.4g", (int)gi, (int)gj, dss.Value());
						} else {
							snprintf(buf, sizeof(buf), " %d-%d=?", (int)gi, (int)gj);
						}
						gaps += buf;
					}
				}
				std::string faces_detail;
				for (size_t ki = 0; ki < kinds.size(); ++ki) {
					faces_detail += (ki ? " | " : "") + std::to_string(ki) + ": " + kinds[ki];
				}
				ifcopenshell::logger::root().error("GEO", 172,
					"Sewing " + std::to_string(solids.size()) + " surfaces of a folded layer boundary gave "
					+ std::to_string(util::count(s, TopAbs_FACE)) + " face(s) that are not a shell"
					+ (uncut ? ", " + std::to_string(uncut) + " of the pairs did not trim each other (" + detail + ")"
					         : ", though every pair trimmed each other")
					+ " [" + faces_detail + "; gaps:" + gaps + "]");
				return false;
			}
		}
	}

	if (shells.Extent() == 0) {

		return false;

	} else if (shells.Extent() == 1) {

		for (std::vector<conversion_result>::const_iterator it = items.begin(); it != items.end(); ++it) {
			TopoDS_Shape a, b;
			if (split_solid_by_shell(std::static_pointer_cast<open_cascade_shape>(it->shape())->shape(), shells.First(), a, b, tol)) {
				opencascade::handle<Geom_Surface> first;
				if (!surfaces.empty() && !surfaces.front().empty()) {
					first = surfaces.front().front();
				}
				const bool a_is_first = nearer_first_boundary(a, b, first);
				const TopoDS_Shape& s0 = a_is_first ? a : b;
				const TopoDS_Shape& s1 = a_is_first ? b : a;
				result.push_back(conversion_result(it->ItemId(), it->placement(), new open_cascade_shape(s0), (!!styles[0] ? styles[0] : it->style_ptr())));
				result.push_back(conversion_result(it->ItemId(), it->placement(), new open_cascade_shape(s1), (!!styles[1] ? styles[1] : it->style_ptr())));
			} else {
				continue;
			}
		}

		return true;

	} else {

		for (std::vector<conversion_result>::const_iterator it = items.begin(); it != items.end(); ++it) {

			const TopoDS_Shape& s = std::static_pointer_cast<open_cascade_shape>(it->shape())->shape();
			TopoDS_Shape sld = close_for_splitting(s, tol);

			std::vector<std::vector<TopoDS_Shape>> slices;
			if (split(sld, shells, tol, slices) && slices.size() == styles.size()) {
				for (size_t i = 0; i < slices.size(); ++i) {
					for (auto& piece : slices[i]) {
						result.push_back(conversion_result(it->ItemId(), it->placement(), new open_cascade_shape(piece), (!!styles[i] ? styles[i] : it->style_ptr())));
					}
				}
			} else {
				return false;
			}
		}

		return true;

	}

}

bool ifcopenshell::geom::util::apply_layerset(const std::vector<conversion_result>& items, const std::vector<opencascade::handle<Geom_Surface>>& surfaces, const std::vector<ifcopenshell::geom::taxonomy::style::ptr>& styles, std::vector<conversion_result>& result, double tol) {
	if (surfaces.size() < 3) {

		return false;

	} else if (surfaces.size() == 3) {

		for (std::vector<conversion_result>::const_iterator it = items.begin(); it != items.end(); ++it) {
			TopoDS_Shape a, b;
			if (split_solid_by_surface(std::static_pointer_cast<open_cascade_shape>(it->shape())->shape(), surfaces[1], a, b, tol)) {
				const bool a_is_first = nearer_first_boundary(a, b, surfaces[0]);
				const TopoDS_Shape& s0 = a_is_first ? a : b;
				const TopoDS_Shape& s1 = a_is_first ? b : a;
				result.push_back(conversion_result(it->ItemId(), it->placement(),new open_cascade_shape(s0), (!!styles[0] ? styles[0] : it->style_ptr())));
				result.push_back(conversion_result(it->ItemId(), it->placement(),new open_cascade_shape(s1), (!!styles[1] ? styles[1] : it->style_ptr())));
			} else {
				// Skipping the item would report success while handing back a
				// result this body is missing from, and the caller has no way
				// to tell that from a body that was genuinely divided.
				ifcopenshell::logger::root().warning("GEO", 334, "Unable to divide a body by its single layer boundary");
				return false;
			}
		}

		return true;

	} else {

		/*
		// Determine whether sequence of surfaces is consistent with surface normal, so that
		// layer operations are applied in the correct order. This seems to be always the case.
		Bnd_Box bb;
		for (std::vector<conversion_result>::const_iterator it = items.begin(); it != items.end(); ++it) {
			BRepBndLib::Add(it->shape(), bb);
		}

		double x1, y1, z1, x2, y2, z2;
		bb.Get(x1, y1, z1, x2, y2, z2);
		gp_Pnt p1(x1, y1, z1);
		gp_Pnt p2(x2, y2, z2);
		gp_Pnt avg = (p1.XYZ() + p2.XYZ()) / 2.;

		ShapeAnalysis_Surface sas1(surfaces[0]);
		ShapeAnalysis_Surface sas2(surfaces[1]);
		const gp_Pnt2d uv = sas1.ValueOfUV(avg, 1e-3);

		gp_Pnt ps1, ps2, mass;
		gp_Vec du1, dv1, du2, dv2;
		surfaces[0]->D1(uv.X(), uv.Y(), ps1, du1, dv1);
		const gp_Vec n1 = dv1.XYZ() ^ du1.XYZ();

		const bool reversed = gp_Dir(ps2.XYZ() - ps1.XYZ()).Dot(n1) < 0.;

		surfaces[surfaces.size() - 1]->D0(uv.X(), uv.Y(), mass);
		mass.ChangeCoord() += n1.XYZ();
		*/

		for (std::vector<conversion_result>::const_iterator it = items.begin(); it != items.end(); ++it) {

			const TopoDS_Shape& s = std::static_pointer_cast<open_cascade_shape>(it->shape())->shape();
			TopoDS_Shape sld = close_for_splitting(s, tol);
			NCollection_List<TopoDS_Shape> operands;
			for (unsigned i = 1; i < surfaces.size() - 1; ++i) {
				double u1, v1, u2, v2;
				if (!project(surfaces[i], sld, u1, v1, u2, v2)) {
					ifcopenshell::logger::root().warning("GEO", 332, "Unable to fit layer boundary " + std::to_string(i) + " of " + std::to_string(surfaces.size()) + " to the body it divides");
					return false;
				}

				TopoDS_Face face = BRepBuilderAPI_MakeFace(surfaces[i], u1, u2, v1, v2, 1.e-7).Face();
				operands.Append(face);
			}

			/*
			// enable this is you want to see how IfcOpenShell has placed the layer surfaces
			for (auto& x : operands) {
				result.push_back(conversion_result(it->ItemId(), it->placement(), x, nullptr));
			}
			*/

			std::vector<std::vector<TopoDS_Shape>> slices;
			if (split(sld, operands, tol, slices) && slices.size() == styles.size()) {
				for (size_t i = 0; i < slices.size(); ++i) {
					for (auto& piece : slices[i]) {
						result.push_back(conversion_result(it->ItemId(), it->placement(), new open_cascade_shape(piece), (!!styles[i] ? styles[i] : it->style_ptr())));
					}
				}
			} else {
				ifcopenshell::logger::root().warning("GEO", 333, "Splitting a body by " + std::to_string(operands.Extent()) + " layer boundaries did not give one piece for each of its " + std::to_string(styles.size()) + " layers");
				return false;
			}
		}

		return true;
	}
}


bool ifcopenshell::geom::util::split_solid_by_surface(const TopoDS_Shape& input, const opencascade::handle<Geom_Surface>& surface, TopoDS_Shape& front, TopoDS_Shape& back, double tol) {
	// Use an unbounded surface, that isolate part of the input shape,
	// to split this shape into two parts. Make sure that the addition
	// of the two result volumes matches that of the input.

	double u1, v1, u2, v2;
	if (!project(surface, input, u1, v1, u2, v2)) {
		return false;
	}

	TopoDS_Face face = BRepBuilderAPI_MakeFace(surface, u1, u2, v1, v2, 1.e-7).Face();
	gp_Pnt p, p1, p2; gp_Vec vu, vv, n;
	surface->D1((u1 + u2) / 2., (v1 + v2) / 2., p, vu, vv);
	n = vu ^ vv;
	p1 = p.Translated(-n);
	TopoDS_Solid solid = BRepPrimAPI_MakeHalfSpace(face, p1).Solid();

	const bool b = split_solid_by_shell(input, solid, front, back, tol);
	return b;
}

bool ifcopenshell::geom::util::split_solid_by_shell(const TopoDS_Shape& input, const TopoDS_Shape& shell, TopoDS_Shape& front, TopoDS_Shape& back, double tol) {
	// Use a shell, typically one or more connected faces, that isolate part
	// of the input shape, to split this shape into two parts. Make sure that
	// the addition of the two result volumes matches that of the input.

	TopoDS_Solid solid;
	if (shell.ShapeType() == TopAbs_SHELL) {
		solid = BRepBuilderAPI_MakeSolid(TopoDS::Shell(shell)).Solid();
	} else if (shell.ShapeType() == TopAbs_SOLID) {
		solid = TopoDS::Solid(shell);
	} else {
		return false;
	}

#if OCC_VERSION_HEX >= 0x70300
	NCollection_List<TopoDS_Shape> shapes;
#else
	BOPCol_ListOfShape shapes;
#endif
	shapes.Append(input);
	shapes.Append(solid);
	BOPAlgo_PaveFiller filler(new NCollection_IncAllocator); // TODO: Does this need to be freed?
	filler.SetArguments(shapes);
	filler.Perform();
	front = BRepAlgoAPI_Cut(input, solid, filler);
	back = BRepAlgoAPI_Common(input, solid, filler);

	bool is_null[2];

	for (int i = 0; i < 2; ++i) {
		TopoDS_Shape& shape = i == 0 ? front : back;
		const bool result_is_null = is_null[i] = shape.IsNull() != 0;
		if (result_is_null) {
			continue;
		}
		try {
			ShapeFix_Shape fix(shape);
			if (fix.Perform()) {
				shape = fix.Shape();
			}
		} catch (const Standard_Failure& e) {
			if (e.GetMessageString() && strlen(e.GetMessageString())) {
				ifcopenshell::logger::root().error("GEO", 173, e.GetMessageString());
			} else {
				ifcopenshell::logger::root().error("GEO", 174, "Unknown error performing fixes");
			}
		} catch (...) {
			ifcopenshell::logger::root().error("GEO", 175, "Unknown error performing fixes");
		}
		BRepCheck_Analyzer analyser(shape);
		bool is_valid = analyser.IsValid() != 0;
		if (!is_valid) {
			return false;
		}
	}

	if (is_null[0] || is_null[1]) {
		ifcopenshell::logger::root().message(ifcopenshell::logger::LOG_ERROR, "GEO", 176, "Null result obtained from layerset slicing");
		if (is_null[0] && is_null[1]) {
			return false;
		}
	}

	const double ab = shape_volume(input);
	const double a = shape_volume(front);
	const double b = shape_volume(back);

	return std::fabs(ab - (a + b)) < tol;
}
