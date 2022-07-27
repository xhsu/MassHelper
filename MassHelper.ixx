module;

// C
#include <cassert>
#include <cfloat>

// C++
#include <algorithm>
#include <array>
#include <concepts>
#include <deque>
#include <functional>
#include <tuple>
#include <unordered_map>
#include <vector>

// C++ Exp
#include <experimental/generator>

// Others
#include <fmt/color.h>
#include <gcem.hpp>

// Platforms
#include <fcntl.h>
#include <io.h>

export module MassHelper;

import UtlConcepts;

import PeriodicTable;

using std::array;
using std::deque;
using std::function;
using std::list;
using std::pair;
using std::string;
using std::tuple;
using std::unordered_map;
using std::vector;

using std::experimental::generator;

export enum AminoAcids_e : char
{
	NOT_AN_AMINO_ACID = '\0',

	Glycine = 'G',
	Alanine = 'A',
	Serine = 'S',
	Proline = 'P',
	Valine = 'V',
	Threonine = 'T',
	Cysteine = 'C',
	Isoleucine = 'I',
	Leucine = 'L',
	Asparagine = 'N',
	Aspartic_Acid = 'D',
	Lysine = 'K',
	Glutamine = 'Q',
	Glutamic_Acid = 'E',
	Methionine = 'M',
	Histidine = 'H',
	Methionine_Sulfoxide = 'm',
	Phenylalanine = 'F',
	Arginine = 'R',
	Carbamidomethyl_Cysteine = 'ç',
	Carboxyethyl_Cysteine = 'c',
	Tyrosine = 'Y',
	Tryptophan = 'W',
};

enum ProteinModifications_e : unsigned char
{
	NO_MODIFICATION = 0,
	GLU_CYCLIZATION,
	METHYLATION,
	ACETYLATION,
	PHOSPHORYLATION,
};

export enum IgnoreLevel_e : unsigned char
{
	IGNORE_ALL_MODIFICATION = 0,
	IGNORE_UNCOMMON,	// don't ignore Methionine_Sulfoxide, Carbamidomethyl_Cysteine and Carboxyethyl_Cysteine
	IGNORE_NOTHING,
};

constexpr auto g_rgflModMassShift = array
{
	0.0,
	-Monoisotopic::MWt<"NH3">,	// Loss during free amino group of glutamic acid or glutamine cyclizes to form a lactam.
	Monoisotopic::MWt<"CH2">,	// Methylation: Plug an -CH3 group onto N, in replace a hydrogen atom.
	Monoisotopic::MWt<"CH2CO">,	// Acetylation: Replace a hydrogen atom with CH3CO-
	Monoisotopic::MWt<"OPOOH">,	// Phosphorylation (O-phos): Replace ??? with -OPOOH
};

constexpr auto g_rgpszModName = array
{
	"",
	"(cy.)",
	"(met.)",
	"(ace.)",
	"(pho.)",
};

struct Molecule_t
{
	double m_Mass = 0;
	const char* m_Name = "none";

	template <StringLiteral FORMULA>
	static consteval Molecule_t Monoisotope(void) noexcept
	{
		return Molecule_t{ Monoisotopic::MWt<FORMULA>, FORMULA };
	}
};

struct AminoAcid_t
{
	double m_ResidueMass = 0;
	int m_ImmoniumIonMass = 0;
	const char* m_3Letters = "err";
	Molecule_t m_CharacteristicNeuLoss{};
};

constexpr auto g_rgCommonNeuLoss = array
{
	Molecule_t::Monoisotope<"H2O">(),
	Molecule_t::Monoisotope<"NH3">(),
	Molecule_t::Monoisotope<"HPO3">(),
	Molecule_t::Monoisotope<"H3PO4">(),
};

constexpr Molecule_t NOTHINGNESS = {};
constexpr Molecule_t WATER = Molecule_t::Monoisotope<"H2O">();
constexpr Molecule_t AMMONIA = Molecule_t::Monoisotope<"NH3">();

const unordered_map<std::underlying_type_t<AminoAcids_e>, AminoAcid_t> g_rgAminoAcidsData =
{
	{AminoAcids_e::NOT_AN_AMINO_ACID,		{}},
	{AminoAcids_e::Glycine,					{57.0215, 30, "Gly", NOTHINGNESS}},
	{AminoAcids_e::Alanine,					{71.0371, 44, "Ala", NOTHINGNESS}},
	{AminoAcids_e::Serine,					{87.0320, 60, "Ser", WATER}},
	{AminoAcids_e::Proline,					{97.0528, 70, "Pro", NOTHINGNESS}},
	{AminoAcids_e::Valine,					{99.0684, 72, "Val", NOTHINGNESS}},
	{AminoAcids_e::Threonine,				{101.0477, 74, "Thr", WATER}},
	{AminoAcids_e::Cysteine,				{103.0092, 76, "Cys", Molecule_t::Monoisotope<"HS2">()}},
	{AminoAcids_e::Isoleucine,				{113.0841, 86, "Ile", NOTHINGNESS}},
	{AminoAcids_e::Leucine,					{113.0841, 86, "Leu", NOTHINGNESS}},
	{AminoAcids_e::Asparagine,				{114.0429, 87, "Asn", AMMONIA}},
	{AminoAcids_e::Aspartic_Acid,			{115.0269, 88, "Asp", WATER}},
	{AminoAcids_e::Lysine,					{128.0950, 101, "Lys", AMMONIA}},
	{AminoAcids_e::Glutamine,				{128.0586, 101, "Gln", AMMONIA}},
	{AminoAcids_e::Glutamic_Acid,			{129.0426, 102, "Glu", WATER}},
	{AminoAcids_e::Methionine,				{131.0405, 104, "Met", Molecule_t::Monoisotope<"CH3SH">()}},
	{AminoAcids_e::Histidine,				{137.0589, 110, "His", NOTHINGNESS}},
	{AminoAcids_e::Methionine_Sulfoxide,	{147.0354, 120, "Mox", Molecule_t::Monoisotope<"CH3SOH">()}},
	{AminoAcids_e::Phenylalanine,			{147.0684, 120, "Phe", NOTHINGNESS}},
	{AminoAcids_e::Arginine,				{156.0684, 129, "Arg", AMMONIA}},
	{AminoAcids_e::Carbamidomethyl_Cysteine,{160.0307, 133, "Cac", Molecule_t::Monoisotope<"HSCH2CONH2">()}},
	{AminoAcids_e::Carboxyethyl_Cysteine,	{161.0147, 134, "Cec", Molecule_t::Monoisotope<"HSCH2COOH">()}},
	{AminoAcids_e::Tyrosine,				{163.0633, 134, "Tyr", NOTHINGNESS}},
	{AminoAcids_e::Tryptophan,				{186.0793, 159, "Trp", NOTHINGNESS}},
};

export struct MassPeak_t
{
	constexpr MassPeak_t(void) noexcept {}
	constexpr MassPeak_t(Arithmetic auto v) noexcept : m_Value(v) {}
	MassPeak_t& operator=(const MassPeak_t& rhs) noexcept { m_Value = rhs.m_Value; m_String = rhs.m_String; return *this; }
	MassPeak_t& operator=(MassPeak_t&& rhs) noexcept { m_Value = rhs.m_Value; std::swap(m_String, rhs.m_String); return *this; }
	MassPeak_t(const MassPeak_t& rhs) noexcept : m_Value(rhs.m_Value), m_String(rhs.m_String) {}
	MassPeak_t(MassPeak_t&& rhs) noexcept : m_Value(rhs.m_Value), m_String(std::move(rhs.m_String)) {}
	~MassPeak_t(void) noexcept {}

	double m_Value = 0;
	string m_String {};

	void Identify(const string& sz) noexcept
	{
		if (!m_String.empty())
			fmt::print(fg(fmt::color::dark_red), "Ion {} was previously identified as {}.\n", sz, m_String);

		m_String = sz;
	}

	constexpr bool Approx(double v) const noexcept { return !(bool)gcem::round(m_Value - v); }

	constexpr auto operator<=> (Arithmetic auto const& rhs) const noexcept { return m_Value <=> rhs; }
	constexpr bool operator== (Arithmetic auto const& rhs) const noexcept { return gcem::abs(m_Value - rhs) <= DBL_EPSILON; }
	constexpr operator double& () noexcept { return m_Value; }
	constexpr operator const double& () const noexcept { return m_Value; }
};

struct Cell_t
{
	constexpr Cell_t(void) noexcept {}
	constexpr Cell_t(AminoAcids_e what, ProteinModifications_e mod, double b = 0, double y = 0) noexcept : m_AminoAcid(what), m_Modification(mod), m_bTypeIon(b), m_yTypeIon(y) {}

	double m_aTypeIon = 0;	// Subtract -CO-
	double m_bNeutralLostIon = 0;
	double m_bAsteriskIon = 0;	// -Ammonia
	double m_bCircleIon = 0;	// -Water
	double m_bTypeIon = 0;	// Including this amino acid Coule be b[n] i.e. [M+1] ion itself.
	AminoAcids_e m_AminoAcid = NOT_AN_AMINO_ACID;
	ProteinModifications_e m_Modification = NO_MODIFICATION;
	double m_yTypeIon = 0;	// Including this amino acid. Could be y[n] i.e. [M+1] ion itself.
	double m_yCircleIon = 0;	// -Water
	double m_yAsteriskIon = 0;	// -Ammonia
	double m_yNeutralLostIon = 0;

	void Deduce(void) noexcept
	{
		auto const fnGetPhosphorylationNeuLoss = [&]
		{
			switch (m_AminoAcid)
			{
			case Serine:
			case Threonine:
				return Monoisotopic::MWt<"H3PO4">;

			case Tyrosine:
			default:
				// "Could also happening to others residues, including Lys and Arg.
				// Typically, the neutral loss of 80Da for these other amino acids is quite labile,
				// so these modifications can be difficult to analyze.
				return Monoisotopic::MWt<"HPO3">;
			}
		};
		auto const iNeutralMass = g_rgAminoAcidsData.at(m_AminoAcid).m_CharacteristicNeuLoss.m_Mass;

		if (!m_aTypeIon)
			m_aTypeIon = m_bTypeIon - Monoisotopic::MWt<"CO">;
		if (!m_bNeutralLostIon && m_AminoAcid != NOT_AN_AMINO_ACID && iNeutralMass && iNeutralMass != 17 && iNeutralMass != 18)
			m_bNeutralLostIon = m_bTypeIon - g_rgAminoAcidsData.at(m_AminoAcid).m_CharacteristicNeuLoss.m_Mass;
		if (!m_bAsteriskIon)
			m_bAsteriskIon = m_bTypeIon - Monoisotopic::MWt<"NH3">;
		if (!m_bCircleIon)
			m_bCircleIon = m_bTypeIon - Monoisotopic::MWt<"H2O">;

		if (!m_yCircleIon)
			m_yCircleIon = m_yTypeIon - Monoisotopic::MWt<"H2O">;
		if (!m_yAsteriskIon)
			m_yAsteriskIon = m_yTypeIon - Monoisotopic::MWt<"NH3">;
		if (!m_yNeutralLostIon && m_AminoAcid != NOT_AN_AMINO_ACID && iNeutralMass && iNeutralMass != 17 && iNeutralMass != 18)
			m_yNeutralLostIon = m_yTypeIon - g_rgAminoAcidsData.at(m_AminoAcid).m_CharacteristicNeuLoss.m_Mass;

		if (m_Modification == PHOSPHORYLATION) [[unlikely]]	// Override
		{
			m_bNeutralLostIon = m_bTypeIon - fnGetPhosphorylationNeuLoss();
			m_yNeutralLostIon = m_yTypeIon - fnGetPhosphorylationNeuLoss();
		}
	}

	void Disambiguate(const auto& rgflMassPeaks) noexcept
	{
		switch (m_AminoAcid)
		{
			// Nothing we can do about Leucine and Isoleucine.
			//case Leucine:
			//case Isoleucine:

			// Nothing we can do about Glutamine and Lysine.
			// their neutral losses are both NH3.
			//case Lysine:
			//case Glutamine:

		case Phenylalanine:
		case Methionine_Sulfoxide:
		{
			auto const itBegin = rgflMassPeaks.cbegin(), itEnd = rgflMassPeaks.cend();
			auto const itYNeuPeak = std::find_if(itBegin, itEnd, std::bind(&MassPeak_t::Approx, std::placeholders::_1, m_yTypeIon - g_rgAminoAcidsData.at(Methionine_Sulfoxide).m_CharacteristicNeuLoss.m_Mass));	// #UPDATE_AT_CPP23 std::bind_front
			auto const itBNeuPeak = std::find_if(itBegin, itEnd, std::bind(&MassPeak_t::Approx, std::placeholders::_1, m_bTypeIon - g_rgAminoAcidsData.at(Methionine_Sulfoxide).m_CharacteristicNeuLoss.m_Mass));
			bool const bIsMetSf = (itYNeuPeak != itEnd) && (itBNeuPeak != itEnd);

			m_AminoAcid = bIsMetSf ? Methionine_Sulfoxide : Phenylalanine;
			m_yNeutralLostIon = bIsMetSf ? itYNeuPeak->m_Value : 0;
			m_bNeutralLostIon = bIsMetSf ? itBNeuPeak->m_Value : 0;
			break;
		}
		default:
			break;
		}
	}

	constexpr bool Filled(void) const noexcept { return m_AminoAcid != NOT_AN_AMINO_ACID && (m_bTypeIon != 0 || m_yTypeIon != 0); }
	constexpr bool Explained(double flPeak) const noexcept { return !(bool)gcem::round(flPeak - m_aTypeIon) || !(bool)gcem::round(flPeak - m_bNeutralLostIon) || !(bool)gcem::round(flPeak - m_bAsteriskIon) || !(bool)gcem::round(flPeak - m_bCircleIon) || !(bool)gcem::round(flPeak - m_bTypeIon) || !(bool)gcem::round(flPeak - m_yTypeIon) || !(bool)gcem::round(flPeak - m_yCircleIon) || !(bool)gcem::round(flPeak - m_yAsteriskIon) || !(bool)gcem::round(flPeak - m_yNeutralLostIon); }

	constexpr bool operator==(const Cell_t& rhs) const noexcept { return m_AminoAcid == rhs.m_AminoAcid; }
};

export struct AlternativeReality_t
{
	AlternativeReality_t(void) noexcept {}
	AlternativeReality_t& operator=(const AlternativeReality_t& rhs) noexcept { m_MPlusOne = rhs.m_MPlusOne; m_PendingPeaks = rhs.m_PendingPeaks; m_Solution = rhs.m_Solution; return *this; }
	AlternativeReality_t& operator=(AlternativeReality_t&& rhs) noexcept { m_MPlusOne = rhs.m_MPlusOne; std::swap(m_PendingPeaks, rhs.m_PendingPeaks); std::swap(m_Solution, rhs.m_Solution); return *this; }
	AlternativeReality_t(const AlternativeReality_t& rhs) noexcept : m_MPlusOne(rhs.m_MPlusOne), m_PendingPeaks(rhs.m_PendingPeaks), m_Solution(rhs.m_Solution) {}
	AlternativeReality_t(AlternativeReality_t&& rhs) noexcept : m_MPlusOne(rhs.m_MPlusOne), m_PendingPeaks(std::move(rhs.m_PendingPeaks)), m_Solution(std::move(rhs.m_Solution)) {}
	~AlternativeReality_t(void) noexcept {}

	double m_MPlusOne = 0;
	vector<MassPeak_t> m_PendingPeaks{};
	deque<Cell_t> m_Solution{};	// Why deque? We need both emplace_front() and index-based access.
};

export template<uint16 iChargeFrom, uint16 iChargeTo>
constexpr decltype(auto) M2ZConversion(std::floating_point auto flFrom) noexcept
{
	using T = decltype(flFrom);

	static constexpr int iChargeDiff = iChargeTo - iChargeFrom;

	flFrom *= T(iChargeFrom);
	flFrom += T(iChargeDiff) * Monoisotopic::Hydrogen;
	flFrom /= T(iChargeTo);

	return flFrom;
}

pair<AminoAcids_e, ProteinModifications_e> TestNumber(double flPeakDiff, IgnoreLevel_e iIgnoreLv = IGNORE_UNCOMMON, double const flTolerance = 0.5) noexcept
{
	using Ty = tuple<decltype(g_rgAminoAcidsData)::key_type, ProteinModifications_e, double>;
	auto const fn = [&]() -> vector<Ty>
	{
		vector<Ty> ret{};

		for (const auto& [AminoAcid, Data] : g_rgAminoAcidsData)
		{
			if (AminoAcid == NOT_AN_AMINO_ACID
				|| ((AminoAcid == Methionine_Sulfoxide || AminoAcid == Carbamidomethyl_Cysteine || AminoAcid == Carboxyethyl_Cysteine) && iIgnoreLv == IGNORE_ALL_MODIFICATION)
				)
				continue;

			// Exact match on residue
			if (double const diff = std::abs(flPeakDiff - Data.m_ResidueMass); diff < flTolerance)
				ret.emplace_back(AminoAcid, NO_MODIFICATION, diff);

			// It doesn't make any sense to have further modifications.
			if (iIgnoreLv <= IGNORE_UNCOMMON || AminoAcid == Methionine_Sulfoxide || AminoAcid == Carboxyethyl_Cysteine || AminoAcid == Carbamidomethyl_Cysteine)	// It doesn't make any sense to have further modifications.
				continue;

			// Common amino acid modifications
			// Skip the first possibility: Pyroglutamic acid formed from N-terminal Gln
			// Parent mass shift (Da): -17
			for (int i = METHYLATION; i < g_rgflModMassShift.size(); ++i)
				if (auto diff = std::abs(flPeakDiff - (Data.m_ResidueMass + g_rgflModMassShift[i])); diff < flTolerance)
					ret.emplace_back(AminoAcid, (ProteinModifications_e)i, diff);
		}

		return ret;
	};

	static constexpr auto fnCompLesser = [](const Ty& lhs, const Ty& rhs) constexpr -> bool
	{
		auto const& [lty, lmod, ldiff] = lhs;
		auto const& [rty, rmod, rdiff] = rhs;

		if (gcem::abs(ldiff - rdiff) < 0.01 && lmod != rmod)
			return lmod == NO_MODIFICATION;

		return ldiff < rdiff;
	};

	vector<Ty> ret = fn();
	if (auto const it = std::min_element(ret.cbegin(), ret.cend(), fnCompLesser); it != ret.cend())
		return std::make_pair((AminoAcids_e)std::get<0>(*it), std::get<1>(*it));

	return std::make_pair(NOT_AN_AMINO_ACID, NO_MODIFICATION);
}

export list<pair<string, double>> SummaryPeakAsList(double flPeakDiff, IgnoreLevel_e iIgnoreLv = IGNORE_UNCOMMON, double const flTolerance = 0.5) noexcept
{
	using ListTy = list<pair<string, double>>;

	if (flPeakDiff <= DBL_EPSILON)
		return { { "SELF", 0.0 } };

	ListTy rgFound{};

	for (const auto& [AminoAcid, Data] : g_rgAminoAcidsData)
	{
		if (AminoAcid == NOT_AN_AMINO_ACID
			|| ((AminoAcid == Methionine_Sulfoxide || AminoAcid == Carbamidomethyl_Cysteine || AminoAcid == Carboxyethyl_Cysteine) && iIgnoreLv == IGNORE_ALL_MODIFICATION)
			)
			continue;

		if (auto diff = flPeakDiff - Data.m_ResidueMass; std::abs(diff) < flTolerance)
			rgFound.emplace_back(fmt::format("{}({})", Data.m_3Letters, AminoAcid), diff);

		if (auto diff = flPeakDiff - Data.m_CharacteristicNeuLoss.m_Mass; Data.m_CharacteristicNeuLoss.m_Mass > 0 && std::abs(diff) < flTolerance)
			rgFound.emplace_back(Data.m_CharacteristicNeuLoss.m_Name, diff);

		// It doesn't make any sense to have further modifications.
		if (iIgnoreLv <= IGNORE_UNCOMMON || AminoAcid == Methionine_Sulfoxide || AminoAcid == Carboxyethyl_Cysteine || AminoAcid == Carbamidomethyl_Cysteine)
			continue;

		// Common amino acid modifications
		// Skip the first possibility: Pyroglutamic acid formed from N-terminal Gln
		// Parent mass shift (Da): -17

		for (int i = METHYLATION; i < g_rgflModMassShift.size(); ++i)
			if (auto diff = flPeakDiff - (Data.m_ResidueMass + g_rgflModMassShift[i]); std::abs(diff) < flTolerance)
				rgFound.emplace_back(fmt::format("{}({}){}", Data.m_3Letters, AminoAcid, g_rgpszModName[i]), diff);
	}

	for (const auto& Molecule : g_rgCommonNeuLoss)
	{
		if (auto diff = flPeakDiff - Molecule.m_Mass; std::abs(diff) < flTolerance)
			rgFound.emplace_back(Molecule.m_Name, diff);
	}

	rgFound.sort([](const ListTy::value_type& lhs, const ListTy::value_type& rhs) -> bool { return std::abs(lhs.second) < std::abs(rhs.second); });
	rgFound.unique([](const ListTy::value_type& lhs, const ListTy::value_type& rhs) -> bool { return lhs.first == rhs.first; });
	return rgFound;
}

export string SummaryPeakAsString(double flPeakDiff, unsigned iDecimal, IgnoreLevel_e iIgnoreLv = IGNORE_UNCOMMON, bool const bOneLine = true) noexcept
{
	if (flPeakDiff <= DBL_EPSILON)
		return "SELF";

	if (auto const rgPossibilities = SummaryPeakAsList(flPeakDiff, iIgnoreLv); !rgPossibilities.empty())
	{
		string sz{};
		for (const auto& [szInfo, flDelta] : rgPossibilities)
			sz += fmt::format("{0}[{2:.{1}f}]{3}", szInfo, iDecimal, flDelta, bOneLine ? ", " : " \n");

		sz.pop_back();
		sz.pop_back();
		return sz;
	}
	else
		return "-";
}

export double CalcMWtBySeq(const string& szSeq) noexcept
{
	double ret = Monoisotopic::MWt<"H2O">;	// Add H to -NH-, Add -OH to -CO-

	for (auto AminoAcid : szSeq)
	{
		switch (AminoAcid)
		{
		case Glycine:
		case Alanine:
		case Serine:
		case Proline:
		case Valine:
		case Threonine:
		case Cysteine:
		case Isoleucine:
		case Leucine:
		case Asparagine:
		case Aspartic_Acid:
		case Lysine:
		case Glutamine:
		case Glutamic_Acid:
		case Methionine:
		case Histidine:
		case Methionine_Sulfoxide:
		case Phenylalanine:
		case Arginine:
		case Carbamidomethyl_Cysteine:
		case Carboxyethyl_Cysteine:
		case Tyrosine:
		case Tryptophan:
			ret += g_rgAminoAcidsData.at(AminoAcid).m_ResidueMass;
			break;

		default:
			break;
		}
	}

	return ret;
}

auto CalcMWtByCells(const auto& Solution) noexcept
{
	auto const szSeq = Conclude(Solution);
	auto flBaseMWt = CalcMWtBySeq(szSeq);

	for (auto const& Cell : Solution)
		flBaseMWt += g_rgflModMassShift[Cell.m_Modification];

	return flBaseMWt;
}

// Get a polypeptide sequence as string.
export template<bool bFirmlySure = false>
string Conclude(const auto& rgCells) noexcept
{
	string ret;
	ret.resize(rgCells.size());
	ret.clear();

	for (const auto& Cell : rgCells)
	{
		if constexpr (!bFirmlySure)
		{
			switch (Cell.m_AminoAcid)
			{
			case NOT_AN_AMINO_ACID:
				continue;
	
			case Isoleucine:
			case Leucine:
				ret.push_back(Isoleucine);
				break;
	
			case Lysine:
			case Glutamine:
				ret.push_back(Glutamine);
				break;
	
			case Methionine_Sulfoxide:
			case Phenylalanine:
				ret.push_back(Phenylalanine);
				break;
	
			default:
				ret.push_back(Cell.m_AminoAcid);
				break;
			}
		}
		else
		{
			ret.push_back(Cell.m_AminoAcid);
		}
	}

	return ret;
}

bool AllPeaksExplained(const auto& rgflMassPeaks, const auto& rgCells) noexcept
{
	size_t iCount = 0;

	fmt::print(fg(fmt::color::cyan), Conclude(rgCells));
	fmt::print(fg(fmt::color::cyan), " checks: ");

	for (const auto& Peak : rgflMassPeaks)
	{
		size_t v = std::find_if(rgCells.cbegin(), rgCells.cend(), std::bind(&Cell_t::Explained, std::placeholders::_1, Peak)) != rgCells.cend();
		if (!v)
			fmt::print("[{}] ", Peak.m_Value);
		iCount += v;
	}

	fmt::print(fg(fmt::color::cyan), "({}/{})\n", iCount, rgflMassPeaks.size());

	return iCount == rgflMassPeaks.size();
}

template <typename T> T FilterUnexplainablePeaks(const T& rgflMassPeaks, const auto& rgCells) noexcept
{
	T ret{};

	for (const auto& Peak : rgflMassPeaks)
	{
		if (std::find_if(rgCells.cbegin(), rgCells.cend(), std::bind(&Cell_t::Explained, std::placeholders::_1, Peak)) == rgCells.cend())
			ret.push_back(Peak);
	}

	return ret;
}

template <template <typename, typename> class T, typename AllocTy1, template <typename, typename> class U, typename AllocTy2>
void RemoveExplainedPeaks(T<MassPeak_t, AllocTy1>* prgflPendingPeaks, const U<Cell_t, AllocTy2>& rgSolutions) noexcept
{
	for (auto iter = prgflPendingPeaks->begin(); iter != prgflPendingPeaks->end(); /* Does nothing. */)
	{
		bool bExplained = false;
		for (const auto& Cell : rgSolutions)
			bExplained = bExplained || Cell.Explained(*iter);

		if (bExplained)
			iter = prgflPendingPeaks->erase(iter);
		else
			++iter;
	}
}

void fnRecursiveFromB1(list<AlternativeReality_t>* prgWorldlines, AlternativeReality_t* pThisWorldline, const function<bool(double)>& pfnShouldCheck, IgnoreLevel_e iIgnoreLv) noexcept
{
	// The counterpart ion here is actually reflecting the previous cell.
	bool bPredictedFound = false;
	auto const& TheHead = pThisWorldline->m_Solution.back();
	double flPredictedCounterpartPeak = TheHead.m_yTypeIon - g_rgAminoAcidsData.at(TheHead.m_AminoAcid).m_ResidueMass - g_rgflModMassShift[TheHead.m_Modification];

	if (pfnShouldCheck(flPredictedCounterpartPeak))
	{
		for (const auto& Peak : pThisWorldline->m_PendingPeaks)
		{
			if (Peak.Approx(flPredictedCounterpartPeak))
			{
				bPredictedFound = true;
				flPredictedCounterpartPeak = Peak.m_Value;	// Set to the observed counterpart ion.
				break;
			}
		}
		// We are not going to have the current worldline collapsed on failure check.
	}

	bool bAlreadyFoundOne = false;
	double const flFrontierPeak = TheHead.m_bTypeIon;
	AlternativeReality_t ThisCopy = *pThisWorldline;

	for (auto iter = pThisWorldline->m_PendingPeaks.rbegin(); iter != pThisWorldline->m_PendingPeaks.rend(); ++iter)
	{
		if (*iter <= flFrontierPeak)
			continue;

		if (auto const [iFound, iModType] = TestNumber(*iter - flFrontierPeak, iIgnoreLv); iFound != NOT_AN_AMINO_ACID)
		{
			if (!bAlreadyFoundOne)
			{
				bAlreadyFoundOne = true;
				pThisWorldline->m_Solution.emplace_back(iFound, iModType, iter->m_Value, flPredictedCounterpartPeak);
			}
			else
			{
				AlternativeReality_t& OtherWorld = prgWorldlines->emplace_back(ThisCopy);
				OtherWorld.m_Solution.emplace_back(iFound, iModType, iter->m_Value, flPredictedCounterpartPeak);
				OtherWorld.m_PendingPeaks.erase(std::find(OtherWorld.m_PendingPeaks.begin(), OtherWorld.m_PendingPeaks.end(), *iter));

				fnRecursiveFromB1(prgWorldlines, &OtherWorld, pfnShouldCheck, iIgnoreLv);
			}
		}
	}

	// Only clear the parsed peak after all branches expanded.
	RemoveExplainedPeaks(&pThisWorldline->m_PendingPeaks, pThisWorldline->m_Solution);

	if (bAlreadyFoundOne && !pThisWorldline->m_PendingPeaks.empty())
		fnRecursiveFromB1(prgWorldlines, pThisWorldline, pfnShouldCheck, iIgnoreLv);
	else
		pThisWorldline->m_Solution.emplace_back();	// Place a sealer dummy.
}

void fnRecursiveFromBn(list<AlternativeReality_t>* prgWorldlines, AlternativeReality_t* pThisWorldline, const function<bool(double)>& pfnShouldCheck, IgnoreLevel_e iIgnoreLv) noexcept
{
	bool bAlreadyFoundOne = false;
	double flFrontierPeak = pThisWorldline->m_Solution.front().m_bTypeIon;
	double flVarificationPeak = pThisWorldline->m_Solution[1].m_yTypeIon;
	AlternativeReality_t ThisCopy = *pThisWorldline;

	for (const auto& Peak : pThisWorldline->m_PendingPeaks)
	{
		if (Peak >= flFrontierPeak)
			continue;

		if (auto const [iFound, iModType] = TestNumber(flFrontierPeak - Peak, iIgnoreLv); iFound != NOT_AN_AMINO_ACID)
		{
			double flPredictedCounterpartPeak = flVarificationPeak + g_rgAminoAcidsData.at(iFound).m_ResidueMass + g_rgflModMassShift[iModType];
			bool bPredictedFound = false;
			if (pfnShouldCheck(flPredictedCounterpartPeak))
			{
				for (const auto& CurPeak : pThisWorldline->m_PendingPeaks)
				{
					if (CurPeak.Approx(flPredictedCounterpartPeak))
					{
						bPredictedFound = true;
						flPredictedCounterpartPeak = CurPeak.m_Value;	// Set to the observed counterpart ion.
						break;
					}
				}

				if (!bPredictedFound)
					continue;
			}

			if (!bAlreadyFoundOne)
			{
				bAlreadyFoundOne = true;
				pThisWorldline->m_Solution.front().m_AminoAcid = iFound;
				pThisWorldline->m_Solution.front().m_Modification = iModType;
				pThisWorldline->m_Solution.front().m_yTypeIon = flPredictedCounterpartPeak;
				pThisWorldline->m_Solution.emplace_front(NOT_AN_AMINO_ACID, NO_MODIFICATION, Peak.m_Value);
			}
			else
			{
				AlternativeReality_t& OtherWorld = prgWorldlines->emplace_back(ThisCopy);
				OtherWorld.m_Solution.front().m_AminoAcid = iFound;
				OtherWorld.m_Solution.front().m_Modification = iModType;
				OtherWorld.m_Solution.front().m_yTypeIon = flPredictedCounterpartPeak;
				OtherWorld.m_Solution.emplace_front(NOT_AN_AMINO_ACID, NO_MODIFICATION, Peak.m_Value);
				OtherWorld.m_PendingPeaks.erase(std::find(OtherWorld.m_PendingPeaks.begin(), OtherWorld.m_PendingPeaks.end(), Peak));

				fnRecursiveFromBn(prgWorldlines, &OtherWorld, pfnShouldCheck, iIgnoreLv);
			}
		}
	}

	// Only clear the parsed peak after all branches expanded.
	RemoveExplainedPeaks(&pThisWorldline->m_PendingPeaks, pThisWorldline->m_Solution);

	if (bAlreadyFoundOne && !pThisWorldline->m_PendingPeaks.empty())
		fnRecursiveFromBn(prgWorldlines, pThisWorldline, pfnShouldCheck, iIgnoreLv);
}

void fnRecursiveFromY1(list<AlternativeReality_t>* prgWorldlines, AlternativeReality_t* pThisWorldline, const function<bool(double)>& pfnShouldCheck, IgnoreLevel_e iIgnoreLv) noexcept
{
	// The counterpart ion here is actually reaflecting the previous cell.
	bool bPredictedFound = false;
	auto const& TheHead = pThisWorldline->m_Solution[0];
	double flPredictedCounterpartPeak = TheHead.m_bTypeIon - g_rgAminoAcidsData.at(TheHead.m_AminoAcid).m_ResidueMass - g_rgflModMassShift[TheHead.m_Modification];
	
	if (pfnShouldCheck(flPredictedCounterpartPeak))
	{
		for (const auto& Peak : pThisWorldline->m_PendingPeaks)
		{
			if (Peak.Approx(flPredictedCounterpartPeak))
			{
				bPredictedFound = true;
				flPredictedCounterpartPeak = Peak.m_Value;	// Set to the observed counterpart ion.
				break;
			}
		}
		// We are not going to have the current worldline collapsed on failure check.
	}

	bool bAlreadyFoundOne = false;
	double flFrontierPeak = TheHead.m_yTypeIon;
	AlternativeReality_t ThisCopy = *pThisWorldline;

	for (auto iter = pThisWorldline->m_PendingPeaks.rbegin(); iter != pThisWorldline->m_PendingPeaks.rend(); ++iter)
	{
		if (*iter <= flFrontierPeak)
			continue;

		if (auto const [iFound, iModType] = TestNumber(*iter - flFrontierPeak, iIgnoreLv); iFound != NOT_AN_AMINO_ACID)
		{
			if (!bAlreadyFoundOne)
			{
				bAlreadyFoundOne = true;
				pThisWorldline->m_Solution.emplace_front(iFound, iModType, flPredictedCounterpartPeak, iter->m_Value);
			}
			else
			{
				AlternativeReality_t& OtherWorld = prgWorldlines->emplace_back(ThisCopy);
				OtherWorld.m_Solution.emplace_front(iFound, iModType, flPredictedCounterpartPeak, iter->m_Value);
				OtherWorld.m_PendingPeaks.erase(std::find(OtherWorld.m_PendingPeaks.begin(), OtherWorld.m_PendingPeaks.end(), *iter));

				fnRecursiveFromY1(prgWorldlines, &OtherWorld, pfnShouldCheck, iIgnoreLv);
			}
		}
	}

	// Only clear the parsed peak after all branches expanded.
	RemoveExplainedPeaks(&pThisWorldline->m_PendingPeaks, pThisWorldline->m_Solution);

	if (bAlreadyFoundOne && !pThisWorldline->m_PendingPeaks.empty())
		fnRecursiveFromY1(prgWorldlines, pThisWorldline, pfnShouldCheck, iIgnoreLv);
	else
		pThisWorldline->m_Solution.emplace_front();	// Place a sealer dummy.
}

void fnRecursiveFromYn(list<AlternativeReality_t>* prgWorldlines, AlternativeReality_t* pThisWorldline, const function<bool(double)>& pfnShouldCheck, IgnoreLevel_e iIgnoreLv) noexcept
{
	bool bAlreadyFoundOne = false;
	double flFrontierPeak = pThisWorldline->m_Solution.back().m_yTypeIon;
	double flVarificationPeak = pThisWorldline->m_Solution[pThisWorldline->m_Solution.size() - 2U].m_bTypeIon;
	AlternativeReality_t ThisCopy = *pThisWorldline;

	for (const auto& Peak : pThisWorldline->m_PendingPeaks)
	{
		if (Peak >= flFrontierPeak)
			continue;

		if (auto const [iFound, iModType] = TestNumber(flFrontierPeak - Peak, iIgnoreLv); iFound != NOT_AN_AMINO_ACID)
		{
			double flPredictedCounterpartPeak = flVarificationPeak + g_rgAminoAcidsData.at(iFound).m_ResidueMass + g_rgflModMassShift[iModType];
			bool bPredictedFound = false;
			if (pfnShouldCheck(flPredictedCounterpartPeak))
			{
				for (const auto& CurPeak : pThisWorldline->m_PendingPeaks)
				{
					if (CurPeak.Approx(flPredictedCounterpartPeak))
					{
						bPredictedFound = true;
						flPredictedCounterpartPeak = CurPeak.m_Value;	// Set to the observed counterpart ion.
						break;
					}
				}

				if (!bPredictedFound)
					continue;
			}

			if (!bAlreadyFoundOne)
			{
				bAlreadyFoundOne = true;
				pThisWorldline->m_Solution.back().m_AminoAcid = iFound;
				pThisWorldline->m_Solution.back().m_Modification = iModType;
				pThisWorldline->m_Solution.back().m_bTypeIon = flPredictedCounterpartPeak;
				pThisWorldline->m_Solution.emplace_back(NOT_AN_AMINO_ACID, NO_MODIFICATION, 0, Peak.m_Value);
			}
			else
			{
				AlternativeReality_t& OtherWorld = prgWorldlines->emplace_back(ThisCopy);
				OtherWorld.m_Solution.back().m_AminoAcid = iFound;
				OtherWorld.m_Solution.back().m_Modification = iModType;
				OtherWorld.m_Solution.back().m_bTypeIon = flPredictedCounterpartPeak;
				OtherWorld.m_Solution.emplace_back(NOT_AN_AMINO_ACID, NO_MODIFICATION, 0, Peak.m_Value);
				OtherWorld.m_PendingPeaks.erase(std::find(OtherWorld.m_PendingPeaks.begin(), OtherWorld.m_PendingPeaks.end(), Peak));

				fnRecursiveFromYn(prgWorldlines, &OtherWorld, pfnShouldCheck, iIgnoreLv);
			}
		}
	}

	// Only clear the parsed peak after all branches expanded.
	RemoveExplainedPeaks(&pThisWorldline->m_PendingPeaks, pThisWorldline->m_Solution);

	if (bAlreadyFoundOne && !pThisWorldline->m_PendingPeaks.empty())
		fnRecursiveFromYn(prgWorldlines, pThisWorldline, pfnShouldCheck, iIgnoreLv);
}

void PrintWorldline(const AlternativeReality_t& Worldline, const auto& rgflOriginalMassData, unsigned iDecimal) noexcept
{
	auto const itBegin = rgflOriginalMassData.cbegin(), itEnd = rgflOriginalMassData.cend();

	for (const auto& Cell : Worldline.m_Solution)
	{
		if (std::find_if(itBegin, itEnd, std::bind(&MassPeak_t::Approx, std::placeholders::_1, Cell.m_bTypeIon)) != itEnd)	// #CPP23_UPGRADE
			fmt::print("{:.{}f}\t", Cell.m_bTypeIon, iDecimal);
		else
			fmt::print(fg(fmt::color::red), "{:.{}f}\t", Cell.m_bTypeIon, iDecimal);
	}

	fmt::print("\n");

	for (const auto& Cell : Worldline.m_Solution)
	{
		bool bbFound = std::find_if(itBegin, itEnd, std::bind(&MassPeak_t::Approx, std::placeholders::_1, Cell.m_bTypeIon)) != itEnd;	// #UPDATE_AT_CPP23 std::bind_front
		bool byFound = std::find_if(itBegin, itEnd, std::bind(&MassPeak_t::Approx, std::placeholders::_1, Cell.m_yTypeIon)) != itEnd;

		if (bbFound && byFound)
			fmt::print(fg(fmt::color::green), "{}{}\t", (std::underlying_type_t<AminoAcids_e>)Cell.m_AminoAcid, g_rgpszModName[Cell.m_Modification]);
		else if (bbFound || byFound)
			fmt::print(fg(fmt::color::dark_golden_rod), "{}{}\t", (std::underlying_type_t<AminoAcids_e>)Cell.m_AminoAcid, g_rgpszModName[Cell.m_Modification]);
		else
			fmt::print(fg(fmt::color::dark_red), "{}{}\t", (std::underlying_type_t<AminoAcids_e>)Cell.m_AminoAcid, g_rgpszModName[Cell.m_Modification]);
	}

	fmt::print("\n");

	for (const auto& Cell : Worldline.m_Solution)
	{
		if (std::find_if(itBegin, itEnd, std::bind(&MassPeak_t::Approx, std::placeholders::_1, Cell.m_yTypeIon)) != itEnd)
			fmt::print("{:.{}f}\t", Cell.m_yTypeIon, iDecimal);
		else
			fmt::print(fg(fmt::color::red), "{:.{}f}\t", Cell.m_yTypeIon, iDecimal);
	}

	fmt::print(fg(fmt::color::dim_gray), "\nLeft peaks({}/{}): ", Worldline.m_PendingPeaks.size(), rgflOriginalMassData.size());
	for (const auto& Peak : Worldline.m_PendingPeaks)
	{
		fmt::print(fg(fmt::color::dim_gray), "{}, ", Peak.m_Value);
	}

	fmt::print("\n");
}

export [[nodiscard]]
list<AlternativeReality_t> Solve(const vector<MassPeak_t>& rgflMassData, double M_Plus_1, unsigned iDecimal, IgnoreLevel_e iIgnoreLv) noexcept
{
	vector<MassPeak_t> rgflMassData2 = rgflMassData;
	int iMPlusOne = (int)std::round(M_Plus_1);
	int iMPlusTwo = (int)std::round(M2ZConversion<1, 2>(M_Plus_1));
	auto const fnCheck = [&rgflMassData2](double flTest) -> bool { return flTest >= rgflMassData2.back() && flTest <= rgflMassData2.front(); };

#pragma region Mass peak preprocessing.
	for (auto iter = rgflMassData2.begin(); iter != rgflMassData2.end(); /* Does nothing. */)
	{
		if (int iMass = (int)std::round(*iter); iMass == iMPlusOne || iMass == iMPlusTwo)
		{
			fmt::print("{} is the full ion peak.\n", iter->m_Value);
			iter = rgflMassData2.erase(iter);
			continue;
		}

		++iter;
	}
	std::sort(rgflMassData2.begin(), rgflMassData2.end(), std::greater<MassPeak_t>());
#pragma endregion Mass peak preprocessing.

#pragma region C-Terminal attempts
	list<AlternativeReality_t> rgCTerminalGuesses{ AlternativeReality_t{} };
	AlternativeReality_t& FirstCWorldline = rgCTerminalGuesses.front();
	FirstCWorldline.m_PendingPeaks = rgflMassData2;
	FirstCWorldline.m_MPlusOne = M_Plus_1;

	for (auto iter = FirstCWorldline.m_PendingPeaks.begin(); iter != FirstCWorldline.m_PendingPeaks.end(); /* Does nothing. */)
	{
		if (int iNeutralLoss = (int)std::round(M_Plus_1 - *iter); iNeutralLoss == 146 || iNeutralLoss == 174)
		{
			fmt::print("{} is the b[n-1] ion and the C-terminal amino acid is {}.\n", iter->m_Value, g_rgAminoAcidsData.at(iNeutralLoss == 146 ? Lysine : Arginine).m_3Letters);
			FirstCWorldline.m_Solution.emplace_front(iNeutralLoss == 146 ? Lysine : Arginine, NO_MODIFICATION, M_Plus_1 - Monoisotopic::MWt<"H2O"> /* b[n] is [M+1] ion with a H2O loss. */, iNeutralLoss + 1);
			FirstCWorldline.m_Solution.emplace_front(NOT_AN_AMINO_ACID, NO_MODIFICATION, iter->m_Value /* b[n-1] */);
			iter = FirstCWorldline.m_PendingPeaks.erase(iter);
			continue;
		}

		++iter;
	}

	if (FirstCWorldline.m_Solution.size() != 2)
	{
		for (auto iter = FirstCWorldline.m_PendingPeaks.begin(); iter != FirstCWorldline.m_PendingPeaks.end(); /* Does nothing. */)
		{
			if (int iMass = (int)std::round(*iter); iMass == 147 || iMass == 175)
			{
				fmt::print("{} is the y1 ion and the C-terminal amino acid is {}.\n", iter->m_Value, g_rgAminoAcidsData.at(iMass == 147 ? Lysine : Arginine).m_3Letters);
				FirstCWorldline.m_Solution.emplace_back(iMass == 147 ? Lysine : Arginine, NO_MODIFICATION, M_Plus_1 - Monoisotopic::MWt<"H2O"> /* b[n] is [M+1] ion with a H2O loss. */, *iter);
				iter = FirstCWorldline.m_PendingPeaks.erase(iter);
			}
			else
				++iter;
		}

		fnRecursiveFromY1(&rgCTerminalGuesses, &FirstCWorldline, fnCheck, iIgnoreLv);
	}
	else
		fnRecursiveFromBn(&rgCTerminalGuesses, &FirstCWorldline, fnCheck, iIgnoreLv);
#pragma endregion C-Terminal attempts

#pragma region N-Terminal attempts
	list<AlternativeReality_t> rgNTerminalGuesses {};
	for (auto iter = rgflMassData2.cbegin(); iter != rgflMassData2.cend(); ++iter)
	{
		if (auto const [iFound, iModType] = TestNumber(M_Plus_1 - *iter, iIgnoreLv); iFound != NOT_AN_AMINO_ACID)
		{
			fmt::print("{} could be the y[n-1] ion and N-terminal amino acid could be {}.\n", iter->m_Value, g_rgAminoAcidsData.at(iFound).m_3Letters);

			AlternativeReality_t& AnotherNWorldline = rgNTerminalGuesses.emplace_back();
			AnotherNWorldline.m_PendingPeaks = rgflMassData2;
			AnotherNWorldline.m_MPlusOne = M_Plus_1;

			AnotherNWorldline.m_Solution.emplace_back(iFound, iModType, g_rgAminoAcidsData.at(iFound).m_ResidueMass + Monoisotopic::Hydrogen, M_Plus_1 /* y[n] is [M+1] ion. */);
			AnotherNWorldline.m_Solution.emplace_back(NOT_AN_AMINO_ACID, iModType, 0, iter->m_Value /* y[n-1] */);
			AnotherNWorldline.m_PendingPeaks.erase(std::find(AnotherNWorldline.m_PendingPeaks.begin(), AnotherNWorldline.m_PendingPeaks.end(), *iter));

			fnRecursiveFromYn(&rgNTerminalGuesses, &AnotherNWorldline, fnCheck, iIgnoreLv);
		}
	}

	if (rgNTerminalGuesses.empty())
	{
		for (auto iter = rgflMassData2.cbegin(); iter != rgflMassData2.cend(); ++iter)
		{
			if (auto const [iFound, iModType] = TestNumber(*iter - Monoisotopic::Hydrogen, iIgnoreLv); iFound != NOT_AN_AMINO_ACID)
			{
				fmt::print("{} could be the b1 ion and N-terminal amino acid could be {}.\n", iter->m_Value, g_rgAminoAcidsData.at(iFound).m_3Letters);

				AlternativeReality_t& AnotherNWorldline = rgNTerminalGuesses.emplace_back();
				AnotherNWorldline.m_PendingPeaks = rgflMassData2;
				AnotherNWorldline.m_MPlusOne = M_Plus_1;

				AnotherNWorldline.m_Solution.emplace_back(iFound, iModType, iter->m_Value, M_Plus_1 /* y[n] is [M+1] ion. */);
				AnotherNWorldline.m_PendingPeaks.erase(std::find(AnotherNWorldline.m_PendingPeaks.begin(), AnotherNWorldline.m_PendingPeaks.end(), *iter));

				fnRecursiveFromB1(&rgNTerminalGuesses, &AnotherNWorldline, fnCheck, iIgnoreLv);
			}
		}
	}
#pragma endregion N-Terminal attempts

#pragma region Merging results
	list<AlternativeReality_t> rgExplanations{};

	for (const auto& Worldline : rgNTerminalGuesses)
	{
		auto szSequence = Conclude(Worldline.m_Solution);
		for (const auto& WorldlineCompareWith : rgCTerminalGuesses)
		{
			auto szSeqCompareWith = Conclude(WorldlineCompareWith.m_Solution);
			auto pszSeqCompW = szSeqCompareWith.c_str();
			auto iLength = szSequence.length();
			for (auto pszSeq = szSequence.c_str(); *pszSeq != '\0' && iLength > 1; ++pszSeq, --iLength)
			{
				if (!strncmp(pszSeq, pszSeqCompW, iLength))
				{
					auto& MergedWorldline = rgExplanations.emplace_back();
					MergedWorldline.m_MPlusOne = M_Plus_1;

					auto& Solution = MergedWorldline.m_Solution;
					Solution.insert(Solution.end(), Worldline.m_Solution.begin(), Worldline.m_Solution.begin() + (szSequence.length() - iLength));
					Solution.insert(Solution.end(), WorldlineCompareWith.m_Solution.begin() + 1, WorldlineCompareWith.m_Solution.end());

					if (std::round(CalcMWtByCells(MergedWorldline.m_Solution) + Monoisotopic::Hydrogen - M_Plus_1))
					{
						rgExplanations.pop_back();
						continue;
					}

					for (auto& Cell : Solution)
					{
						Cell.Disambiguate(rgflMassData2);
						Cell.Deduce();
					}

					MergedWorldline.m_PendingPeaks = FilterUnexplainablePeaks(rgflMassData2, Solution);
				}
			}
		}
	}

	rgExplanations.sort([](const AlternativeReality_t& lhs, const AlternativeReality_t& rhs) -> bool { return lhs.m_PendingPeaks.size() < rhs.m_PendingPeaks.size(); });
#pragma endregion Merging results

	if (rgExplanations.empty())
	{
		fmt::print(fg(fmt::color::dark_golden_rod), "No sensible solution successfully combined.\n");

		if (!rgNTerminalGuesses.empty())
		{
			fmt::print(fg(fmt::color::cyan), "Listing N-terminal potential solution(s).\n");
			for (const auto& Worldline : rgNTerminalGuesses)
			{
				PrintWorldline(Worldline, rgflMassData2, iDecimal);
				fmt::print("\n\n");
			}
		}
		else
			fmt::print(fg(fmt::color::dark_red), "No N-terminal potential solution.\n");

		if (!rgCTerminalGuesses.empty())
		{
			fmt::print(fg(fmt::color::magenta), "Listing C-terminal potential solution(s).\n");
			for (const auto& Worldline : rgCTerminalGuesses)
			{
				PrintWorldline(Worldline, rgflMassData2, iDecimal);
				fmt::print("\n\n");
			}
		}
		else
			fmt::print(fg(fmt::color::dark_red), "No C-terminal potential solution.\n");
	}
	else
	{
		for (const auto& Explanation : rgExplanations)
		{
			PrintWorldline(Explanation, rgflMassData2, iDecimal);
			fmt::print("\n\n");
		}
	}

	return rgExplanations;
}

export void MarkPeaks(const AlternativeReality_t& Worldline, auto* prgflMassPeaks)
{
	short iCount = 1;
	short const iTotalCount = (short)Worldline.m_Solution.size() + 1;
	for (const auto& Cell : Worldline.m_Solution)
	{
		for (auto& Peak : *prgflMassPeaks)
		{
			// a
			if (Cell.m_aTypeIon && !gcem::round(Peak - Cell.m_aTypeIon))
				Peak.Identify(fmt::format("a{}", iCount));

			// b
			else if (Cell.m_bNeutralLostIon && !gcem::round(Peak - Cell.m_bNeutralLostIon))
				Peak.Identify(fmt::format("b{}\'", iCount));
			else if (Cell.m_bAsteriskIon && !gcem::round(Peak - Cell.m_bAsteriskIon))
				Peak.Identify(fmt::format("b{}*", iCount));
			else if (Cell.m_bCircleIon && !gcem::round(Peak - Cell.m_bCircleIon))
				//Peak.Identify(fmt::vformat((const char*)u8"b{}°", fmt::make_format_args(iCount)));	// Fuck C++20
				Peak.Identify(fmt::format(u8"b{}°", iCount));	// Fuck C++20
			else if (Cell.m_bTypeIon && !gcem::round(Peak - Cell.m_bTypeIon))
				Peak.Identify(fmt::format("b{}", iCount));

			// y
			else if (Cell.m_yNeutralLostIon && !gcem::round(Peak - Cell.m_yNeutralLostIon))
				Peak.Identify(fmt::format("y{}\'", iTotalCount - iCount));
			else if (Cell.m_yAsteriskIon && !gcem::round(Peak - Cell.m_yAsteriskIon))
				Peak.Identify(fmt::format("y{}*", iTotalCount - iCount));
			else if (Cell.m_yCircleIon && !gcem::round(Peak - Cell.m_yCircleIon))
				Peak.Identify(fmt::format(u8"y{}°", iTotalCount - iCount));
			else if (Cell.m_yTypeIon && !gcem::round(Peak - Cell.m_yTypeIon))
				Peak.Identify(fmt::format("y{}", iTotalCount - iCount));
		}

		++iCount;
	}

	double const MPlusTwo = M2ZConversion<1, 2>(Worldline.m_MPlusOne);
	double const MPlusThree = M2ZConversion<1, 3>(Worldline.m_MPlusOne);
	for (auto& Peak : *prgflMassPeaks)
	{
		if (Peak.m_String.empty())
		{
			if (!gcem::round(Peak - Worldline.m_MPlusOne))
				Peak.Identify("[M+H]");
			else if (!gcem::round(Peak - MPlusTwo))
				Peak.Identify("[M+2H]");
			else if (!gcem::round(Peak - MPlusThree))	[[unlikely]]
				Peak.Identify("[M+3H]");
		}
	}
}

export void ResetPeaks(auto* prgflMassPeaks)
{
	for (auto& Peak : *prgflMassPeaks)
		Peak.m_String.clear();
}
