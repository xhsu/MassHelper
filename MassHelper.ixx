module;

// C
#include <cassert>
#include <cfloat>

// C++
#include <algorithm>
#include <array>
#include <concepts>
#include <deque>
#include <format>
#include <functional>
#include <iostream>
#include <tuple>
#include <unordered_map>
#include <vector>

// C++ Exp
#include <experimental/generator>

// Others
#include <gcem.hpp>

// Platforms
#include <fcntl.h>
#include <io.h>

export module MassHelper;

import UtlConcepts;
import UtlWinConsole;

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
using std::wstring;

using std::experimental::generator;

export enum AminoAcids_e : wchar_t
{
	NOT_AN_AMINO_ACID = L'\0',

	Glycine = L'G',
	Alanine = L'A',
	Serine = L'S',
	Proline = L'P',
	Valine = L'V',
	Threonine = L'T',
	Cysteine = L'C',
	Isoleucine = L'I',
	Leucine = L'L',
	Asparagine = L'N',
	Aspartic_Acid = L'D',
	Lysine = L'K',
	Glutamine = L'Q',
	Glutamic_Acid = L'E',
	Methionine = L'M',
	Histidine = L'H',
	Methionine_Sulfoxide = L'm',
	Phenylalanine = L'F',
	Arginine = L'R',
	Carbamidomethyl_Cysteine = L'ç',
	Carboxyethyl_Cysteine = L'c',
	Tyrosine = L'Y',
	Tryptophan = L'W',
};

enum ProteinModifications_e : unsigned char
{
	NO_MODIFICATION = 0,
	GLU_CYCLIZATION,
	METHYLATION,
	ACETYLATION,
	PHOSPHORYLATION,
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
	L"",
	L"(cy.)",
	L"(met.)",
	L"(ace.)",
	L"(pho.)",
};

struct Molecule_t
{
	double m_Mass = 0;
	const wchar_t* m_Name = L"none";

	template <StringLiteral FORMULA>
	static consteval decltype(auto) Monoisotope(void) noexcept
	{
		return Molecule_t{ Monoisotopic::MWt<FORMULA>, FORMULA };
	}
};

struct AminoAcid_t
{
	double m_ResidueMass = 0;
	int m_ImmoniumIonMass = 0;
	const wchar_t* m_3Letters = L"err";
	Molecule_t m_CharacteristicNeuLoss{};
};

constexpr auto g_rgCommonNeuLoss = array
{
	Molecule_t::Monoisotope<L"H2O">(),
	Molecule_t::Monoisotope<L"NH3">(),
	Molecule_t::Monoisotope<L"HPO3">(),
	Molecule_t::Monoisotope<L"H3PO4">(),
};

constexpr Molecule_t NOTHINGNESS = {};
constexpr Molecule_t WATER = Molecule_t::Monoisotope<L"H2O">();
constexpr Molecule_t AMMONIA = Molecule_t::Monoisotope<L"NH3">();

const unordered_map<std::underlying_type_t<AminoAcids_e>, AminoAcid_t> g_rgAminoAcidsData =
{
	{AminoAcids_e::NOT_AN_AMINO_ACID,		{}},
	{AminoAcids_e::Glycine,					{57.0215, 30, L"Gly", NOTHINGNESS}},
	{AminoAcids_e::Alanine,					{71.0371, 44, L"Ala", NOTHINGNESS}},
	{AminoAcids_e::Serine,					{87.0320, 60, L"Ser", WATER}},
	{AminoAcids_e::Proline,					{97.0528, 70, L"Pro", NOTHINGNESS}},
	{AminoAcids_e::Valine,					{99.0684, 72, L"Val", NOTHINGNESS}},
	{AminoAcids_e::Threonine,				{101.0477, 74, L"Thr", WATER}},
	{AminoAcids_e::Cysteine,				{103.0092, 76, L"Cys", Molecule_t::Monoisotope<L"HS2">()}},
	{AminoAcids_e::Isoleucine,				{113.0841, 86, L"Ile", NOTHINGNESS}},
	{AminoAcids_e::Leucine,					{113.0841, 86, L"Leu", NOTHINGNESS}},
	{AminoAcids_e::Asparagine,				{114.0429, 87, L"Asn", AMMONIA}},
	{AminoAcids_e::Aspartic_Acid,			{115.0269, 88, L"Asp", WATER}},
	{AminoAcids_e::Lysine,					{128.0950, 101, L"Lys", AMMONIA}},
	{AminoAcids_e::Glutamine,				{128.0586, 101, L"Gln", AMMONIA}},
	{AminoAcids_e::Glutamic_Acid,			{129.0426, 102, L"Glu", WATER}},
	{AminoAcids_e::Methionine,				{131.0405, 104, L"Met", Molecule_t::Monoisotope<L"CH3SH">()}},
	{AminoAcids_e::Histidine,				{137.0589, 110, L"His", NOTHINGNESS}},
	{AminoAcids_e::Methionine_Sulfoxide,	{147.0354, 120, L"Mox", Molecule_t::Monoisotope<L"CH3SOH">()}},
	{AminoAcids_e::Phenylalanine,			{147.0684, 120, L"Phe", NOTHINGNESS}},
	{AminoAcids_e::Arginine,				{156.0684, 129, L"Arg", AMMONIA}},
	{AminoAcids_e::Carbamidomethyl_Cysteine,{160.0307, 133, L"Cac", Molecule_t::Monoisotope<L"HSCH2CONH2">()}},
	{AminoAcids_e::Carboxyethyl_Cysteine,	{161.0147, 134, L"Cec", Molecule_t::Monoisotope<L"HSCH2COOH">()}},
	{AminoAcids_e::Tyrosine,				{163.0633, 134, L"Tyr", NOTHINGNESS}},
	{AminoAcids_e::Tryptophan,				{186.0793, 159, L"Trp", NOTHINGNESS}},
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
			cout_r() << std::format("Ion {} was previously identified as {}.\n", sz, m_String) << white_text;

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

pair<AminoAcids_e, ProteinModifications_e> TestNumber(double flPeakDiff, double const flTolerance = 0.5) noexcept
{
	using Ty = tuple<decltype(g_rgAminoAcidsData)::key_type, ProteinModifications_e, double>;
	auto const fn = [&]() -> generator<Ty>
	{
		for (const auto& [AminoAcid, Data] : g_rgAminoAcidsData)
		{
			// Exact match on residue
			if (double const diff = std::abs(flPeakDiff - Data.m_ResidueMass); diff < flTolerance)
				co_yield std::make_tuple(AminoAcid, NO_MODIFICATION, diff);

			// It doesn't make any sense to have further modifications.
			if (AminoAcid == Methionine_Sulfoxide || AminoAcid == Carboxyethyl_Cysteine || AminoAcid == Carbamidomethyl_Cysteine)	// It doesn't make any sense to have further modifications.
				continue;

			// Common amino acid modifications
			// Skip the first possibility: Pyroglutamic acid formed from N-terminal Gln
			// Parent mass shift (Da): -17
			for (int i = METHYLATION; i < g_rgflModMassShift.size(); ++i)
				if (auto diff = std::abs(flPeakDiff - (Data.m_ResidueMass + g_rgflModMassShift[i])); diff < flTolerance)
					co_yield std::make_tuple(AminoAcid, (ProteinModifications_e)i, diff);
		}

		co_return;
	};

	auto const fnCompLesser = [](const Ty& lhs, const Ty& rhs) -> bool
	{
		auto const& [lty, lmod, ldiff] = lhs;
		auto const& [rty, rmod, rdiff] = rhs;

		if (std::abs(ldiff - rdiff) < 0.01 && lmod != rmod)
			return lmod == NO_MODIFICATION;

		return ldiff < rdiff;
	};

	std::vector<Ty> ret(fn().begin(), fn().end());
	if (auto const it = std::min_element(ret.cbegin(), ret.cend(), fnCompLesser); it != ret.cend())
		return std::make_pair((AminoAcids_e)std::get<0>(*it), std::get<1>(*it));

	return std::make_pair(NOT_AN_AMINO_ACID, NO_MODIFICATION);
}

export list<pair<wstring, double>> SummaryPeakAsList(double flPeakDiff, double const flTolerance = 0.5) noexcept
{
	using ListTy = list<pair<wstring, double>>;

	if (flPeakDiff <= DBL_EPSILON)
		return { { L"SELF", 0.0 } };

	ListTy rgFound{};

	for (const auto& [AminoAcid, Data] : g_rgAminoAcidsData)
	{
		if (AminoAcid == NOT_AN_AMINO_ACID)
			continue;

		if (auto diff = flPeakDiff - Data.m_ResidueMass; std::abs(diff) < flTolerance)
			rgFound.emplace_back(Data.m_3Letters + std::format(L"({})", AminoAcid), diff);

		if (auto diff = flPeakDiff - Data.m_CharacteristicNeuLoss.m_Mass; Data.m_CharacteristicNeuLoss.m_Mass > 0 && std::abs(diff) < flTolerance)
			rgFound.emplace_back(Data.m_CharacteristicNeuLoss.m_Name, diff);

		// It doesn't make any sense to have further modifications.
		if (AminoAcid == Methionine_Sulfoxide || AminoAcid == Carboxyethyl_Cysteine || AminoAcid == Carbamidomethyl_Cysteine)
			continue;

		// Common amino acid modifications
		// Skip the first possibility: Pyroglutamic acid formed from N-terminal Gln
		// Parent mass shift (Da): -17

		for (int i = METHYLATION; i < g_rgflModMassShift.size(); ++i)
			if (auto diff = flPeakDiff - (Data.m_ResidueMass + g_rgflModMassShift[i]); std::abs(diff) < flTolerance)
				rgFound.emplace_back(Data.m_3Letters + std::format(L"({}){}", AminoAcid, g_rgpszModName[i]), diff);
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

export wstring SummaryPeakAsString(double flPeakDiff, bool const bOneLine = true) noexcept
{
	if (flPeakDiff <= DBL_EPSILON)
		return L"SELF";

	if (auto const rgPossibilities = SummaryPeakAsList(flPeakDiff); !rgPossibilities.empty())
	{
		wstring sz{};
		for (const auto& [szInfo, flDelta] : rgPossibilities)
			sz += std::format(L"{}[{:.4f}]{}", szInfo, flDelta, bOneLine ? L", " : L" \n");

		sz.pop_back();
		sz.pop_back();
		return sz;
	}
	else
		return L"-";
}

export double CalcMWtBySeq(const wstring& szSeq) noexcept
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
wstring Conclude(const auto& rgCells) noexcept
{
	wstring ret;
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

	cout_cyan() << Conclude(rgCells) << " checks: ";

	for (const auto& Peak : rgflMassPeaks)
	{
		size_t v = std::find_if(rgCells.cbegin(), rgCells.cend(), std::bind(&Cell_t::Explained, std::placeholders::_1, Peak)) != rgCells.cend();
		if (!v)
			cout_cyan() << std::format("[{}] ", Peak.m_Value);
		iCount += v;
	}

	cout_cyan() << '(' << iCount << '/' << rgflMassPeaks.size() << ")\n";

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

void fnRecursiveFromB1(list<AlternativeReality_t>* prgWorldlines, AlternativeReality_t* pThisWorldline, const function<bool(double)>& pfnShouldCheck) noexcept
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

		if (auto const [iFound, iModType] = TestNumber(*iter - flFrontierPeak); iFound != NOT_AN_AMINO_ACID)
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

				fnRecursiveFromB1(prgWorldlines, &OtherWorld, pfnShouldCheck);
			}
		}
	}

	// Only clear the parsed peak after all branches expanded.
	RemoveExplainedPeaks(&pThisWorldline->m_PendingPeaks, pThisWorldline->m_Solution);

	if (bAlreadyFoundOne && !pThisWorldline->m_PendingPeaks.empty())
		fnRecursiveFromB1(prgWorldlines, pThisWorldline, pfnShouldCheck);
	else
		pThisWorldline->m_Solution.emplace_back();	// Place a sealer dummy.
}

void fnRecursiveFromBn(list<AlternativeReality_t>* prgWorldlines, AlternativeReality_t* pThisWorldline, const function<bool(double)>& pfnShouldCheck) noexcept
{
	bool bAlreadyFoundOne = false;
	double flFrontierPeak = pThisWorldline->m_Solution.front().m_bTypeIon;
	double flVarificationPeak = pThisWorldline->m_Solution[1].m_yTypeIon;
	AlternativeReality_t ThisCopy = *pThisWorldline;

	for (const auto& Peak : pThisWorldline->m_PendingPeaks)
	{
		if (Peak >= flFrontierPeak)
			continue;

		if (auto const [iFound, iModType] = TestNumber(flFrontierPeak - Peak); iFound != NOT_AN_AMINO_ACID)
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

				fnRecursiveFromBn(prgWorldlines, &OtherWorld, pfnShouldCheck);
			}
		}
	}

	// Only clear the parsed peak after all branches expanded.
	RemoveExplainedPeaks(&pThisWorldline->m_PendingPeaks, pThisWorldline->m_Solution);

	if (bAlreadyFoundOne && !pThisWorldline->m_PendingPeaks.empty())
		fnRecursiveFromBn(prgWorldlines, pThisWorldline, pfnShouldCheck);
}

void fnRecursiveFromY1(list<AlternativeReality_t>* prgWorldlines, AlternativeReality_t* pThisWorldline, const function<bool(double)>& pfnShouldCheck) noexcept
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

		if (auto const [iFound, iModType] = TestNumber(*iter - flFrontierPeak); iFound != NOT_AN_AMINO_ACID)
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

				fnRecursiveFromY1(prgWorldlines, &OtherWorld, pfnShouldCheck);
			}
		}
	}

	// Only clear the parsed peak after all branches expanded.
	RemoveExplainedPeaks(&pThisWorldline->m_PendingPeaks, pThisWorldline->m_Solution);

	if (bAlreadyFoundOne && !pThisWorldline->m_PendingPeaks.empty())
		fnRecursiveFromY1(prgWorldlines, pThisWorldline, pfnShouldCheck);
	else
		pThisWorldline->m_Solution.emplace_front();	// Place a sealer dummy.
}

void fnRecursiveFromYn(list<AlternativeReality_t>* prgWorldlines, AlternativeReality_t* pThisWorldline, const function<bool(double)>& pfnShouldCheck) noexcept
{
	bool bAlreadyFoundOne = false;
	double flFrontierPeak = pThisWorldline->m_Solution.back().m_yTypeIon;
	double flVarificationPeak = pThisWorldline->m_Solution[pThisWorldline->m_Solution.size() - 2U].m_bTypeIon;
	AlternativeReality_t ThisCopy = *pThisWorldline;

	for (const auto& Peak : pThisWorldline->m_PendingPeaks)
	{
		if (Peak >= flFrontierPeak)
			continue;

		if (auto const [iFound, iModType] = TestNumber(flFrontierPeak - Peak); iFound != NOT_AN_AMINO_ACID)
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

				fnRecursiveFromYn(prgWorldlines, &OtherWorld, pfnShouldCheck);
			}
		}
	}

	// Only clear the parsed peak after all branches expanded.
	RemoveExplainedPeaks(&pThisWorldline->m_PendingPeaks, pThisWorldline->m_Solution);

	if (bAlreadyFoundOne && !pThisWorldline->m_PendingPeaks.empty())
		fnRecursiveFromYn(prgWorldlines, pThisWorldline, pfnShouldCheck);
}

void PrintWorldline(const AlternativeReality_t& Worldline, const auto& rgflOriginalMassData) noexcept
{
	auto const itBegin = rgflOriginalMassData.cbegin(), itEnd = rgflOriginalMassData.cend();

	for (const auto& Cell : Worldline.m_Solution)
	{
		if (std::find_if(itBegin, itEnd, std::bind(&MassPeak_t::Approx, std::placeholders::_1, Cell.m_bTypeIon)) != itEnd)	// #CPP23_UPGRADE
			cout_w() << Cell.m_bTypeIon << '\t';
		else
			cout_pink() << Cell.m_bTypeIon << '\t';
	}

	cout_w() << '\n';

	for (const auto& Cell : Worldline.m_Solution)
	{
		bool bbFound = std::find_if(itBegin, itEnd, std::bind(&MassPeak_t::Approx, std::placeholders::_1, Cell.m_bTypeIon)) != itEnd;	// #UPDATE_AT_CPP23 std::bind_front
		bool byFound = std::find_if(itBegin, itEnd, std::bind(&MassPeak_t::Approx, std::placeholders::_1, Cell.m_yTypeIon)) != itEnd;

		if (bbFound && byFound)
			wcout_g() << (std::underlying_type_t<AminoAcids_e>)Cell.m_AminoAcid << g_rgpszModName[Cell.m_Modification] << '\t';
		else if (bbFound || byFound)
			wcout_gold() << (std::underlying_type_t<AminoAcids_e>)Cell.m_AminoAcid << g_rgpszModName[Cell.m_Modification] << '\t';
		else
			wcout_r() << (std::underlying_type_t<AminoAcids_e>)Cell.m_AminoAcid << g_rgpszModName[Cell.m_Modification] << '\t';
	}

	cout_w() << '\n';

	for (const auto& Cell : Worldline.m_Solution)
	{
		if (std::find_if(itBegin, itEnd, std::bind(&MassPeak_t::Approx, std::placeholders::_1, Cell.m_yTypeIon)) != itEnd)
			cout_w() << Cell.m_yTypeIon << '\t';
		else
			cout_pink() << Cell.m_yTypeIon << '\t';
	}

	cout_gray() << std::format("\nLeft peaks({}/{}): ", Worldline.m_PendingPeaks.size(), rgflOriginalMassData.size());
	for (const auto& Peak : Worldline.m_PendingPeaks)
	{
		cout_gray() << Peak.m_Value << ", ";
	}

	cout_w() << '\n';
}

export [[nodiscard]]
list<AlternativeReality_t> Solve(const vector<MassPeak_t>& rgflMassData, double M_Plus_1) noexcept
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
			std::cout << std::format("{} is the full ion peak.\n", iter->m_Value);
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
			std::wcout << std::format(L"{} is the b[n-1] ion and the C-terminal amino acid is {}.\n", iter->m_Value, g_rgAminoAcidsData.at(iNeutralLoss == 146 ? Lysine : Arginine).m_3Letters);
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
				std::wcout << std::format(L"{} is the y1 ion and the C-terminal amino acid is {}.\n", iter->m_Value, g_rgAminoAcidsData.at(iMass == 147 ? Lysine : Arginine).m_3Letters);
				FirstCWorldline.m_Solution.emplace_back(iMass == 147 ? Lysine : Arginine, NO_MODIFICATION, M_Plus_1 - Monoisotopic::MWt<"H2O"> /* b[n] is [M+1] ion with a H2O loss. */, *iter);
				iter = FirstCWorldline.m_PendingPeaks.erase(iter);
			}
			else
				++iter;
		}

		fnRecursiveFromY1(&rgCTerminalGuesses, &FirstCWorldline, fnCheck);
	}
	else
		fnRecursiveFromBn(&rgCTerminalGuesses, &FirstCWorldline, fnCheck);
#pragma endregion C-Terminal attempts

#pragma region N-Terminal attempts
	list<AlternativeReality_t> rgNTerminalGuesses {};
	for (auto iter = rgflMassData2.cbegin(); iter != rgflMassData2.cend(); ++iter)
	{
		if (auto const [iFound, iModType] = TestNumber(M_Plus_1 - *iter); iFound != NOT_AN_AMINO_ACID)
		{
			std::wcout << std::format(L"{} could be the y[n-1] ion and N-terminal amino acid could be {}.\n", iter->m_Value, g_rgAminoAcidsData.at(iFound).m_3Letters);

			AlternativeReality_t& AnotherNWorldline = rgNTerminalGuesses.emplace_back();
			AnotherNWorldline.m_PendingPeaks = rgflMassData2;
			AnotherNWorldline.m_MPlusOne = M_Plus_1;

			AnotherNWorldline.m_Solution.emplace_back(iFound, iModType, g_rgAminoAcidsData.at(iFound).m_ResidueMass + Monoisotopic::Hydrogen, M_Plus_1 /* y[n] is [M+1] ion. */);
			AnotherNWorldline.m_Solution.emplace_back(NOT_AN_AMINO_ACID, iModType, 0, iter->m_Value /* y[n-1] */);
			AnotherNWorldline.m_PendingPeaks.erase(std::find(AnotherNWorldline.m_PendingPeaks.begin(), AnotherNWorldline.m_PendingPeaks.end(), *iter));

			fnRecursiveFromYn(&rgNTerminalGuesses, &AnotherNWorldline, fnCheck);
		}
	}

	if (rgNTerminalGuesses.empty())
	{
		for (auto iter = rgflMassData2.cbegin(); iter != rgflMassData2.cend(); ++iter)
		{
			if (auto const [iFound, iModType] = TestNumber(*iter - Monoisotopic::Hydrogen); iFound != NOT_AN_AMINO_ACID)
			{
				std::wcout << std::format(L"{} could be the b1 ion and N-terminal amino acid could be {}.\n", iter->m_Value, g_rgAminoAcidsData.at(iFound).m_3Letters);

				AlternativeReality_t& AnotherNWorldline = rgNTerminalGuesses.emplace_back();
				AnotherNWorldline.m_PendingPeaks = rgflMassData2;
				AnotherNWorldline.m_MPlusOne = M_Plus_1;

				AnotherNWorldline.m_Solution.emplace_back(iFound, iModType, iter->m_Value, M_Plus_1 /* y[n] is [M+1] ion. */);
				AnotherNWorldline.m_PendingPeaks.erase(std::find(AnotherNWorldline.m_PendingPeaks.begin(), AnotherNWorldline.m_PendingPeaks.end(), *iter));

				fnRecursiveFromB1(&rgNTerminalGuesses, &AnotherNWorldline, fnCheck);
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
			for (const wchar_t* pszSeq = szSequence.c_str(); *pszSeq != '\0' && iLength > 1; ++pszSeq, --iLength)
			{
				if (!wcsncmp(pszSeq, pszSeqCompW, iLength))
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
		cout_gold() << "No sensible solution successfully combined.\n";

		if (!rgNTerminalGuesses.empty())
		{
			cout_cyan() << "Listing N-terminal potential solution(s).\n";
			for (const auto& Worldline : rgNTerminalGuesses)
			{
				PrintWorldline(Worldline, rgflMassData2);
				cout_w() << "\n\n";
			}
		}
		else
			cout_r() << "No N-terminal potential solution.\n";

		if (!rgCTerminalGuesses.empty())
		{
			cout_magenta() << "Listing C-terminal potential solution(s).\n";
			for (const auto& Worldline : rgCTerminalGuesses)
			{
				PrintWorldline(Worldline, rgflMassData2);
				cout_w() << "\n\n";
			}
		}
		else
			cout_r() << "No C-terminal potential solution.\n";
	}
	else
	{
		for (const auto& Explanation : rgExplanations)
		{
			PrintWorldline(Explanation, rgflMassData2);
			cout_w() << "\n\n";
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
				Peak.Identify(std::format("a{}", iCount));

			// b
			else if (Cell.m_bNeutralLostIon && !gcem::round(Peak - Cell.m_bNeutralLostIon))
				Peak.Identify(std::format("b{}\'", iCount));
			else if (Cell.m_bAsteriskIon && !gcem::round(Peak - Cell.m_bAsteriskIon))
				Peak.Identify(std::format("b{}*", iCount));
			else if (Cell.m_bCircleIon && !gcem::round(Peak - Cell.m_bCircleIon))
				Peak.Identify(std::vformat((const char*)u8"b{}°", std::make_format_args(iCount)));	// Fuck C++20
			else if (Cell.m_bTypeIon && !gcem::round(Peak - Cell.m_bTypeIon))
				Peak.Identify(std::format("b{}", iCount));

			// y
			else if (Cell.m_yNeutralLostIon && !gcem::round(Peak - Cell.m_yNeutralLostIon))
				Peak.Identify(std::format("y{}\'", iTotalCount - iCount));
			else if (Cell.m_yAsteriskIon && !gcem::round(Peak - Cell.m_yAsteriskIon))
				Peak.Identify(std::format("y{}*", iTotalCount - iCount));
			else if (Cell.m_yCircleIon && !gcem::round(Peak - Cell.m_yCircleIon))
				Peak.Identify(std::vformat((const char*)u8"y{}°", std::make_format_args(iTotalCount - iCount)));
			else if (Cell.m_yTypeIon && !gcem::round(Peak - Cell.m_yTypeIon))
				Peak.Identify(std::format("y{}", iTotalCount - iCount));
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
