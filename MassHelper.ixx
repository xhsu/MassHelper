module;

#include <algorithm>
#include <array>
#include <cassert>
#include <cfloat>
#include <concepts>
#include <format>
#include <functional>
#include <iostream>
#include <unordered_map>
#include <vector>

#include <gcem.hpp>

export module MassHelper;

import UtlConcepts;
import UtlWinConsole;

import PeriodicTable;

using std::array;
using std::function;
using std::list;
using std::pair;
using std::string;
using std::unordered_map;
using std::vector;

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
	Carboxyethyl_Cysteine = 'c',
	Tyrosine = 'Y',
	Tryptophan = 'W',
};

struct Molecule_t
{
	int m_Mass = 0;
	const char* m_Name = "none";
};

struct AminoAcid_t
{
	double m_ResidueMass = 0;
	int m_ImmoniumIonMass = 0;
	const char* m_3Letters = "err";
	Molecule_t m_NeutralLoss{};
};

constexpr Molecule_t NOTHINGNESS = {};
constexpr Molecule_t WATER = { 18, "H2O" };
constexpr Molecule_t AMMONIA = { 17, "NH3" };

unordered_map<char, AminoAcid_t> g_rgAminoAcidsData =
{
	{AminoAcids_e::NOT_AN_AMINO_ACID,		{}},
	{AminoAcids_e::Glycine,					{57.0215, 30, "Gly", NOTHINGNESS}},
	{AminoAcids_e::Alanine,					{71.0371, 44, "Ala", NOTHINGNESS}},
	{AminoAcids_e::Serine,					{87.0320, 60, "Ser", WATER}},
	{AminoAcids_e::Proline,					{97.0528, 70, "Pro", NOTHINGNESS}},
	{AminoAcids_e::Valine,					{99.0684, 72, "Val", NOTHINGNESS}},
	{AminoAcids_e::Threonine,				{101.0477, 74, "Thr", WATER}},
	{AminoAcids_e::Cysteine,				{103.0092, 76, "Cys", { 34, "SH2" }}},
	{AminoAcids_e::Isoleucine,				{113.0841, 86, "Ile", NOTHINGNESS}},
	{AminoAcids_e::Leucine,					{113.0841, 86, "Leu", NOTHINGNESS}},
	{AminoAcids_e::Asparagine,				{114.0429, 87, "Asn", AMMONIA}},
	{AminoAcids_e::Aspartic_Acid,			{115.0269, 88, "Asp", WATER}},
	{AminoAcids_e::Lysine,					{128.0950, 101, "Lys", AMMONIA}},
	{AminoAcids_e::Glutamine,				{128.0586, 101, "Gln", AMMONIA}},
	{AminoAcids_e::Glutamic_Acid,			{129.0426, 102, "Glu", WATER}},
	{AminoAcids_e::Methionine,				{131.0405, 104, "Met", { 48, "CH3SH" }}},
	{AminoAcids_e::Histidine,				{137.0589, 110, "His", NOTHINGNESS}},
	{AminoAcids_e::Methionine_Sulfoxide,	{147.0354, 120, "Mox", { 64, "CH3SOH" }}},
	{AminoAcids_e::Phenylalanine,			{147.0684, 120, "Phe", NOTHINGNESS}},
	{AminoAcids_e::Arginine,				{156.0684, 129, "Arg", AMMONIA}},
	{AminoAcids_e::Carboxyethyl_Cysteine,	{161.0147, 134, "Cec", { 92, "HSCH2COOH" }}},
	{AminoAcids_e::Tyrosine,				{163.0633, 134, "Tyr", NOTHINGNESS}},
	{AminoAcids_e::Tryptophan,				{186.0793, 159, "Trp", NOTHINGNESS}},
};

export enum class IonType : char
{
	UNKNOWN_TYPE = '\0',

	// Full polypeptide
	M_PLUS_H = '1',
	M_PLUS_2H = '2',
	
	// N-terminal
	a = 'a', // C-C, alpha C and carboxylic C, counterpart x
	b = 'b', // C-N, peptide bond, counterpart y
	c = 'c', // N-C, amine group from next amino acid and its alpha C, counterpart z

	// C-terminal
	x = 'x',
	y = 'y',
	z = 'z',
};

export enum IonTrait : uint16
{
	NO_LOSS = 0U,
	LOSS_H2O,
	LOSS_NH3,
	LOSS_SPEC_NEUTRAL,
};

export struct MassPeak_t
{
	constexpr MassPeak_t() noexcept {}
	constexpr MassPeak_t(Arithmetic auto v) noexcept : m_Value(v) {}

	double m_Value = 0;
	string m_String {};

	void Identify(const string& sz) noexcept
	{
		if (!m_String.empty())
			cout_r() << std::format("Ion {} was previously identified as {}.\n", sz, m_String);

		m_String = sz;
	}

	constexpr bool Approx(double v) const noexcept { return !(bool)gcem::round(m_Value - v); }

	constexpr auto operator<=> (Arithmetic auto rhs) const noexcept { return m_Value <=> rhs; }
	constexpr bool operator== (Arithmetic auto rhs) const noexcept { return gcem::abs(m_Value - rhs) <= DBL_EPSILON; }
	constexpr operator double& () noexcept { return m_Value; }
	constexpr operator const double& () const noexcept { return m_Value; }
};

struct Cell_t
{
	constexpr Cell_t(void) noexcept {}
	constexpr Cell_t(AminoAcids_e what, double b = 0, double y = 0) noexcept : m_AminoAcid(what), m_bTypeIon(b), m_yTypeIon(y) {}

	double m_aTypeIon = 0;	// Subtract -CO-
	double m_bNeutralLostIon = 0;
	double m_bAsteriskIon = 0;	// -Ammonia
	double m_bCircleIon = 0;	// -Water
	double m_bTypeIon = 0;	// Including this amino acid Coule be b[n] i.e. [M+1] ion itself.
	AminoAcids_e m_AminoAcid = NOT_AN_AMINO_ACID;
	double m_yTypeIon = 0;	// Including this amino acid. Could be y[n] i.e. [M+1] ion itself.
	double m_yCircleIon = 0;	// -Water
	double m_yAsteriskIon = 0;	// -Ammonia
	double m_yNeutralLostIon = 0;

	void Deduce(void) noexcept
	{
		auto const iNeutralMass = g_rgAminoAcidsData[m_AminoAcid].m_NeutralLoss.m_Mass;

		if (!m_aTypeIon)
			m_aTypeIon = m_bTypeIon - amu::Carbon - amu::Oxygen;
		if (!m_bNeutralLostIon && m_AminoAcid != NOT_AN_AMINO_ACID && iNeutralMass && iNeutralMass != 17 && iNeutralMass != 18)
			m_bNeutralLostIon = m_bTypeIon - g_rgAminoAcidsData[m_AminoAcid].m_NeutralLoss.m_Mass;
		if (!m_bAsteriskIon)
			m_bAsteriskIon = m_bTypeIon - amu::Hydrogen * 3 - amu::Nitrogen;
		if (!m_bCircleIon)
			m_bCircleIon = m_bTypeIon - amu::Hydrogen * 2 - amu::Oxygen;

		if (!m_yCircleIon)
			m_yCircleIon = m_yTypeIon - amu::Hydrogen * 2 - amu::Oxygen;
		if (!m_yAsteriskIon)
			m_yAsteriskIon = m_yTypeIon - amu::Hydrogen * 3 - amu::Nitrogen;
		if (!m_yNeutralLostIon && m_AminoAcid != NOT_AN_AMINO_ACID && iNeutralMass && iNeutralMass != 17 && iNeutralMass != 18)
			m_yNeutralLostIon = m_yTypeIon - g_rgAminoAcidsData[m_AminoAcid].m_NeutralLoss.m_Mass;
	}

	template<typename T> void Disambiguate(const T& rgflMassPeaks) noexcept
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
			auto const itYNeuPeak = std::find_if(itBegin, itEnd, std::bind(&MassPeak_t::Approx, std::placeholders::_1, m_yTypeIon - g_rgAminoAcidsData[Methionine_Sulfoxide].m_NeutralLoss.m_Mass));
			auto const itBNeuPeak = std::find_if(itBegin, itEnd, std::bind(&MassPeak_t::Approx, std::placeholders::_1, m_bTypeIon - g_rgAminoAcidsData[Methionine_Sulfoxide].m_NeutralLoss.m_Mass));
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
	vector<Cell_t> m_Solution{};
};

export template<uint16 iChargeFrom, uint16 iChargeTo>
constexpr decltype(auto) M2ZConversion(std::floating_point auto flFrom) noexcept
{
	using T = decltype(flFrom);

	constexpr int iChargeDiff = iChargeTo - iChargeFrom;

	flFrom *= T(iChargeFrom);
	flFrom += T(iChargeDiff) * amu::Hydrogen;
	flFrom /= T(iChargeTo);

	return flFrom;
}

AminoAcids_e TestNumber(double flPeakDiff) noexcept
{
	double flMinDiff = 9999.0;
	AminoAcids_e iFound = NOT_AN_AMINO_ACID;

	for (const auto& [AminoAcid, Data] : g_rgAminoAcidsData)
	{
		const double diff = std::abs(flPeakDiff - Data.m_ResidueMass);
		if (diff >= 0.5 || diff > flMinDiff)
			continue;

		iFound = (AminoAcids_e)AminoAcid;
		flMinDiff = diff;
	}

	return iFound;
}

export list<pair<string, double>> TestNumber2(double flPeakDiff) noexcept
{
	if (flPeakDiff <= DBL_EPSILON)
		return { { "SELF", 0.0 } };

	list<pair<string, double>> rgFound{};

	for (const auto& [AminoAcid, Data] : g_rgAminoAcidsData)
	{
		if (std::abs(flPeakDiff - Data.m_ResidueMass) < 1.0)
			rgFound.emplace_back(Data.m_3Letters + std::format("({})", AminoAcid), flPeakDiff - Data.m_ResidueMass);

		if (std::abs(flPeakDiff - Data.m_NeutralLoss.m_Mass) < 1.0)
			rgFound.emplace_back(Data.m_NeutralLoss.m_Name, flPeakDiff - Data.m_NeutralLoss.m_Mass);
	}

	rgFound.sort([](const pair<string, double>& lhs, const pair<string, double>& rhs) -> bool { return std::abs(lhs.second) < std::abs(rhs.second); });
	rgFound.unique([](const pair<string, double>& lhs, const pair<string, double>& rhs) -> bool { return lhs.first == rhs.first; });
	return rgFound;
}

export string TestNumber3(double flPeakDiff) noexcept
{
	if (flPeakDiff <= DBL_EPSILON)
		return "SELF";

	if (auto rgPossibilities = TestNumber2(flPeakDiff); !rgPossibilities.empty())
	{
		string sz;
		for (const auto& [szInfo, flDelta] : rgPossibilities)
			sz += std::format("{}[{:.4f}], ", szInfo, flDelta);

		sz.pop_back();
		sz.pop_back();
		return sz;
	}
	else
		return "-";
}

export double MolecularWeight(const string& szSeq) noexcept
{
	double ret = amu::Hydrogen * 2 + amu::Oxygen;	// Add H to -NH-, Add -OH to -CO-

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
		case Carboxyethyl_Cysteine:
		case Tyrosine:
		case Tryptophan:
			ret += g_rgAminoAcidsData[AminoAcid].m_ResidueMass;
			break;

		default:
			break;
		}
	}

	return ret;
}

// Get a polypeptide sequence as string.
export template<bool bFirmlySure = false, typename T = vector<Cell_t>>
string Conclude(const T& rgCells) noexcept
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

template<typename T, typename U> bool AllPeaksExplained(const T& rgflMassPeaks, const U& rgCells) noexcept
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

template<typename T, typename U> T FilterUnexplainablePeaks(const T& rgflMassPeaks, const U& rgCells) noexcept
{
	T ret;

	for (const auto& Peak : rgflMassPeaks)
	{
		if (std::find_if(rgCells.cbegin(), rgCells.cend(), std::bind(&Cell_t::Explained, std::placeholders::_1, Peak)) == rgCells.cend())
			ret.push_back(Peak);
	}

	return ret;
}

void fnRecursiveFromB1(list<AlternativeReality_t>& rgWorldlines, AlternativeReality_t& ThisWorldline, const function<bool(double)>& pfnShouldCheck) noexcept
{
	// The counterpart ion here is actually reaflecting the previous cell.
	bool bPredictedFound = false;
	double flPredictedCounterpartPeak = ThisWorldline.m_Solution.back().m_yTypeIon - g_rgAminoAcidsData[ThisWorldline.m_Solution.back().m_AminoAcid].m_ResidueMass;

	if (pfnShouldCheck(flPredictedCounterpartPeak))
	{
		for (const auto& Peak : ThisWorldline.m_PendingPeaks)
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
	double flFrontierPeak = ThisWorldline.m_Solution.back().m_bTypeIon;
	AlternativeReality_t ThisCopy = ThisWorldline;

	for (auto iter = ThisWorldline.m_PendingPeaks.rbegin(); iter != ThisWorldline.m_PendingPeaks.rend(); /* Does nothing. */)
	{
		if (*iter <= flFrontierPeak)
			goto LAB_CONTINUE;

		for (const auto& Cell : ThisWorldline.m_Solution)
		{
			if (Cell.Explained(*iter))
			{
				iter = decltype(iter)(ThisWorldline.m_PendingPeaks.erase(std::next(iter).base()));
				goto LAB_NEXT_CYCLE;
			}
		}

		if (AminoAcids_e iFound = TestNumber(*iter - flFrontierPeak); iFound != NOT_AN_AMINO_ACID)
		{
			if (!bAlreadyFoundOne)
			{
				bAlreadyFoundOne = true;
				ThisWorldline.m_Solution.emplace_back(iFound, iter->m_Value, flPredictedCounterpartPeak);

				iter = decltype(iter)(ThisWorldline.m_PendingPeaks.erase(std::next(iter).base()));
				continue;
			}
			else
			{
				AlternativeReality_t& OtherWorld = rgWorldlines.emplace_back(ThisCopy);
				OtherWorld.m_Solution.emplace_back(iFound, iter->m_Value, flPredictedCounterpartPeak);
				OtherWorld.m_PendingPeaks.erase(std::find(OtherWorld.m_PendingPeaks.begin(), OtherWorld.m_PendingPeaks.end(), *iter));

				fnRecursiveFromB1(rgWorldlines, OtherWorld, pfnShouldCheck);
			}
		}

	LAB_CONTINUE:;
		++iter;
	LAB_NEXT_CYCLE:;
	}

	if (bAlreadyFoundOne && !ThisWorldline.m_PendingPeaks.empty())
		fnRecursiveFromB1(rgWorldlines, ThisWorldline, pfnShouldCheck);
	else
		ThisWorldline.m_Solution.emplace_back();	// Place a sealer dummy.
}

void fnRecursiveFromBn(list<AlternativeReality_t>& rgWorldlines, AlternativeReality_t& ThisWorldline, const function<bool(double)>& pfnShouldCheck) noexcept
{
	bool bAlreadyFoundOne = false;
	double flFrontierPeak = ThisWorldline.m_Solution.front().m_bTypeIon;
	double flVarificationPeak = ThisWorldline.m_Solution[1].m_yTypeIon;
	AlternativeReality_t ThisCopy = ThisWorldline;

	for (auto iter = ThisWorldline.m_PendingPeaks.begin(); iter != ThisWorldline.m_PendingPeaks.end(); /* Does nothing. */)
	{
		if (*iter >= flFrontierPeak)
			goto LAB_CONTINUE;

		for (const auto& Cell : ThisWorldline.m_Solution)
		{
			if (Cell.Explained(*iter))
			{
				iter = ThisWorldline.m_PendingPeaks.erase(iter);
				goto LAB_NEXT_CYCLE;
			}
		}

		if (AminoAcids_e iFound = TestNumber(flFrontierPeak - *iter); iFound != NOT_AN_AMINO_ACID)
		{
			double flPredictedCounterpartPeak = flVarificationPeak + g_rgAminoAcidsData[iFound].m_ResidueMass;
			bool bPredictedFound = false;
			if (pfnShouldCheck(flPredictedCounterpartPeak))
			{
				for (const auto& Peak : ThisWorldline.m_PendingPeaks)
				{
					if (Peak.Approx(flPredictedCounterpartPeak))
					{
						bPredictedFound = true;
						flPredictedCounterpartPeak = Peak.m_Value;	// Set to the observed counterpart ion.
						break;
					}
				}

				if (!bPredictedFound)
					goto LAB_CONTINUE;
			}

			if (!bAlreadyFoundOne)
			{
				bAlreadyFoundOne = true;
				ThisWorldline.m_Solution.front().m_AminoAcid = iFound;
				ThisWorldline.m_Solution.front().m_yTypeIon = flPredictedCounterpartPeak;
				ThisWorldline.m_Solution.emplace(ThisWorldline.m_Solution.begin(), NOT_AN_AMINO_ACID, iter->m_Value);

				iter = ThisWorldline.m_PendingPeaks.erase(iter);
				continue;
			}
			else
			{
				AlternativeReality_t& OtherWorld = rgWorldlines.emplace_back(ThisCopy);
				OtherWorld.m_Solution.front().m_AminoAcid = iFound;
				OtherWorld.m_Solution.front().m_yTypeIon = flPredictedCounterpartPeak;
				OtherWorld.m_Solution.emplace(OtherWorld.m_Solution.begin(), NOT_AN_AMINO_ACID, iter->m_Value);
				OtherWorld.m_PendingPeaks.erase(std::find(OtherWorld.m_PendingPeaks.begin(), OtherWorld.m_PendingPeaks.end(), *iter));

				fnRecursiveFromBn(rgWorldlines, OtherWorld, pfnShouldCheck);
			}
		}

	LAB_CONTINUE:;
		++iter;
	LAB_NEXT_CYCLE:;
	}

	if (bAlreadyFoundOne && !ThisWorldline.m_PendingPeaks.empty())
		fnRecursiveFromBn(rgWorldlines, ThisWorldline, pfnShouldCheck);
}

void fnRecursiveFromY1(list<AlternativeReality_t>& rgWorldlines, AlternativeReality_t& ThisWorldline, const function<bool(double)>& pfnShouldCheck) noexcept
{
	// The counterpart ion here is actually reaflecting the previous cell.
	bool bPredictedFound = false;
	double flPredictedCounterpartPeak = ThisWorldline.m_Solution[0].m_bTypeIon - g_rgAminoAcidsData[ThisWorldline.m_Solution[0].m_AminoAcid].m_ResidueMass;
	
	if (pfnShouldCheck(flPredictedCounterpartPeak))
	{
		for (const auto& Peak : ThisWorldline.m_PendingPeaks)
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
	double flFrontierPeak = ThisWorldline.m_Solution[0].m_yTypeIon;
	AlternativeReality_t ThisCopy = ThisWorldline;

	for (auto iter = ThisWorldline.m_PendingPeaks.rbegin(); iter != ThisWorldline.m_PendingPeaks.rend(); /* Does nothing. */)
	{
		if (*iter <= flFrontierPeak)
			goto LAB_CONTINUE;

		for (const auto& Cell : ThisWorldline.m_Solution)
		{
			if (Cell.Explained(*iter))
			{
				iter = decltype(iter)(ThisWorldline.m_PendingPeaks.erase(std::next(iter).base()));
				goto LAB_NEXT_CYCLE;
			}
		}

		if (AminoAcids_e iFound = TestNumber(*iter - flFrontierPeak); iFound != NOT_AN_AMINO_ACID)
		{
			if (!bAlreadyFoundOne)
			{
				bAlreadyFoundOne = true;
				ThisWorldline.m_Solution.emplace(ThisWorldline.m_Solution.begin(), iFound, flPredictedCounterpartPeak, iter->m_Value);

				iter = decltype(iter)(ThisWorldline.m_PendingPeaks.erase(std::next(iter).base()));
				continue;
			}
			else
			{
				AlternativeReality_t& OtherWorld = rgWorldlines.emplace_back(ThisCopy);
				OtherWorld.m_Solution.emplace(OtherWorld.m_Solution.begin(), iFound, flPredictedCounterpartPeak, iter->m_Value);
				OtherWorld.m_PendingPeaks.erase(std::find(OtherWorld.m_PendingPeaks.begin(), OtherWorld.m_PendingPeaks.end(), *iter));

				fnRecursiveFromY1(rgWorldlines, OtherWorld, pfnShouldCheck);
			}
		}

	LAB_CONTINUE:;
		++iter;
	LAB_NEXT_CYCLE:;
	}

	if (bAlreadyFoundOne && !ThisWorldline.m_PendingPeaks.empty())
		fnRecursiveFromY1(rgWorldlines, ThisWorldline, pfnShouldCheck);
	else
		ThisWorldline.m_Solution.emplace(ThisWorldline.m_Solution.begin());	// Place a sealer dummy.
}

void fnRecursiveFromYn(list<AlternativeReality_t>& rgWorldlines, AlternativeReality_t& ThisWorldline, const function<bool(double)>& pfnShouldCheck) noexcept
{
	bool bAlreadyFoundOne = false;
	double flFrontierPeak = ThisWorldline.m_Solution.back().m_yTypeIon;
	double flVarificationPeak = ThisWorldline.m_Solution[ThisWorldline.m_Solution.size() - 2U].m_bTypeIon;
	AlternativeReality_t ThisCopy = ThisWorldline;

	for (auto iter = ThisWorldline.m_PendingPeaks.begin(); iter != ThisWorldline.m_PendingPeaks.end(); /* Does nothing. */)
	{
		if (*iter >= flFrontierPeak)
			goto LAB_CONTINUE;

		for (const auto& Cell : ThisWorldline.m_Solution)
		{
			if (Cell.Explained(*iter))
			{
				iter = ThisWorldline.m_PendingPeaks.erase(iter);
				goto LAB_NEXT_CYCLE;
			}
		}

		if (AminoAcids_e iFound = TestNumber(flFrontierPeak - *iter); iFound != NOT_AN_AMINO_ACID)
		{
			double flPredictedCounterpartPeak = flVarificationPeak + g_rgAminoAcidsData[iFound].m_ResidueMass;
			bool bPredictedFound = false;
			if (pfnShouldCheck(flPredictedCounterpartPeak))
			{
				for (const auto& Peak : ThisWorldline.m_PendingPeaks)
				{
					if (Peak.Approx(flPredictedCounterpartPeak))
					{
						bPredictedFound = true;
						flPredictedCounterpartPeak = Peak.m_Value;	// Set to the observed counterpart ion.
						break;
					}
				}

				if (!bPredictedFound)
					goto LAB_CONTINUE;
			}

			if (!bAlreadyFoundOne)
			{
				bAlreadyFoundOne = true;
				ThisWorldline.m_Solution.back().m_AminoAcid = iFound;
				ThisWorldline.m_Solution.back().m_bTypeIon = flPredictedCounterpartPeak;
				ThisWorldline.m_Solution.emplace_back(NOT_AN_AMINO_ACID, 0, iter->m_Value);

				iter = ThisWorldline.m_PendingPeaks.erase(iter);
				continue;
			}
			else
			{
				AlternativeReality_t& OtherWorld = rgWorldlines.emplace_back(ThisCopy);
				OtherWorld.m_Solution.back().m_AminoAcid = iFound;
				OtherWorld.m_Solution.back().m_bTypeIon = flPredictedCounterpartPeak;
				OtherWorld.m_Solution.emplace_back(NOT_AN_AMINO_ACID, 0, iter->m_Value);
				OtherWorld.m_PendingPeaks.erase(std::find(OtherWorld.m_PendingPeaks.begin(), OtherWorld.m_PendingPeaks.end(), *iter));

				fnRecursiveFromYn(rgWorldlines, OtherWorld, pfnShouldCheck);
			}
		}

	LAB_CONTINUE:;
		++iter;
	LAB_NEXT_CYCLE:;
	}

	if (bAlreadyFoundOne && !ThisWorldline.m_PendingPeaks.empty())
		fnRecursiveFromYn(rgWorldlines, ThisWorldline, pfnShouldCheck);
}

template<typename T> void PrintWorldline(const AlternativeReality_t& Worldline, const T& rgflOriginalMassData) noexcept
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
		bool bbFound = std::find_if(itBegin, itEnd, std::bind(&MassPeak_t::Approx, std::placeholders::_1, Cell.m_bTypeIon)) != itEnd;
		bool byFound = std::find_if(itBegin, itEnd, std::bind(&MassPeak_t::Approx, std::placeholders::_1, Cell.m_yTypeIon)) != itEnd;

		if (bbFound && byFound)
			cout_g() << Cell.m_AminoAcid << '\t';
		else if (bbFound || byFound)
			cout_gold() << Cell.m_AminoAcid << '\t';
		else
			cout_r() << Cell.m_AminoAcid << '\t';
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
			std::cout << std::format("{} is the b[n-1] ion and the C-terminal amino acid is {}.\n", iter->m_Value, g_rgAminoAcidsData[iNeutralLoss == 146 ? Lysine : Arginine].m_3Letters);
			FirstCWorldline.m_Solution.emplace(FirstCWorldline.m_Solution.begin(), iNeutralLoss == 146 ? Lysine : Arginine, M_Plus_1 - 18 /* b[n] is [M+1] ion with a H2O loss. */, iNeutralLoss + 1);
			FirstCWorldline.m_Solution.emplace(FirstCWorldline.m_Solution.begin(), NOT_AN_AMINO_ACID, iter->m_Value /* b[n-1] */);
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
				std::cout << std::format("{} is the y1 ion and the C-terminal amino acid is {}.\n", iter->m_Value, g_rgAminoAcidsData[iMass == 147 ? Lysine : Arginine].m_3Letters);
				FirstCWorldline.m_Solution.emplace_back(iMass == 147 ? Lysine : Arginine, M_Plus_1 - 18 /* b[n] is [M+1] ion with a H2O loss. */, *iter);
				iter = FirstCWorldline.m_PendingPeaks.erase(iter);
			}
			else
				++iter;
		}

		fnRecursiveFromY1(rgCTerminalGuesses, FirstCWorldline, fnCheck);
	}
	else
		fnRecursiveFromBn(rgCTerminalGuesses, FirstCWorldline, fnCheck);
#pragma endregion C-Terminal attempts

#pragma region N-Terminal attempts
	list<AlternativeReality_t> rgNTerminalGuesses {};
	for (auto iter = rgflMassData2.cbegin(); iter != rgflMassData2.cend(); ++iter)
	{
		if (AminoAcids_e iFound = TestNumber(M_Plus_1 - *iter); iFound != NOT_AN_AMINO_ACID)
		{
			std::cout << std::format("{} could be the y[n-1] ion and N-terminal amino acid could be {}.\n", iter->m_Value, g_rgAminoAcidsData[iFound].m_3Letters);

			AlternativeReality_t& AnotherNWorldline = rgNTerminalGuesses.emplace_back();
			AnotherNWorldline.m_PendingPeaks = rgflMassData2;
			AnotherNWorldline.m_MPlusOne = M_Plus_1;

			AnotherNWorldline.m_Solution.emplace_back(iFound, g_rgAminoAcidsData[iFound].m_ResidueMass + amu::Hydrogen, M_Plus_1 /* y[n] is [M+1] ion. */);
			AnotherNWorldline.m_Solution.emplace_back(NOT_AN_AMINO_ACID, 0, iter->m_Value /* y[n-1] */);
			AnotherNWorldline.m_PendingPeaks.erase(std::find(AnotherNWorldline.m_PendingPeaks.begin(), AnotherNWorldline.m_PendingPeaks.end(), *iter));

			fnRecursiveFromYn(rgNTerminalGuesses, AnotherNWorldline, fnCheck);
		}
	}

	if (rgNTerminalGuesses.empty())
	{
		for (auto iter = rgflMassData2.cbegin(); iter != rgflMassData2.cend(); ++iter)
		{
			if (AminoAcids_e iFound = TestNumber(*iter - amu::Hydrogen); iFound != NOT_AN_AMINO_ACID)
			{
				std::cout << std::format("{} could be the b1 ion and N-terminal amino acid could be {}.\n", iter->m_Value, g_rgAminoAcidsData[iFound].m_3Letters);

				AlternativeReality_t& AnotherNWorldline = rgNTerminalGuesses.emplace_back();
				AnotherNWorldline.m_PendingPeaks = rgflMassData2;
				AnotherNWorldline.m_MPlusOne = M_Plus_1;

				AnotherNWorldline.m_Solution.emplace_back(iFound, iter->m_Value, M_Plus_1 /* y[n] is [M+1] ion. */);
				AnotherNWorldline.m_PendingPeaks.erase(std::find(AnotherNWorldline.m_PendingPeaks.begin(), AnotherNWorldline.m_PendingPeaks.end(), *iter));

				fnRecursiveFromB1(rgNTerminalGuesses, AnotherNWorldline, fnCheck);
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
			for (const char* pszSeq = szSequence.c_str(); *pszSeq != '\0'; ++pszSeq)
			{
				if (auto iLength = strlen(pszSeq); iLength > 1 && !strncmp(pszSeq, pszSeqCompW, iLength))
				{
					string s = szSequence;
					s.erase(s.length() - iLength);
					s += szSeqCompareWith;

					if (std::round(MolecularWeight(s) + amu::Hydrogen - M_Plus_1))
						continue;

					auto& MergedWorldline = rgExplanations.emplace_back();
					MergedWorldline.m_MPlusOne = M_Plus_1;

					auto& Solution = MergedWorldline.m_Solution;
					Solution.insert(Solution.end(), Worldline.m_Solution.begin(), Worldline.m_Solution.begin() + (szSequence.length() - iLength));
					Solution.insert(Solution.end(), WorldlineCompareWith.m_Solution.begin() + 1, WorldlineCompareWith.m_Solution.end());

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

export template<typename T>
void MarkPeaks(const AlternativeReality_t& Worldline, T& rgflMassPeaks)
{
	short iCount = 1;
	short const iTotalCount = (short)Worldline.m_Solution.size() + 1;
	for (const auto& Cell : Worldline.m_Solution)
	{
		for (auto& Peak : rgflMassPeaks)
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
	for (auto& Peak : rgflMassPeaks)
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

export template<typename T>
void ResetPeaks(T& rgflMassPeaks)
{
	for (auto& Peak : rgflMassPeaks)
		Peak.m_String.clear();
}
