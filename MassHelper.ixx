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

using std::array;
using std::function;
using std::list;
using std::pair;
using std::string;
using std::unordered_map;
using std::vector;

export constexpr double HYDROGEN_AMU = 1.00784;

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

constexpr Molecule_t NO_LOSS = {};
constexpr Molecule_t WATER = { 18, "H2O" };
constexpr Molecule_t AMMONIA = { 17, "NH3" };

unordered_map<char, AminoAcid_t> g_rgAminoAcidsData =
{
	{AminoAcids_e::NOT_AN_AMINO_ACID,		{}},
	{AminoAcids_e::Glycine,					{57.0215, 30, "Gly", NO_LOSS}},
	{AminoAcids_e::Alanine,					{71.0371, 44, "Ala", NO_LOSS}},
	{AminoAcids_e::Serine,					{87.0320, 60, "Ser", WATER}},
	{AminoAcids_e::Proline,					{97.0528, 70, "Pro", NO_LOSS}},
	{AminoAcids_e::Valine,					{99.0684, 72, "Val", NO_LOSS}},
	{AminoAcids_e::Threonine,				{101.0477, 74, "Thr", WATER}},
	{AminoAcids_e::Cysteine,				{103.0092, 76, "Cys", { 34, "SH2" }}},
	{AminoAcids_e::Isoleucine,				{113.0841, 86, "Ile", NO_LOSS}},
	{AminoAcids_e::Leucine,					{113.0841, 86, "Leu", NO_LOSS}},
	{AminoAcids_e::Asparagine,				{114.0429, 87, "Asn", AMMONIA}},
	{AminoAcids_e::Aspartic_Acid,			{115.0269, 88, "Asp", WATER}},
	{AminoAcids_e::Lysine,					{128.0950, 101, "Lys", AMMONIA}},
	{AminoAcids_e::Glutamine,				{128.0586, 101, "Gln", AMMONIA}},
	{AminoAcids_e::Glutamic_Acid,			{129.0426, 102, "Glu", WATER}},
	{AminoAcids_e::Methionine,				{131.0405, 104, "Met", { 48, "CH3SH" }}},
	{AminoAcids_e::Histidine,				{137.0589, 110, "His", NO_LOSS}},
	{AminoAcids_e::Methionine_Sulfoxide,	{147.0354, 120, "Mox", { 64, "CH3SOH" }}},
	{AminoAcids_e::Phenylalanine,			{147.0684, 120, "Phe", NO_LOSS}},
	{AminoAcids_e::Arginine,				{156.0684, 129, "Arg", AMMONIA}},
	{AminoAcids_e::Carboxyethyl_Cysteine,	{161.0147, 134, "Cec", { 92, "HSCH2COOH" }}},
	{AminoAcids_e::Tyrosine,				{163.0633, 134, "Tyr", NO_LOSS}},
	{AminoAcids_e::Tryptophan,				{186.0793, 159, "Trp", NO_LOSS}},
};

export enum IonType : char
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

export struct MassPeak_t
{
	constexpr MassPeak_t() noexcept {}
	constexpr MassPeak_t(Arithmetic auto v) noexcept : m_Value(v) {}

	double m_Value = 0;
	bool m_Identified = false;
	IonType m_Type = IonType::UNKNOWN_TYPE;
	short m_Count = 0;

	string ToString(void) const noexcept
	{
		if (!m_Identified)
			return "Unknown";

		switch (m_Type)
		{
		case IonType::M_PLUS_2H:
			return "[M+2H]";

		case IonType::M_PLUS_H:
			return "[M+H]";

		default:
			return (m_Count > 0) ? std::format("{}{}", (char)m_Type, m_Count) : std::format("{}[n{}]", (char)m_Type, m_Count);
		}
	}

	void Identify(IonType iType, short iCount = 0) noexcept
	{
		if (m_Identified)
			std::cout << std::format("Ion {} was previously identified as {}.\n", m_Value, ToString());

		m_Type = iType;
		m_Count = iCount;
		m_Identified = true;
	}

	constexpr auto operator<=> (const MassPeak_t& rhs) noexcept { return m_Value <=> rhs.m_Value; }
	constexpr auto operator<=> (Arithmetic auto rhs) noexcept { return m_Value <=> rhs; }
	constexpr bool operator== (const MassPeak_t& rhs) noexcept { return gcem::abs(m_Value - rhs.m_Value) <= DBL_EPSILON; }
	constexpr bool operator== (Arithmetic auto rhs) noexcept { return gcem::abs(m_Value - rhs) <= DBL_EPSILON; }
	constexpr operator double& () noexcept { return m_Value; }
	constexpr operator const double& () const noexcept { return m_Value; }
};

struct Cell_t
{
	constexpr Cell_t() noexcept {}
	constexpr Cell_t(AminoAcids_e what, double b = 0, double y = 0) noexcept : m_AminoAcid(what), m_bTypeIon(b), m_yTypeIon(y) {}

	double m_aTypeIon = 0;
	double m_bAsteriskIon = 0;
	double m_bCircleIon = 0;
	double m_bTypeIon = 0;	// Including this amino acid Coule be b[n] i.e. [M+1] ion itself.
	AminoAcids_e m_AminoAcid = NOT_AN_AMINO_ACID;
	double m_yTypeIon = 0;	// Including this amino acid. Could be y[n] i.e. [M+1] ion itself.
	double m_yCircleIon = 0;
	double m_yAsteriskIon = 0;

	constexpr bool Filled(void) const noexcept { return m_AminoAcid != NOT_AN_AMINO_ACID && (m_bTypeIon != 0 || m_yTypeIon != 0); }

	constexpr bool operator==(const Cell_t& rhs) const noexcept { return m_AminoAcid == rhs.m_AminoAcid; }
};

export struct AlternativeReality_t
{
	double m_MPlusOne = 0;
	vector<MassPeak_t> m_PendingPeaks{};
	list<Cell_t> m_Solution{};
};

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

export void IdentifyBorderIons(vector<MassPeak_t>& rgflMassData, double M_plus_1) noexcept
{
	int iMPlusOne = (int)std::round(M_plus_1);
	int iMPlusTwo = (int)std::round((M_plus_1 + HYDROGEN_AMU) / 2.0);

	for (auto& Peak : rgflMassData)
	{
		// Test 1: y1 ion?
		switch ((int)std::round(Peak))
		{
		case 147:	// residue mass + H3O[+]
			std::cout << std::format("{} is the y1 ion and the C-terminal amino acid is {}.\n", Peak.m_Value, g_rgAminoAcidsData[Lysine].m_3Letters);
			Peak.Identify(IonType::y, 1);
			continue;

		case 175:
			std::cout << std::format("{} is the y1 ion and the C-terminal amino acid is {}.\n", Peak.m_Value, g_rgAminoAcidsData[Arginine].m_3Letters);
			Peak.Identify(IonType::y, 1);
			continue;

		default:
			break;
		}

		for (const auto& [AminoAcid, Data] : g_rgAminoAcidsData)
		{
			// Test 2: a1 ion?
			if ((int)std::round(Peak) == Data.m_ImmoniumIonMass)
			{
				std::cout << std::format("Peak {} was identified as a1 ion and the N-terminal amino acid is {}.\n", Peak.m_Value, Data.m_3Letters);
				Peak.Identify(IonType::a, 1);
				break;
			}

			// Test 3: b1 ion?
			// Note: Chemically speaking, you won't see b1 ion for some reason.
			else if ((int)std::round(Peak) == Data.m_ImmoniumIonMass + 28)
			{
				std::cout << std::format("Peak {} was identified as b1 ion and the N-terminal amino acid is {}.\n", Peak.m_Value, Data.m_3Letters);
				Peak.Identify(IonType::b, 1);
				break;
			}
		}

		// Test 4: y[n-1] ion?
		if (auto AminoAcid = TestNumber(M_plus_1 - Peak); AminoAcid != NOT_AN_AMINO_ACID)
		{
			std::cout << std::format("Peak {} was identified as y[n-1] ion and the N-terminal amino acid is {}.\n", Peak.m_Value, g_rgAminoAcidsData[AminoAcid].m_3Letters);
			Peak.Identify(IonType::y, -1);
			continue;
		}

		// Test 5: b[n-1] ion?
		switch ((int)std::round(M_plus_1 - Peak))
		{
		case 146:	// residue mass + H2O (Consider as natural loss).
			std::cout << std::format("{} is the b[n-1] ion and the C-terminal amino acid is {}.\n", Peak.m_Value, g_rgAminoAcidsData[Lysine].m_3Letters);
			Peak.Identify(IonType::b, -1);
			continue;

		case 174:
			std::cout << std::format("{} is the b[n-1] ion and the C-terminal amino acid is {}.\n", Peak.m_Value, g_rgAminoAcidsData[Arginine].m_3Letters);
			Peak.Identify(IonType::b, -1);
			continue;

		default:
			break;
		}

		// Test 6: [M+H] or [M+2H]?
		if (int iPeakMass = (int)std::round(Peak); iPeakMass == iMPlusOne || iPeakMass == iMPlusTwo)
		{
			Peak.Identify(iPeakMass == iMPlusOne ? IonType::M_PLUS_H : IonType::M_PLUS_2H);
			continue;
		}
	}
}

void fnRecursiveB(list<AlternativeReality_t>& rgWorldlines, AlternativeReality_t& ThisWorldline, const function<bool(double)>& pfnShouldCheck) noexcept
{
	bool bAlreadyFoundOne = false;
	double flLastBIon = ThisWorldline.m_Solution.front().m_bTypeIon;
	AlternativeReality_t ThisCopy = ThisWorldline;

	for (auto iter = ThisWorldline.m_PendingPeaks.begin(); iter != ThisWorldline.m_PendingPeaks.end(); /* Does nothing. */)
	{
		if (*iter >= flLastBIon)
			goto LAB_CONTINUE;

		if (AminoAcids_e iFound = TestNumber(flLastBIon - *iter); iFound != NOT_AN_AMINO_ACID)
		{
			double flAccumulatedMass = HYDROGEN_AMU + g_rgAminoAcidsData[iFound].m_ResidueMass;
			for (const auto& Cell : ThisWorldline.m_Solution)
				flAccumulatedMass += g_rgAminoAcidsData[Cell.m_AminoAcid].m_ResidueMass;

			bool bPredictedFound = false;
			//if (pfnShouldCheck(flAccumulatedMass))
			//{
				for (const auto& Peak : ThisWorldline.m_PendingPeaks)
				{
					if ((int)std::round(flAccumulatedMass) == (int)std::round(Peak))
					{
						bPredictedFound = true;
						flAccumulatedMass = Peak.m_Value;	// Set to counterpart Y ion.
						break;
					}
				}

			//	if (!bPredictedFound)
			//		goto LAB_CONTINUE;
			//}

			if (!bAlreadyFoundOne)
			{
				bAlreadyFoundOne = true;
				ThisWorldline.m_Solution.front().m_AminoAcid = iFound;
				ThisWorldline.m_Solution.emplace_front(NOT_AN_AMINO_ACID, iter->m_Value, flAccumulatedMass);

				iter = ThisWorldline.m_PendingPeaks.erase(iter);
				continue;
			}
			else
			{
				rgWorldlines.emplace_back(ThisCopy);

				AlternativeReality_t& OtherWorld = rgWorldlines.back();
				OtherWorld.m_Solution.front().m_AminoAcid = iFound;
				OtherWorld.m_Solution.emplace_front(NOT_AN_AMINO_ACID, iter->m_Value, flAccumulatedMass);
				OtherWorld.m_PendingPeaks.erase(std::find(OtherWorld.m_PendingPeaks.begin(), OtherWorld.m_PendingPeaks.end(), *iter));

				fnRecursiveB(rgWorldlines, OtherWorld, pfnShouldCheck);
			}
		}

	LAB_CONTINUE:;
		++iter;
	}

	if (bAlreadyFoundOne && !ThisWorldline.m_PendingPeaks.empty())
		fnRecursiveB(rgWorldlines, ThisWorldline, pfnShouldCheck);
}

void fnRecursiveY(list<AlternativeReality_t>& rgWorldlines, AlternativeReality_t& ThisWorldline, const function<bool(double)>& pfnShouldCheck) noexcept
{
	bool bAlreadyFoundOne = false;
	double flLastYIon = ThisWorldline.m_Solution.back().m_yTypeIon;
	AlternativeReality_t ThisCopy = ThisWorldline;

	for (auto iter = ThisWorldline.m_PendingPeaks.begin(); iter != ThisWorldline.m_PendingPeaks.end(); /* Does nothing. */)
	{
		if (*iter >= flLastYIon)
			goto LAB_CONTINUE;

		if (AminoAcids_e iFound = TestNumber(flLastYIon - *iter); iFound != NOT_AN_AMINO_ACID)
		{
			double flAccumulatedMass = HYDROGEN_AMU + g_rgAminoAcidsData[iFound].m_ResidueMass;
			for (const auto& Cell : ThisWorldline.m_Solution)
				flAccumulatedMass += g_rgAminoAcidsData[Cell.m_AminoAcid].m_ResidueMass;

			bool bPredictedFound = false;
			//if (pfnShouldCheck(flAccumulatedMass))
			//{
				for (const auto& Peak : ThisWorldline.m_PendingPeaks)
				{
					if ((int)std::round(flAccumulatedMass) == (int)std::round(Peak))
					{
						bPredictedFound = true;
						flAccumulatedMass = Peak.m_Value;	// Set to counterpart B ion.
						break;
					}
				}

			//	if (!bPredictedFound)
			//		goto LAB_CONTINUE;
			//}

			if (!bAlreadyFoundOne)
			{
				bAlreadyFoundOne = true;
				ThisWorldline.m_Solution.back().m_AminoAcid = iFound;
				ThisWorldline.m_Solution.emplace_back(NOT_AN_AMINO_ACID, flAccumulatedMass, iter->m_Value);

				iter = ThisWorldline.m_PendingPeaks.erase(iter);
				continue;
			}
			else
			{
				rgWorldlines.emplace_back(ThisCopy);

				AlternativeReality_t& OtherWorld = rgWorldlines.back();
				OtherWorld.m_Solution.back().m_AminoAcid = iFound;
				OtherWorld.m_Solution.emplace_back(NOT_AN_AMINO_ACID, flAccumulatedMass, iter->m_Value);
				OtherWorld.m_PendingPeaks.erase(std::find(OtherWorld.m_PendingPeaks.begin(), OtherWorld.m_PendingPeaks.end(), *iter));

				fnRecursiveY(rgWorldlines, OtherWorld, pfnShouldCheck);
			}
		}

	LAB_CONTINUE:;
		++iter;
	}

	if (bAlreadyFoundOne && !ThisWorldline.m_PendingPeaks.empty())
		fnRecursiveY(rgWorldlines, ThisWorldline, pfnShouldCheck);
}

export void Solve(const vector<MassPeak_t>& rgflMassData, double M_Plus_1) noexcept
{
	vector<MassPeak_t> rgflMassData2 = rgflMassData;
	int iMPlusOne = (int)std::round(M_Plus_1);
	int iMPlusTwo = (int)std::round((M_Plus_1 + HYDROGEN_AMU) / 2.0);

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

	list<AlternativeReality_t> rgWorldlines{ AlternativeReality_t{} };
	AlternativeReality_t& FirstBWorldline = rgWorldlines.front();
	FirstBWorldline.m_PendingPeaks = rgflMassData2;
	FirstBWorldline.m_MPlusOne = M_Plus_1;

	for (auto iter = FirstBWorldline.m_PendingPeaks.begin(); iter != FirstBWorldline.m_PendingPeaks.end(); /* Does nothing. */)
	{
		if (int iNeutralLoss = (int)std::round(M_Plus_1 - *iter); iNeutralLoss == 146 || iNeutralLoss == 174)
		{
			std::cout << std::format("{} is the b[n-1] ion and the C-terminal amino acid is {}.\n", iter->m_Value, g_rgAminoAcidsData[iNeutralLoss == 146 ? Lysine : Arginine].m_3Letters);
			FirstBWorldline.m_Solution.emplace_front(iNeutralLoss == 146 ? Lysine : Arginine, M_Plus_1 - 18 /* b[n] is [M+1] ion with a H2O loss. */, *iter);
			FirstBWorldline.m_Solution.emplace_front(NOT_AN_AMINO_ACID, iter->m_Value /* b[n-1] */);
			iter = FirstBWorldline.m_PendingPeaks.erase(iter);
			continue;
		}

		++iter;
	}

	if (FirstBWorldline.m_Solution.empty())
	{
		for (auto iter = FirstBWorldline.m_PendingPeaks.begin(); iter != FirstBWorldline.m_PendingPeaks.end(); /* Does nothing. */)
		{
			if (int iMass = (int)std::round(*iter); iMass == 147 || iMass == 175)
			{
				std::cout << std::format("{} is the y1 ion and the C-terminal amino acid is {}.\n", iter->m_Value, g_rgAminoAcidsData[iMass == 147 ? Lysine : Arginine].m_3Letters);
				FirstBWorldline.m_Solution.emplace_front(iMass == 147 ? Lysine : Arginine, M_Plus_1 - 18 /* b[n] is [M+1] ion with a H2O loss. */, *iter);
				iter = FirstBWorldline.m_PendingPeaks.erase(iter);
				continue;
			}

			++iter;
		}
	}

	assert(FirstBWorldline.m_Solution.size() == 1 || FirstBWorldline.m_Solution.size() == 2);

	fnRecursiveB(rgWorldlines, FirstBWorldline, [rgflMassData2](double flAccumulatedMass) -> bool { return flAccumulatedMass >= rgflMassData2.back(); });

	for (auto iter = rgflMassData2.begin(); iter != rgflMassData2.end(); ++iter)
	{
		if (AminoAcids_e iFound = TestNumber(M_Plus_1 - *iter); iFound != NOT_AN_AMINO_ACID)
		{
			std::cout << std::format("{} could be the y[n-1] ion and N-terminal amino acid could be {}.\n", iter->m_Value, g_rgAminoAcidsData[iFound].m_3Letters);

			rgWorldlines.emplace_back();
			AlternativeReality_t& AnotherYWorldline = rgWorldlines.back();
			AnotherYWorldline.m_PendingPeaks = rgflMassData2;
			AnotherYWorldline.m_MPlusOne = M_Plus_1;

			AnotherYWorldline.m_Solution.emplace_back(iFound, 0, M_Plus_1 /* y[n] is [M+1] ion. */);
			AnotherYWorldline.m_Solution.emplace_back(NOT_AN_AMINO_ACID, 0, iter->m_Value /* y[n-1] */);
			AnotherYWorldline.m_PendingPeaks.erase(std::find(AnotherYWorldline.m_PendingPeaks.begin(), AnotherYWorldline.m_PendingPeaks.end(), *iter));

			fnRecursiveY(rgWorldlines, AnotherYWorldline, [rgflMassData2](double flAccumulatedMass) -> bool { return flAccumulatedMass >= rgflMassData2.back(); });
		}
	}

	for (const auto& Worldline : rgWorldlines)
	{
		for (const auto& Cell : Worldline.m_Solution)
		{
			if (std::find(rgflMassData2.cbegin(), rgflMassData2.cend(), Cell.m_bTypeIon) != rgflMassData2.cend())
				cout_w() << Cell.m_bTypeIon << '\t';
			else
				cout_pink() << Cell.m_bTypeIon << '\t';
		}

		cout_w() << '\n';

		for (const auto& Cell : Worldline.m_Solution)
		{
			if (std::find(rgflMassData2.cbegin(), rgflMassData2.cend(), Cell.m_bTypeIon) != rgflMassData2.cend()
				&& std::find(rgflMassData2.cbegin(), rgflMassData2.cend(), Cell.m_yTypeIon) != rgflMassData2.cend())
				cout_g() << Cell.m_AminoAcid << '\t';
			else
				cout_gold() << Cell.m_AminoAcid << '\t';
		}

		cout_w() << '\n';
		
		for (const auto& Cell : Worldline.m_Solution)
		{
			if (std::find(rgflMassData2.cbegin(), rgflMassData2.cend(), Cell.m_yTypeIon) != rgflMassData2.cend())
				cout_w() << Cell.m_yTypeIon << '\t';
			else
				cout_pink() << Cell.m_yTypeIon << '\t';
		}

		cout_gray() << "\nLeft peaks: ";
		for (const auto& Peak : Worldline.m_PendingPeaks)
		{
			cout_gray() << Peak.m_Value << ", ";
		}

		cout_w() << '\n';
		cout_w() << '\n';
	}
}

export template<IonType iType>
void ParseSpectrum(vector<MassPeak_t>& rgflMassData, double M_Plus_1) noexcept
{
	list<AminoAcids_e> rg;
	double flLast = M_Plus_1;

	for (auto iter = rgflMassData.crbegin(); iter != rgflMassData.crend(); ++iter)
	{
		if (!iter->m_Identified)
			continue;

		if (iter->m_Type != iType)
			continue;

		if constexpr (iType == IonType::x || iType == IonType::y || iType == IonType::z)
			rg.push_back(TestNumber(std::abs(flLast - *iter)));
		else if constexpr (iType == IonType::a || iType == IonType::b || iType == IonType::c)
			rg.push_front(TestNumber(std::abs(flLast - *iter)));

		flLast = *iter;
	}

	for (const auto& AminoAcid : rg)
	{
		std::cout << AminoAcid;
	}

	std::cout << '\n';
}
