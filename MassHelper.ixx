module;

#include <unordered_map>
#include <array>
#include <vector>
#include <format>
#include <cfloat>
#include <concepts>

#include <iostream>

#include <gcem.hpp>

export module MassHelper;

using std::unordered_map;
using std::array;
using std::pair;
using std::vector;
using std::string;
using std::list;

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
	constexpr MassPeak_t(std::floating_point auto v) noexcept : m_Value(v) {}

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
	constexpr auto operator<=> (std::floating_point auto rhs) noexcept { return m_Value <=> rhs; }
	constexpr bool operator== (const MassPeak_t& rhs) noexcept { return gcem::abs(m_Value - rhs.m_Value) <= DBL_EPSILON; }
	constexpr operator double& () noexcept { return m_Value; }
	constexpr operator const double& () const noexcept { return m_Value; }
};

AminoAcids_e TestNumber(double flPeakDiff) noexcept
{
	double flMinDiff = 9999.0;
	AminoAcids_e iFound = NOT_AN_AMINO_ACID;

	for (const auto& [AminoAcid, Data] : g_rgAminoAcidsData)
	{
		const double diff = std::abs(flPeakDiff - Data.m_ResidueMass);
		if (diff >= 1.0 || diff > flMinDiff)
			continue;

		iFound = (AminoAcids_e)AminoAcid;
		flMinDiff = diff;
	}

	return iFound;
}

export void IdentifyBorderIons(vector<MassPeak_t>& rgflMassData, double M_plus_1) noexcept
{
	int iMPlusOne = (int)std::round(M_plus_1);
	int iMPlusTwo = (int)std::round((M_plus_1 + 1.00784) / 2.0);

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

export void RecursiveIdentify(vector<MassPeak_t>& rgflMassData, IonType iType, short iIndex) noexcept
{
	if (!iIndex)
	{
		std::cout << "iIndex parameter reaches 0.\n";
		return;
	}

	short iLastIdentified = iIndex < 0 ? (iIndex + 1) : (iIndex - 1);
	const MassPeak_t* pLastIdentified = nullptr;

	for (const auto& Peak : rgflMassData)
	{
		if (Peak.m_Count == iLastIdentified && Peak.m_Type == iType)
		{
			pLastIdentified = &Peak;
			break;
		}
	}

	if (!pLastIdentified)
		goto LAB_EXIT_IDENTIFY;

	for (auto& Peak : rgflMassData)
	{
		if (Peak.m_Identified)
			continue;

		if ((iIndex < 0 && Peak >= *pLastIdentified) || (iIndex > 0 && Peak <= *pLastIdentified))
			continue;

		if (AminoAcids_e iFound = TestNumber(std::abs(Peak - *pLastIdentified)); iFound != NOT_AN_AMINO_ACID)
		{
			Peak.Identify(iType, iIndex);
			std::cout << std::format("Peak {0} was identified as {1} ion. Loss {2}[{3}] from {4} to {1}.\n", Peak.m_Value, Peak.ToString(), g_rgAminoAcidsData[iFound].m_3Letters, (char)iFound, pLastIdentified->ToString());
			return RecursiveIdentify(rgflMassData, iType, iIndex < 0 ? (iIndex - 1) : (iIndex + 1));
		}
	}

LAB_EXIT_IDENTIFY:;
	std::cout << std::format("No {}{}{}{} ion found.\n", (char)iType, iIndex < 0 ? "[n" : "", iIndex, iIndex < 0 ? "]" : "");
}

export
template<IonType iType>
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
