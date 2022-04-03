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

export enum IonType : unsigned short
{
	UNKNOWN_TYPE = 0,
	
	// N-terminal
	a, // C-C, alpha C and carboxylic C, counterpart x
	b, // C-N, peptide bond, counterpart y
	c, // N-C, amine group from next amino acid and its alpha C, counterpart z

	// C-terminal
	x,
	y,
	z,
};

export struct MassPeak_t
{
	constexpr MassPeak_t() noexcept {}
	constexpr MassPeak_t(std::floating_point auto v) noexcept : m_Value(v) {}

	double m_Value = 0;
	bool m_Identified = false;
	IonType m_Type = IonType::UNKNOWN_TYPE;
	unsigned short m_Count = 0;

	constexpr auto operator<=> (const MassPeak_t& rhs) noexcept { return m_Value <=> rhs.m_Value; }
	constexpr auto operator<=> (std::floating_point auto rhs) noexcept { return m_Value <=> rhs; }
	constexpr bool operator== (const MassPeak_t& rhs) noexcept { return gcem::abs(m_Value - rhs.m_Value) <= DBL_EPSILON; }
};

// Y1 ion.
export AminoAcids_e IdentifyCTerminal(const vector<double>& rgflMassData, double M_plus_1) noexcept
{
	for (const auto& mass : rgflMassData)
	{
		switch ((int)std::round(M_plus_1 - mass))
		{
		case 146:	// residue mass + H2O (Consider as natural loss).
			std::cout << std::format("{} is the b[n-1] ion and the C-terminal amino acid is {}.\n", mass, g_rgAminoAcidsData[Lysine].m_3Letters);
			return Lysine;
		case 174:
			std::cout << std::format("{} is the b[n-1] ion and the C-terminal amino acid is {}.\n", mass, g_rgAminoAcidsData[Arginine].m_3Letters);
			return Arginine;
		default:
			continue;
		}
	}

	std::cout << std::format("b[n-1] ion no found.\n");
	return NOT_AN_AMINO_ACID;
}