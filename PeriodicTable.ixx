module;

#include <algorithm>

export module PeriodicTable;

import UtlString;

export namespace Monoisotopic
{
	constexpr double Hydrogen = 1.007825031898;	// H-1

	constexpr double Lithium = 7.016003434;	// Li-7
	constexpr double Beryllium = 9.01218306;	// Be-9
	constexpr double Boron = 11.009305167;	// B-11
	constexpr double Carbon = 12;	// C-12
	constexpr double Nitrogen = 14.003074004;	// N-14
	constexpr double Oxygen = 15.99491461956;	// O-16
	constexpr double Fluorine = 18.998403162067;	// F-19

	constexpr double Sodium = 22.9897692820;	// Na-23
	constexpr double Magnesium = 23.985041697;	// Mg-24
	constexpr double Aluminium = 26.98153841;	// Al-27
	constexpr double Silicon = 27.9769265350;	// Si-28
	constexpr double Phosphorus = 30.9737619986;	// P-31
	constexpr double Sulfur = 31.972071;	// S-32
	constexpr double Chlorine = 34.96885269;	// Cl-35

	constexpr double Iodine = 126.904473;	// I-127

	template<StringLiteral STR, char cLast = '\0', double flAccumulatedWt = 0.0, double flAtomWt = 0.0, uint32 iAtomCt = 0>
	requires (STR.size >= 1)
	consteval double _impl_MWt(void) noexcept
	{
		using Shortened = StringLiteral<char, STR.size - 1>;

#define SIMP_ELEM(elem)	else if constexpr (STR[0] == #elem[0])	\
							return _impl_MWt<Shortened(&STR[1]), STR[0], flAccumulatedWt + flAtomWt * std::max(1U, iAtomCt), elem, 0>()
#define COMP_ELEM(Upper, Lower, elem)	else if constexpr (cLast == Upper && STR[0] == Lower)	\
											return _impl_MWt<Shortened(&STR[1]), STR[0], flAccumulatedWt, elem, 0>()

		// Error handle.
		if constexpr (STR[0] && !isdigit_c(STR[0]) && !islower_c(STR[0]) && !isupper_c(STR[0]))
			return -1.0;

		// Period 1
		SIMP_ELEM(Hydrogen);

		// Period 2
		COMP_ELEM('L', 'i', Lithium);
		COMP_ELEM('B', 'e', Beryllium);
		SIMP_ELEM(Boron);
		SIMP_ELEM(Carbon);
		SIMP_ELEM(Nitrogen);
		SIMP_ELEM(Oxygen);
		SIMP_ELEM(Fluorine);

		// Period 3
		COMP_ELEM('N', 'a', Sodium);
		COMP_ELEM('M', 'g', Magnesium);
		COMP_ELEM('A', 'l', Aluminium);
		COMP_ELEM('S', 'i', Silicon);
		SIMP_ELEM(Phosphorus);
		SIMP_ELEM(Sulfur);
		COMP_ELEM('C', 'l', Chlorine);

		// Period 5
		SIMP_ELEM(Iodine);

#undef SIMP_ELEM
#undef COMP_ELEM

		// Digit handling.
		else if constexpr (isdigit_c(STR[0]))
			return _impl_MWt<Shortened(&STR[1]), STR[0], flAccumulatedWt, flAtomWt, iAtomCt * 10 + (STR[0] - '0')>();

		// General case of letters
		else if constexpr (isupper_c(STR[0]))
			return _impl_MWt<Shortened(&STR[1]), STR[0], flAccumulatedWt + flAtomWt * std::max(1U, iAtomCt), 0.0, 0>();
		else if constexpr (islower_c(STR[0]))
			return _impl_MWt<Shortened(&STR[1]), STR[0], flAccumulatedWt, 0.0, 0>();

		// Final.
		else
			return flAccumulatedWt + flAtomWt * std::max(1U, iAtomCt);
	}

	template<StringLiteral MolecularFormula>
	constexpr double MWt = _impl_MWt<MolecularFormula>();
};
