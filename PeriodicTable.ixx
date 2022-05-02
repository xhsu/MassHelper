module;

#include <algorithm>
#include <array>
#include <tuple>

export module PeriodicTable;

import UtlString;

namespace Isotope
{
	using std::array;
	using std::pair;
	using std::tuple;

	// Period 1
	constexpr auto Hydrogen = array
	{
		std::make_tuple(1.0078250322, 0.99972, 0.99999),
		std::make_tuple(2.0141017781, 0.00001, 0.00028),
	};
	constexpr auto Helium = array
	{
		std::make_pair(3.016029321967, 0.000002),
		std::make_pair(4.002603254130, 0.999998),
	};

	// Period 2
	constexpr auto Lithium = array
	{
		std::make_tuple(6.01512289, 0.019, 0.078),
		std::make_tuple(7.01600344, 0.922, 0.981),
	};
	constexpr auto Beryllium = array
	{
		std::make_pair(9.0121831, 1.0),
	};
	constexpr auto Boron = array
	{
		std::make_tuple(10.0129369, 0.189, 0.204),
		std::make_tuple(11.00930517, 0.796, 0.811),
	};
	constexpr auto Carbon = array
	{
		std::make_tuple(12.0, 0.9884, 0.9904),
		std::make_tuple(13.003354835, 0.0096, 0.0116),
	};
	constexpr auto Nitrogen = array
	{
		std::make_tuple(14.003074004, 0.99578, 0.99663),
		std::make_tuple(15.000108899, 0.00337, 0.00422),
	};
	constexpr auto Oxygen = array
	{
		std::make_tuple(15.994914619, 0.99738, 0.99776),
		std::make_tuple(16.999131757, 0.000367, 0.000400),
		std::make_tuple(17.999159613, 0.00187, 0.00222),
	};
	constexpr auto Fluorine = array
	{
		std::make_pair(18.998403162, 1.0),
	};
	constexpr auto Neon = array
	{
		std::make_pair(19.99244018, 0.9048),
		std::make_pair(20.9938467, 0.0027),
		std::make_pair(21.9913851, 0.0925),
	};

	// Period 3
	constexpr auto Sodium = array
	{
		std::make_pair(22.98976928, 1.0),
	};
	constexpr auto Magnesium = array
	{
		std::make_tuple(23.98504170, 0.7888, 0.7905),
		std::make_tuple(24.9858370, 0.09988, 0.10034),
		std::make_tuple(25.9825930, 0.1096, 0.1109),
	};
	constexpr auto Aluminium = array
	{
		std::make_pair(26.9815384, 1.0),
	};
	constexpr auto Silicon = array
	{
		std::make_tuple(27.976926535, 0.92191, 0.92318),
		std::make_tuple(28.976494665, 0.04645, 0.04699),
		std::make_tuple(29.9737701, 0.03037, 0.03110),
	};
	constexpr auto Phosphorus = array
	{
		std::make_pair(30.973761998, 1.0),
	};
	constexpr auto Sulfur = array
	{
		std::make_tuple(31.972071174, 0.9441, 0.9529),
		std::make_tuple(32.97145891, 0.00729, 0.00797),
		std::make_tuple(33.9678670, 0.0396, 0.0477),
		std::make_tuple(35.967081, 0.000129, 0.000187),
	};
	constexpr auto Chlorine = array
	{
		std::make_tuple(34.9688527, 0.755, 0.761),
		std::make_tuple(36.9659026, 0.239, 0.245),
	};
	constexpr auto Argon = array	// #FIXME Argon is a special case, for it has a 100% situation.
	{
		std::make_tuple(35.9675451, 0.0000, 0.0207),
		std::make_tuple(37.962732, 0.000, 0.043),
		std::make_tuple(39.96238312, 0.936, 1.000),
	};
}

export namespace ArStd
{
	using std::array;
	using std::pair;
	using std::tuple;

	template<size_t N>
	consteval double _impl_CalcRelativeAtomicMass(const array<pair<double, double>, N>& rgfl) noexcept
	{
		double ret = 0;

		for (const auto& [mass, abundance] : rgfl)
			ret += mass * abundance;

		return ret;
	}

	template<size_t N> requires(N > 1)
	consteval pair<double, double> _impl_CalcRelativeAtomicMass(const array<tuple<double, double, double>, N>& rgfl) noexcept
	{
		enum
		{
			MASS,
			LOW,
			HIGH,
		};

		return [&]<size_t... I>(std::index_sequence<I...>) -> pair<double, double>
		{
			return std::make_pair
			(
				std::get<MASS>(rgfl[0]) * (1.0 - (std::get<LOW>(rgfl[I + 1]) + ...))	// Initial value. Assuming the lightest isotope is at maxium.
				+ ((std::get<MASS>(rgfl[I + 1]) * std::get<LOW>(rgfl[I + 1])) + ...),	// Everything else is at minium.

				std::get<MASS>(rgfl[0]) * (1.0 - (std::get<HIGH>(rgfl[I + 1]) + ...))	// This is just the inversed case of 'min'. Assuming the lightest isotope is at minium.
				+ ((std::get<MASS>(rgfl[I + 1]) * std::get<HIGH>(rgfl[I + 1])) + ...)	// Everything else is at maxium.
			);
		}
		(std::make_index_sequence<N - 1>{});
	}

	// Period 1
	constexpr auto Hydrogen = _impl_CalcRelativeAtomicMass(Isotope::Hydrogen);
	constexpr auto Helium = _impl_CalcRelativeAtomicMass(Isotope::Helium);

	// Period 2
	constexpr auto Lithium = _impl_CalcRelativeAtomicMass(Isotope::Lithium);
	constexpr auto Beryllium = _impl_CalcRelativeAtomicMass(Isotope::Beryllium);
	constexpr auto Boron = _impl_CalcRelativeAtomicMass(Isotope::Boron);
	constexpr auto Carbon = _impl_CalcRelativeAtomicMass(Isotope::Carbon);
	constexpr auto Nitrogen = _impl_CalcRelativeAtomicMass(Isotope::Nitrogen);
	constexpr auto Oxygen = _impl_CalcRelativeAtomicMass(Isotope::Oxygen);
	constexpr auto Fluorine = _impl_CalcRelativeAtomicMass(Isotope::Fluorine);
	constexpr auto Neon = _impl_CalcRelativeAtomicMass(Isotope::Neon);

	// Period 3
	constexpr auto Sulfur = _impl_CalcRelativeAtomicMass(Isotope::Sulfur);
	constexpr auto Chlorine = _impl_CalcRelativeAtomicMass(Isotope::Chlorine);
	constexpr auto Argon = _impl_CalcRelativeAtomicMass(Isotope::Argon);
}

namespace Monoisotopic
{
	export constexpr double Hydrogen = 1.007825031898;	// H-1

	export constexpr double Lithium = 7.016003434;	// Li-7
	export constexpr double Beryllium = 9.01218306;	// Be-9
	export constexpr double Boron = 11.009305167;	// B-11
	export constexpr double Carbon = 12;	// C-12
	export constexpr double Nitrogen = 14.003074004;	// N-14
	export constexpr double Oxygen = 15.99491461956;	// O-16
	export constexpr double Fluorine = 18.998403162067;	// F-19

	export constexpr double Sodium = 22.9897692820;	// Na-23
	export constexpr double Magnesium = 23.985041697;	// Mg-24
	export constexpr double Aluminium = 26.98153841;	// Al-27
	export constexpr double Silicon = 27.9769265350;	// Si-28
	export constexpr double Phosphorus = 30.9737619986;	// P-31
	export constexpr double Sulfur = 31.972071;	// S-32
	export constexpr double Chlorine = 34.96885269;	// Cl-35

	export constexpr double Iodine = 126.904473;	// I-127

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

	export template<StringLiteral MolecularFormula>
	constexpr double MWt = _impl_MWt<MolecularFormula>();
};
