#include "dicer.h"
#include <iostream>
#include <mutex>
#include <ctime>
#include <complex>

#include "sus.hpp"

using namespace spaceless;

static void DiceTest() {
	std::cout << "=== Basic Dice Test ===\n";

	// Create a simple d6
	auto symbols = make_simple_static_symbol_set<"N">();
	auto d6 = make_simple_dice(symbols, "1,2,3,4,5,6");

	std::cout << "d6 generating function: " << d6.generating_function() << '\n';

	// Roll 2d6
	auto two_d6 = d6 * 2;
	std::cout << "2d6 generating function: " << two_d6.generating_function() << '\n';

	// Get statistics
	auto stats = two_d6.get_statistics(0);
	std::cout << "2d6 stats: min=" << stats.min << ", max=" << stats.max
		<< ", E=" << stats.E << ", V=" << stats.V << ", D=" << stats.D << '\n';

	// Sample a few rolls
	std::cout << "Sample rolls: ";
	for (int i = 0; i < 5; ++i)
		std::cout << two_d6.sample()[0] << " ";
	std::cout << '\n';
}

static void SymbolDiceTest() {
	std::cout << "\n=== Symbol Dice Test ===\n";

	// Create combat dice with multiple symbol types
	auto combat_symbols = make_simple_static_symbol_set<"SMD">();
	auto combat_die = make_simple_dice(combat_symbols, "0,0,SS,SD,S,MM,M,D");

	std::cout << "Combat die generating function: " << combat_die.generating_function() << '\n';

	// Roll 3 combat dice
	auto three_combat = combat_die * 3;

	// Get statistics for each symbol
	auto sword_stats = three_combat.get_statistics(0);
	auto magic_stats = three_combat.get_statistics(1);
	auto defense_stats = three_combat.get_statistics(2);

	std::cout << "3 combat dice - E[swords]=" << sword_stats.E
		<< ", E[magic]=" << magic_stats.E
		<< ", E[defense]=" << defense_stats.E << '\n';

	// Sample a roll
	auto roll = three_combat.sample();
	std::cout << "Sample roll: " << roll[0] << " swords, "
		<< roll[1] << " magic, " << roll[2] << " defense\n";
}

static void CombinationDiceTest() {
	std::cout << "\n=== Combination Dice Test ===\n";

	// Create a simple d6
	auto symbols = make_simple_static_symbol_set<"N">();
	auto d6 = make_simple_dice(symbols, "1,2,3,4,5,6");

	// Create combination symbol set
	auto comb_symbols = make_combination_symbol_set(std::array{d6});
	std::cout << "Combination symbol set size: " << comb_symbols->size() << '\n';

	// Convert to combination dice
	auto comb_d6 = d6.to_combination_dice(comb_symbols);

	// Roll 10d6 as combination dice (tracks which faces appeared)
	auto ten_d6 = comb_d6 * 10;

	std::cout << "10d6 combination dice - number of outcomes: "
		<< ten_d6.generating_function().terms().size() << '\n';

	// Sample and show face distribution
	auto roll = ten_d6.sample();
	std::cout << "Sample roll (face counts): ";
	for (int i = 0; i < comb_symbols->size(); ++i) {
		int count = roll[i];
		if (count > 0) {
			auto face_value = comb_symbols->from_index(i).face[0];
			std::cout << count << "x" << face_value << " ";
		}
	}
	std::cout << '\n';
}

static void PolynomialTest() {
	{
		polynomial<DYNAMIC> a, b;
		a += polynomial(monomial<DYNAMIC>({ 2 }), 1);
		a += polynomial(monomial<DYNAMIC>({ 2, 1 }), 2);
		a += polynomial(monomial<DYNAMIC>({ 1, 1 }), -2);
		a += polynomial(monomial<DYNAMIC>({ 1, 2 }), 2);
		a += polynomial(monomial<DYNAMIC>({ 0, 2 }), -3);

		b += polynomial(monomial<DYNAMIC>({ 1 }), 1);
		b += polynomial(monomial<DYNAMIC>({ 0, 1 }), 1);

		auto c = a / b, d = a % b;
		std::cout << a << '\n';
		std::cout << b << '\n';
		std::cout << c << '\n';
		std::cout << d << '\n';
		auto e = c * b;
		std::cout << e << '\n';
		auto f = e + d;
		std::cout << f << '\n';
		auto g = f - a;
		g.trim();
		std::cout << g << '\n';
	}
	{
		polynomial<DYNAMIC> a, b;
		a += polynomial(monomial<DYNAMIC>({ 3 }), 1);
		a += polynomial(monomial<DYNAMIC>({ 2 }), -2);
		a -= 4;

		b += polynomial(monomial<DYNAMIC>({ 1 }), 1);
		b -= 3;

		auto c = a / b, d = a % b;
		std::cout << a << '\n';
		std::cout << b << '\n';
		std::cout << c << '\n';
		std::cout << d << '\n';
		auto e = c * b;
		std::cout << e << '\n';
		auto f = e + d;
		std::cout << f << '\n';
		auto g = f - a;
		std::cout << g << '\n';
	}
}

int main(int argc, char **argv) {
	init(argc, argv);

	DiceTest();
	SymbolDiceTest();
	CombinationDiceTest();

	//PolynomialTest();

	return 0;
}
