#include "dicer.hpp"
#include <mutex>
#include <ctime>
#include <complex>

#include "sus.hpp"

using namespace spaceless;

static void DiceTest() {
	LOG(INFO)("=== Basic Dice Test ===");

	// Create a simple d6
	auto symbols = make_simple_static_symbol_set<"N">();
	auto d6 = make_simple_dice(symbols, "1,2,3,4,5,6");

	LOG(INFO)("d6 generating function: {}", d6.generating_function());

	// Roll 2d6
	auto two_d6 = d6 * 2;
	LOG(INFO)("2d6 generating function: {}", two_d6.generating_function());

	// Get statistics
	auto stats = two_d6.get_statistics(0);
	LOG(INFO)("2d6 stats: min={}, max={}, E={}, V={}, D={}",
		stats.min, stats.max, stats.E, stats.V, stats.D);

	// Sample a few rolls
	std::string rolls;
	for (int i = 0; i < 5; ++i)
		rolls += std::to_string(int(two_d6.sample()[0])) + " ";
	LOG(INFO)("Sample rolls: {}", rolls);
}

static void SymbolDiceTest() {
	LOG(INFO)("=== Symbol Dice Test ===");

	// Create combat dice with multiple symbol types
	auto combat_symbols = make_simple_static_symbol_set<"SMD">();
	auto combat_die = make_simple_dice(combat_symbols, "0,0,SS,SD,S,MM,M,D");

	LOG(INFO)("Combat die generating function: {}", combat_die.generating_function());

	// Roll 3 combat dice
	auto three_combat = combat_die * 3;

	// Get statistics for each symbol
	auto sword_stats = three_combat.get_statistics(0);
	auto magic_stats = three_combat.get_statistics(1);
	auto defense_stats = three_combat.get_statistics(2);

	LOG(INFO)("3 combat dice - E[swords]={}, E[magic]={}, E[defense]={}",
		sword_stats.E, magic_stats.E, defense_stats.E);

	// Sample a roll
	auto roll = three_combat.sample();
	LOG(INFO)("Sample roll: {} swords, {} magic, {} defense",
		int(roll[0]), int(roll[1]), int(roll[2]));
}

static void CombinationDiceTest() {
	LOG(INFO)("=== Combination Dice Test ===");

	// Create a simple d6
	auto symbols = make_simple_static_symbol_set<"N">();
	auto d6 = make_simple_dice(symbols, "1,2,3,4,5,6");

	// Create combination symbol set
	auto comb_symbols = make_combination_symbol_set(std::array{d6});
	LOG(INFO)("Combination symbol set size: {}", comb_symbols->size());

	// Convert to combination dice
	auto comb_d6 = d6.to_combination_dice(comb_symbols);

	// Roll 10d6 as combination dice (tracks which faces appeared)
	auto ten_d6 = comb_d6 * 10;

	LOG(INFO)("10d6 combination dice - number of outcomes: {}",
		ten_d6.generating_function().terms().size());

	// Sample and show face distribution
	auto roll = ten_d6.sample();
	std::string face_counts;
	for (int i = 0; i < comb_symbols->size(); ++i) {
		int count = roll.get_exponent(i);
		if (count > 0) {
			auto face_value = comb_symbols->from_index(i).face[0];
			face_counts += std::to_string(count) + "x" + std::to_string(int(face_value)) + " ";
		}
	}
	LOG(INFO)("Sample roll (face counts): {}", face_counts);
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
		LOG(INFO)("{}", a);
		LOG(INFO)("{}", b);
		LOG(INFO)("{}", c);
		LOG(INFO)("{}", d);
		auto e = c * b;
		LOG(INFO)("{}", e);
		auto f = e + d;
		LOG(INFO)("{}", f);
		auto g = f - a;
		g.trim();
		LOG(INFO)("{}", g);
	}
	{
		polynomial<DYNAMIC> a, b;
		a += polynomial(monomial<DYNAMIC>({ 3 }), 1);
		a += polynomial(monomial<DYNAMIC>({ 2 }), -2);
		a -= 4;

		b += polynomial(monomial<DYNAMIC>({ 1 }), 1);
		b -= 3;

		auto c = a / b, d = a % b;
		LOG(INFO)("{}", a);
		LOG(INFO)("{}", b);
		LOG(INFO)("{}", c);
		LOG(INFO)("{}", d);
		auto e = c * b;
		LOG(INFO)("{}", e);
		auto f = e + d;
		LOG(INFO)("{}", f);
		auto g = f - a;
		LOG(INFO)("{}", g);
	}
}

int main(int argc, char **argv) {
	init(argc, argv);
	
	DiceTest();
	SymbolDiceTest();
	CombinationDiceTest();

	PolynomialTest();

	return 0;
}
