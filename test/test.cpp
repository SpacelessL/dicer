#include "mimalloc-new-delete.h"
#include "dicer.h"
#include "dicer/aeon_trespass/aeon_trespass.h"
#include <iostream>
#include <mutex>
#include <ctime>
#include "omp.h"
#include <complex>

using namespace spaceless;

static void AeonTrespassTest() {
	using namespace aeon_trespass;
	titan_profile titan;
	titan.red = 1;
	weapon_profile weapon;
	weapon.hit_bonus = 1;
	weapon.hit_dice_count = 3;
	weapon.perhit_red = 2;
	weapon.perhit_black = 0;
	weapon.power_bonus = 1;
	kratos_pool pool;
	pool.break_count = 2;
	pool.opening_count = 2;
	pool.hope_count = 0;
	reroll_profile reroll;
	reroll.hit_reroll = 3;
	reroll.power_reroll = 2;
	primordial_profile primordial;
	primordial.at_fields = { 8, 8, 8, 9, 10, 11 }; // Bird Lv III
	primordial.to_hit = 7;
	auto prob = titan_attack_primordial(titan, weapon, pool, reroll, primordial);
	std::cout << "overall succ prob for Lv III Brother Bird is " << prob << std::endl;
}

static void ExplodeTest() {
	auto dist = d6.explode([](const number_dice_face &face) -> std::optional<number_dice> {
		if (count_symbol(face) == 1)
			return d6;
		return std::nullopt;
	}, 10);
	dist.debug_log();
}

static void PolynomialTest() {
	{
		polynomial a, b;
		a.terms.emplace(monomial({ 2, 0 }), 1);
		a.terms.emplace(monomial({ 2, 1 }), 2);
		a.terms.emplace(monomial({ 1, 1 }), -2);
		a.terms.emplace(monomial({ 1, 2 }), 2);
		a.terms.emplace(monomial({ 0, 2 }), -3);

		b.terms.emplace(monomial({ 1, 0 }), 1);
		b.terms.emplace(monomial({ 0, 1 }), 1);

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
	{
		polynomial a, b;
		a.terms.emplace(monomial({ 3 }), 1);
		a.terms.emplace(monomial({ 2 }), -2);
		a.terms.emplace(monomial({ 0 }), -4);

		b.terms.emplace(monomial({ 1 }), 1);
		b.terms.emplace(monomial({ 0 }), -3);

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

static void Playground() {
	auto indexer = number_dice::get_combination_indexer(make_reference_range(d6), true);
	auto dd = d6.to_combination_dice(indexer);
	auto output = [&](combination_dice_face_symbol<number_dice> symbol) {
		std::string ret;
		for (auto &&[s, count] : symbol.first) {
			ret += std::to_string(count) + " " + std::to_string((unsigned long long)((void **)symbol.second));
		}
		return ret;
	};
	std::cout << (unsigned long long)((void **)&d6) << std::endl;
	dd.debug_log(output);
}

int main(int argc, char **argv) {
	int begin = clock(), end = 0;

	for (int i = 0; i < 10; i++)
		AeonTrespassTest();
	//ExplodeTest();
	//PolynomialTest();
	//Playground();

	end = clock();
	std::cout << "cost time : " << (end - begin) / 1000.0 << "s" << std::endl;
	return 0;
}
