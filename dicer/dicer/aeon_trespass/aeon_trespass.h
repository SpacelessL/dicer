#pragma once

#include "dicer.h"
#include <map>

namespace spaceless::aeon_trespass {

enum class power_dice_symbol : uint8_t {
	power = 0,
	potential,
	dot,
	count,
};

static const char *ToString(power_dice_symbol symbol) {
	switch (symbol)
	{
	case power_dice_symbol::power:
		return "<";
	case power_dice_symbol::potential:
		return "^";
	case power_dice_symbol::dot:
		return ".";
	default:
		return "UNKNOWN";
	}
}

using power_dice = dice<power_dice_symbol>;

static const auto red_power_dice = make_dice<power_dice_symbol>(
	make_dice_face(power_dice_symbol::power, 2, power_dice_symbol::potential, 0, power_dice_symbol::dot, 1),
	make_dice_face(power_dice_symbol::power, 1, power_dice_symbol::potential, 1, power_dice_symbol::dot, 0),
	make_dice_face(power_dice_symbol::power, 1, power_dice_symbol::potential, 1, power_dice_symbol::dot, 0),
	make_dice_face(power_dice_symbol::power, 1, power_dice_symbol::potential, 0, power_dice_symbol::dot, 0),
	make_dice_face(power_dice_symbol::power, 0, power_dice_symbol::potential, 1, power_dice_symbol::dot, 0),
	make_dice_face(power_dice_symbol::power, 0, power_dice_symbol::potential, 0, power_dice_symbol::dot, 1)
);

static const auto black_power_dice = make_dice<power_dice_symbol>(
	make_dice_face(power_dice_symbol::power, 2, power_dice_symbol::potential, 2, power_dice_symbol::dot, 1),
	make_dice_face(power_dice_symbol::power, 2, power_dice_symbol::potential, 1, power_dice_symbol::dot, 0),
	make_dice_face(power_dice_symbol::power, 2, power_dice_symbol::potential, 0, power_dice_symbol::dot, 1),
	make_dice_face(power_dice_symbol::power, 1, power_dice_symbol::potential, 1, power_dice_symbol::dot, 0),
	make_dice_face(power_dice_symbol::power, 1, power_dice_symbol::potential, 1, power_dice_symbol::dot, 0),
	make_dice_face(power_dice_symbol::power, 0, power_dice_symbol::potential, 1, power_dice_symbol::dot, 1)
);

static const auto white_power_dice = make_dice<power_dice_symbol>(
	make_dice_face(power_dice_symbol::power, 3, power_dice_symbol::potential, 2, power_dice_symbol::dot, 1),
	make_dice_face(power_dice_symbol::power, 2, power_dice_symbol::potential, 3, power_dice_symbol::dot, 0),
	make_dice_face(power_dice_symbol::power, 2, power_dice_symbol::potential, 1, power_dice_symbol::dot, 1),
	make_dice_face(power_dice_symbol::power, 1, power_dice_symbol::potential, 3, power_dice_symbol::dot, 0),
	make_dice_face(power_dice_symbol::power, 1, power_dice_symbol::potential, 2, power_dice_symbol::dot, 1),
	make_dice_face(power_dice_symbol::power, 1, power_dice_symbol::potential, 1, power_dice_symbol::dot, 1)
);

static const auto mortal_power_dice = make_dice<power_dice_symbol>(
	make_dice_face(power_dice_symbol::power, 3, power_dice_symbol::potential, 4, power_dice_symbol::dot, 1),
	make_dice_face(power_dice_symbol::power, 3, power_dice_symbol::potential, 2, power_dice_symbol::dot, 0),
	make_dice_face(power_dice_symbol::power, 2, power_dice_symbol::potential, 4, power_dice_symbol::dot, 0),
	make_dice_face(power_dice_symbol::power, 2, power_dice_symbol::potential, 3, power_dice_symbol::dot, 1),
	make_dice_face(power_dice_symbol::power, 1, power_dice_symbol::potential, 3, power_dice_symbol::dot, 1),
	make_dice_face(power_dice_symbol::power, 1, power_dice_symbol::potential, 2, power_dice_symbol::dot, 1)
);

struct kratos_pool {
	int break_count = 0;
	int opening_count = 0;
	int black_count = 0;
	int hope_count = 0;
};

struct reroll_profile {
	int power_reroll = 0, hit_reroll = 0;
};

struct weapon_profile {
	int perhit_red = 0, perhit_black = 0, perhit_white = 0;
	int hit_dice_count = 0;
	int hit_bonus = 0, power_bonus = 0;
};

struct titan_profile {
	int red = 0, black = 0, white = 0;
};

struct primordial_profile {
	std::vector<std::pair<int, double>> to_at_field_with_prob() const {
		unordered_map<int, int> m;
		for (int at : at_fields)
			m[at]++;
		std::vector<std::pair<int, double>> ret;
		for (auto &&[at, count] : m)
			ret.emplace_back(at, double(count) / at_fields.size());
		return ret;
	}

	std::vector<int> at_fields;
	int to_hit = 0;
};

static number_dice power_distribution(int red, int black, int white, int power_bonus, int break_count, int hope_count) {
	auto dist = red * red_power_dice + black * black_power_dice + white * white_power_dice;
	return dist.transformed([break_count, hope_count](const dice_face<power_dice_symbol> &face) {
		int power = count_symbol(face, power_dice_symbol::power);
		int potential = count_symbol(face, power_dice_symbol::potential);
		int dot = count_symbol(face, power_dice_symbol::dot);
		// process break
		int converted_break = std::min(potential, break_count);
		potential -= converted_break;
		// process hope
		int converted_hope = std::min(potential + dot, hope_count);
		return make_dice_face(power_dice_symbol::power, power + converted_break + converted_hope);
	}).to_number_dice(power_dice_symbol::power) + power_bonus;
}

static number_dice potential_distribution(int red, int black, int white, int bonus) {
	auto dist = red * red_power_dice.to_number_dice(power_dice_symbol::potential) + black * black_power_dice.to_number_dice(power_dice_symbol::potential) + white * white_power_dice.to_number_dice(power_dice_symbol::potential);
	return dist + bonus;
}

static number_dice armor_distribution(int red, int black, int white, int bonus) {
	return potential_distribution(red, black, white, bonus);
}

static number_dice hit_count_distribution(int count, int to_hit, int bonus, int reroll_count) {
	std::array<int, 10> numbers;
	for (int i = 1; i <= 10; i++)
		if (i == 1)
			numbers[i - 1] = 0;
		else if (i == 10)
			numbers[i - 1] = 1;
		else
			numbers[i - 1] = i + bonus >= to_hit ? 1 : 0;
	auto ori_dist = make_number_dice(numbers.begin(), numbers.end());
	auto dist = ori_dist * count;
	std::vector<std::pair<number_dice, double>> vec;
	for (auto &&[face, prob] : dist.to_faces_with_prob()) {
		int a = count_symbol(face);
		int b = std::min(reroll_count, count - a);
		vec.emplace_back(ori_dist * b + a, prob);
	}
	return number_dice::merge_dice(vec);
}
// power distribution
static number_dice titan_attack(titan_profile titan, weapon_profile weapon, kratos_pool kratos_pool, reroll_profile reroll, int to_hit) {
	auto hit = hit_count_distribution(weapon.hit_dice_count, to_hit, weapon.hit_bonus + kratos_pool.opening_count, reroll.hit_reroll);
	std::vector<std::pair<number_dice, double>> vec;
	for (auto &&[face, prob] : hit.to_faces_with_prob()) {
		int hit_count = count_symbol(face);
		int red = titan.red + hit_count * weapon.perhit_red;
		int black = titan.black + hit_count * weapon.perhit_black;
		int white = titan.white + hit_count * weapon.perhit_white;
		if (hit_count == 0)
			vec.emplace_back(make_number_dice(0), prob);
		else
			vec.emplace_back(power_distribution(red, black, white, weapon.power_bonus, kratos_pool.break_count, kratos_pool.hope_count), prob);
	}
	return number_dice::merge_dice(vec);
}

static double titan_attack_primordial(titan_profile titan, weapon_profile weapon, kratos_pool kratos_pool, reroll_profile reroll, primordial_profile primordial) {
	auto hit = hit_count_distribution(weapon.hit_dice_count, primordial.to_hit, weapon.hit_bonus + kratos_pool.opening_count, reroll.hit_reroll);
	static const auto combination_indexer = power_dice::get_combination_indexer(make_reference_range(red_power_dice, black_power_dice, white_power_dice), true);
	static const auto red_comb_dice = red_power_dice.to_combination_dice(combination_indexer);
	static const auto black_comb_dice = black_power_dice.to_combination_dice(combination_indexer);
	static const auto white_comb_dice = white_power_dice.to_combination_dice(combination_indexer);
	double ret = 0;
	unordered_map<int, double> succ_for_at;
	for (auto &&[face, prob] : hit.to_faces_with_prob()) {
		int hit_count = count_symbol(face);
		int red = titan.red + hit_count * weapon.perhit_red;
		int black = titan.black + hit_count * weapon.perhit_black;
		int white = titan.white + hit_count * weapon.perhit_white;
		auto comb = red * red_comb_dice + black * black_comb_dice + white * white_comb_dice;
		for (auto &&[at, at_prob] : primordial.to_at_field_with_prob()) {
			std::map<std::tuple<int, int, int, int, int>, number_dice> succ_prob_cache;
			auto succ_prob = [&succ_prob_cache](int red, int black, int white, int break_count, int hope_count, int at) {
				number_dice dist;
				auto it = succ_prob_cache.find(std::make_tuple(red, black, white, break_count, hope_count));
				if (it != succ_prob_cache.end())
					dist = it->second;
				else {
					dist = power_distribution(red, black, white, 0, break_count, hope_count);
					succ_prob_cache.emplace(std::make_tuple(red, black, white, break_count, hope_count), dist);
				}
				accumulator ret;
				for (auto &&[face, prob] : dist.get_poly().terms)
					if (face[0] >= at)
						ret += prob;
				return double(ret);
			};
			for (auto &&[comb_face, comb_prob] : comb.to_faces_with_prob()) {
				// dfs comb_face, reroll reroll.power_reroll at max, calculate remain break & hope & at, get max succ_prob, add max_succ_prob * comb_prob * prob to final result
				double max_succ = 0;
				std::vector current_combo_face(comb_face.begin(), comb_face.end());
				auto dfs = [&](int idx, int remain_reroll, int reroll_red, int reroll_black, int reroll_white, auto &&dfs) {
					if (idx == current_combo_face.size()) {
						int power = 0, potential = 0, dot = 0;
						for (auto &&[comb_face, count] : current_combo_face) {
							auto face = comb_face.first;
							power += count_symbol(face, power_dice_symbol::power) * count;
							potential += count_symbol(face, power_dice_symbol::potential) * count;
							dot += count_symbol(face, power_dice_symbol::dot) * count;
						}
						int converted_break = std::min(potential, kratos_pool.break_count);
						potential -= converted_break;
						int converted_hope = std::min(potential + dot, kratos_pool.hope_count);

						int remain_break = kratos_pool.break_count - converted_break, remain_hope = kratos_pool.hope_count - converted_hope;
						int current_power = weapon.power_bonus + power + converted_break + converted_hope;
						
						double succ = succ_prob(reroll_red, reroll_black, reroll_white, remain_break, remain_hope, at - current_power);
						if (succ > max_succ)
							max_succ = succ;
						return;
					}
					for (int i = 0; i <= std::min(remain_reroll, current_combo_face[idx].second); i++) {
						auto current_dice = current_combo_face[idx].first.second;
						current_combo_face[idx].second -= i;
						dfs(idx + 1, remain_reroll - i, reroll_red + (current_dice == &red_power_dice ? i : 0)
							, reroll_black + (current_dice == &black_power_dice ? i : 0)
							, reroll_white + (current_dice == &white_power_dice ? i : 0), dfs);
						current_combo_face[idx].second += i;
					}
				};
				dfs(0, reroll.power_reroll, 0, 0, 0, dfs);
				ret += max_succ * comb_prob * at_prob * prob;
				succ_for_at[at] += max_succ * comb_prob * prob;
			}
		}
	}
	for (auto &&[at, prob] : succ_for_at)
		std::cout << "succ rate for AT " << at << " is " << prob << std::endl;
	return ret;
}

}
