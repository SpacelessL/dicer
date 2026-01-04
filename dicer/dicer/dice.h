#pragma once

#include "symbol_set.h"
#include "polynomial.h"
#include "sus/math.hpp"

#include <variant>

namespace spaceless {

template<SymbolSet S, std::signed_integral T = int8_t>
class dice {
public:
	using symbol_set = S;
	static constexpr auto N = symbol_set::N;
	using mono = monomial<N, T>;
	using poly = polynomial<N, T>;

	dice() = default;
	dice(const dice &) = default;
	dice(dice &&) noexcept = default;
	dice &operator = (const dice &) = default;
	dice &operator = (dice &&) noexcept = default;

	dice(std::shared_ptr<symbol_set> s, poly pl) : s(std::move(s)), pl(pl.terms().size() ? std::move(pl) : poly::one()) { normalize(); }

	struct statistics {
		int min = std::numeric_limits<int>::max(), max = std::numeric_limits<int>::min();
		double E = 0, V = 0, D = 0, skewness = 0;
	};

	const symbol_set &symbols() const { return *s; }

	const poly &generating_function() const { return pl; }

	static dice zero(std::shared_ptr<symbol_set> s) {
		return{ std::move(s), poly::identity() };
	}

	dice &clamp(double eps = spaceless::eps) {
		pl.trim(eps);
		normalize();
		return *this;
	}
	[[nodiscard]]
	dice clamped(double eps = spaceless::eps) const {
		dice ret = *this;
		ret.clamp(eps);
		return ret;
	}

	void normalize() {
		accumulator<> sum;
		for (const auto &coef : pl.terms() | std::views::values)
			sum += coef;
		for (auto &coef : pl.terms() | std::views::values)
			coef /= sum;
	}

	// merge range of <distribution, Prob> into a new distribution
	template<std::ranges::input_range Range>
	static dice merge_dice(const Range &dice_with_prob) {
		unordered_map<mono, accumulator<>> map;
		dice ret;
		for (auto &&[dice, prob] : dice_with_prob) {
			for (auto &&[mn, coef] : dice.pl.terms)
				map[mn] += coef * prob;
			ret.s = dice.s;
		}
		ret.pl.terms().reserve(map.size());
		for (auto &&[mn, prob] : map)
			ret.pl.terms().emplace(mn, prob);
		ret.normalize();
		return ret;
	}

	dice operator + (const dice &rhs) const &{
		dice ret(s, pl * rhs.pl);
		ret.normalize();
		return ret;
	}
	dice operator + (dice &&rhs) const &{ return std::move(rhs += *this); }
	dice operator + (const dice &rhs) &&{ return std::move(*this += rhs); }
	dice &operator += (const dice &rhs) {
		pl *= rhs.pl;
		normalize();
		return *this;
	}

	dice operator * (int times) const {
		dice ret(s, pl.power(times));
		ret.normalize();
		return ret;
	}
	friend dice operator * (int times, const dice &d) { return d * times; }

	dice negation() const {
		dice ret;
		ret.s = s;
		for (auto &&[mn, p] : pl.terms())
			ret.pl[mn.inverse()] = p;
		return ret;
	}

	template<typename F, typename NS = S>
	auto transformed(F &&func, std::shared_ptr<NS> ns = {}) const
	requires requires(const mono &mn) { { func(mn) } -> MonomialType; } {
		using new_mono = std::remove_cvref_t<std::invoke_result_t<F, mono>>;
		if constexpr (std::is_same_v<NS, S>) if (!ns) ns = s;
		ENSURE(ns, "Missing new symbol set");
		unordered_map<new_mono, accumulator<>> terms;
		for (auto &&[mono, coef] : pl.terms())
			terms[func(mono)] += coef;
		return distribution<NS, typename new_mono::underlying_type>(std::move(ns), std::move(terms));
	}

	auto sample() const {
		double s = get_random();
		accumulator<> sum;
		for (auto &&[mn, coef] : pl.terms())
			if (double(sum += coef) - s >= -eps)
				return mn;
		// shouldn't be able to reach here
		return pl.terms().begin()->first;
	}

	auto get_statistics(int index) const {
		statistics ret;
		accumulator<> E, E2, E3;
		for (auto &&[mn, coef] : pl.terms()) {
			int value = mn[index];
			ret.min = std::min(ret.min, value);
			ret.max = std::max(ret.max, value);
			E += value * coef;
			E2 += squared(value) * coef;
			E3 += cubed(value) * coef;
		}
		ret.E = E;
		ret.V = E2 - squared(E);
		ret.D = std::sqrt(ret.V);
		ret.skewness = (E3 - 3 * E * ret.V - cubed(E)) / cubed(ret.D);
		return ret;
	}

private:
	std::shared_ptr<symbol_set> s;
	poly pl;
};

template<SymbolSet S, std::signed_integral T = int8_t>
auto make_simple_dice(std::shared_ptr<S> s, std::string_view text) requires std::is_same_v<typename S::symbol, char> {
	polynomial<S::N, T> pl;
	for (auto word : text | std::views::split(',')) {
		monomial<S::N, T> mn;
		auto get_num = [](std::string_view num) -> std::optional<T> {
			T ret = 0;
			auto ed = num.data() + num.length();
			auto [ptr, ec] = std::from_chars(num.data(), ed, ret);
			if (ec == std::errc{} && ptr == ed)
				return ret;
			return{};
		};
		auto process = [&](std::string_view part) {
			if (auto idx = s->get_index(part.back()); idx != -1) {
				if (part.length() == 1)
					mn.add_exponent(idx, 1);
				else {
					auto num = get_num(part.substr(0, part.length() - 1));
					ENSURE(num);
					mn.add_exponent(idx, num.value());
				}
			}
			else {
				ENSURE(S::N == 1);
				auto num = get_num(part);
				ENSURE(num);
				mn.add_exponent(0, num.value());
			}
		};
		for (int i = 0, j = -1; i < word.size(); i++)
			if (s->contains(word[i]) || i == word.size() - 1) {
				process(std::string_view(word).substr(j + 1, i - j));
				j = i;
			}
		pl += mn;
	}
	return dice(s, pl);
}

}
