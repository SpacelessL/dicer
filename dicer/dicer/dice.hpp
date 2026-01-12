#pragma once

#include "symbol_set.hpp"
#include "polynomial.hpp"
#include "sus/math.hpp"
#include "sus/random.hpp"
#include "sus/debug.hpp"
#include "sus/logging.hpp"

#include <variant>

namespace spaceless {

template<size_t N, std::signed_integral T>
struct combination_dice_face {
	polynomial<N, T> source;
	monomial<N, T> face;

	bool operator==(const combination_dice_face &) const = default;
};

}

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
	[[nodiscard]] dice clamped(double eps = spaceless::eps) const & {
		dice ret = *this;
		ret.clamp(eps);
		return ret;
	}
	[[nodiscard]] dice clamped(double eps = spaceless::eps) && {
		clamp(eps);
		return std::move(*this);
	}

	dice &normalize() {
		accumulator<> sum;
		for (const auto &coef : pl.terms() | std::views::values)
			sum += coef;
		for (auto &coef : pl.terms() | std::views::values)
			coef /= sum;
		return *this;
	}
	[[nodiscard]] dice normalized() const & {
		dice ret = *this;
		ret.normalize();
		return ret;
	}
	[[nodiscard]] dice normalized() && {
		normalize();
		return std::move(*this);
	}

	// merge range of <distribution, Prob> into a new distribution
	template<std::ranges::input_range Range>
	static dice merge_dice(const Range &dice_with_prob);

	dice operator + (const dice &rhs) const & {
		dice ret(s, pl * rhs.pl);
		ret.normalize();
		return ret;
	}
	dice operator + (dice &&rhs) const & { return std::move(rhs += *this); }
	dice operator + (const dice &rhs) && { return std::move(*this += rhs); }
	dice &operator += (const dice &rhs) {
		pl *= rhs.pl;
		normalize();
		return *this;
	}

	dice operator * (T times) const {
		dice ret(s, pl.power(times));
		ret.normalize();
		return ret;
	}
	friend dice operator * (T times, const dice &d) { return d * times; }

	[[nodiscard]] dice negated() const {
		dice ret;
		ret.s = s;
		for (auto &&[mn, p] : pl.terms())
			ret.pl[mn.inverse()] = p;
		return ret;
	}
	dice &negate() { return *this = negated(); }

	template<typename F, typename NS = S>
	auto transformed(F &&func, std::shared_ptr<NS> ns = {}) const
	requires requires(const mono &mn) { { func(mn) } -> MonomialType; };

	template<SymbolSet CombSymbolSet>
	auto to_combination_dice(std::shared_ptr<CombSymbolSet> comb_symbol_set) const;

	auto sample() const {
		double s = get_random(0.0, 1.0);
		accumulator<> sum;
		for (auto &&[mn, coef] : pl.terms())
			if (double(sum += coef) - s >= -eps)
				return mn;
		UNREACHABLE("Dice need to be normalized.", s, double(sum));
	}

	auto get_statistics(int index) const -> statistics;

private:
	std::shared_ptr<symbol_set> s;
	poly pl;
};

template<SymbolSet S, std::signed_integral T>
template<std::ranges::input_range Range>
dice<S, T> dice<S, T>::merge_dice(const Range &dice_with_prob) {
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

template<SymbolSet S, std::signed_integral T>
template<typename F, typename NS>
auto dice<S, T>::transformed(F &&func, std::shared_ptr<NS> ns) const
requires requires(const mono &mn) { { func(mn) } -> MonomialType; } {
	using new_mono = std::remove_cvref_t<std::invoke_result_t<F, mono>>;
	using new_dice = dice<NS, typename new_mono::underlying_type>;
	if constexpr (std::is_same_v<NS, S>) if (!ns) ns = s;
	ASSERT(ns, "Missing new symbol set");
	unordered_map<new_mono, accumulator<>> terms;
	for (auto &&[mono, coef] : pl.terms())
		terms[func(mono)] += coef;
	unordered_map<new_mono, double> new_terms;
	for (auto &&[mono, coef] : terms)
		new_terms[mono] = coef;
	return new_dice(std::move(ns), std::move(new_terms));
}

template<SymbolSet S, std::signed_integral T>
template<SymbolSet CombSymbolSet>
auto dice<S, T>::to_combination_dice(std::shared_ptr<CombSymbolSet> comb_symbol_set) const {
	bool include_source_dice = !comb_symbol_set->from_index(0).source.terms().empty();

	return transformed([this, include_source_dice, &comb_symbol_set](const mono &face) {
		poly source = include_source_dice ? pl : poly{};
		monomial<CombSymbolSet::N, T> result;
		result.add_exponent(comb_symbol_set->get_index({source, face}), 1);
		return result;
	}, comb_symbol_set);
}

template<SymbolSet S, std::signed_integral T>
auto dice<S, T>::get_statistics(int index) const -> statistics {
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

template<typename Dice>
using combination_dice_symbol = combination_dice_face<Dice::N, typename Dice::mono::underlying_type>;

template<typename Dice, SymbolSet S = dynamic_symbol_set<combination_dice_symbol<Dice>>>
using combination_dice = dice<S, typename Dice::mono::underlying_type>;

template<std::ranges::input_range Range>
auto make_combination_symbol_set(const Range &dice_range, bool include_source_dice = true) {
	using Dice = std::remove_cvref_t<std::ranges::range_value_t<Range>>;
	using symbol_type = combination_dice_symbol<Dice>;

	auto symbol_set = std::make_shared<dynamic_symbol_set<symbol_type>>();

	for (auto &&d : dice_range) {
		auto source = include_source_dice ? d.generating_function() : typename Dice::poly{};
		for (auto &&[mn, coef] : d.generating_function().terms())
			symbol_set->add_symbol({source, mn});
	}

	return symbol_set;
}

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
					ASSERT(num);
					mn.add_exponent(idx, num.value());
				}
			}
			else {
				auto num = get_num(part);
				ASSERT(num && (S::N == 1 || *num == 0), "Invalid input", text, word, part);
				if (*num) mn.add_exponent(0, num.value());
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

template<spaceless::SymbolSet S, std::signed_integral T>
struct ankerl::unordered_dense::hash<spaceless::dice<S, T>> {
	using is_avalanching = void;

	auto operator()(const spaceless::dice<S, T> &d) const noexcept -> uint64_t {
		return spaceless::get_hash(d.generating_function());
	}
};

template<size_t N, std::signed_integral T>
struct ankerl::unordered_dense::hash<spaceless::combination_dice_face<N, T>> {
	using is_avalanching = void;

	auto operator()(const spaceless::combination_dice_face<N, T> &f) const noexcept -> uint64_t {
		return spaceless::hash_combine(spaceless::get_hash(f.source), spaceless::get_hash(f.face));
	}
};
