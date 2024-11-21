#pragma once

#include "polynomial.h"

#include <concepts>
#include <ranges>
#include <type_traits>

namespace spaceless {
// discrete distribution
template<size_t N, std::signed_integral T = int8_t>
class distribution final {
public:
	struct statistics {
		int min = std::numeric_limits<int>::max(), max = std::numeric_limits<int>::min();
		double E = 0, V = 0, D = 0, skewness = 0;
	};

	using mono = monomial<N, T>;
	using poly = polynomial<N, T>;

	distribution() = default;
	distribution(const distribution &) = default;
	distribution(distribution &&) noexcept = default;
	distribution &operator = (const distribution &) = default;
	distribution &operator = (distribution &&) noexcept = default;

	explicit distribution(poly pl) : pl(pl.terms().size() ? std::move(pl) : poly::one()) { normalize(); }

	auto &&generating_function() const { return pl; }

	static distribution zero() {
		distribution ret;
		ret.pl = poly::identity();
		return ret;
	}

	distribution &clamp(double eps = ::spaceless::eps) {
		pl.trim(eps);
		normalize();
		return *this;
	}
	[[nodiscard]]
	distribution clamped(double eps = ::spaceless::eps) const {
		distribution ret = *this;
		ret.clamp(eps);
		return ret;
	}

	void normalize() {
		accumulator<> sum;
		for (auto &coef : pl.terms_ | std::views::values)
			sum += coef;
		for (auto &coef : pl.terms_ | std::views::values)
			coef /= sum;
	}

	// merge range of <distribution, Prob> into a new distribution
	template<std::ranges::input_range Range>
	static distribution merge_distribution(const Range &distribution_with_prob) {
		unordered_map<mono, accumulator<>> map;
		for (auto &&[distribution, prob] : distribution_with_prob)
			for (auto &&[mn, coef] : distribution.pl.terms)
				map[mn] += coef * prob;
		distribution ret;
		ret.pl.terms_.reserve(map.size());
		for (auto &&[mn, prob] : map)
			ret.pl.terms_.emplace(mn, prob);
		ret.normalize();
		return ret;
	}

	distribution operator + (const distribution &rhs) const & {
		distribution ret(pl * rhs.pl);
		ret.normalize();
		return ret;
	}
	distribution operator + (distribution &&rhs) const & { return std::move(rhs += *this); }
	distribution operator + (const distribution &rhs) && { return std::move(*this += rhs); }
	distribution &operator += (const distribution &rhs) {
		pl *= rhs.pl;
		normalize();
		return *this;
	}

	distribution operator * (int times) const {
		distribution ret(pl.power(times));
		ret.normalize();
		return ret;
	}
	friend distribution operator * (int times, const distribution &distribution) { return distribution * times; }

	template<typename F>
	auto transformed(F &&func) const requires MonomialTransformer<F, mono> {
		using new_mono = std::invoke_result<F, mono>;
		unordered_map<new_mono, accumulator<>> terms;
		for (auto &&[mono, coef] : pl.terms())
			terms[func(mono)] += coef;
		return distribution<new_mono::num_symbol, typename new_mono::underlying_type>(std::move(terms));
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
	poly pl;
};

}
