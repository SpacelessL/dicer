#pragma once

#include "monomial.h"
#include "utils.h"
#include "fft.h"

#include <iostream>
#include <limits>
#include <map>

namespace spaceless {

struct polynomial final {
	using mono = monomial;

	polynomial() = default;
	polynomial(const polynomial &) = default;
	polynomial(polynomial &&) noexcept = default;
	polynomial &operator = (const polynomial &) = default;
	polynomial &operator = (polynomial &&) noexcept = default;

	polynomial(const mono &mono, double coef = 1) { terms.emplace(mono, coef); }
	explicit polynomial(double x) : polynomial({}, x) {}

	static polynomial identity() { return { {}, 1 }; }

	polynomial operator - () const & { return -polynomial(*this); }
	polynomial operator + () const & { return *this; }
	polynomial operator - () && { for (auto &coef : terms | std::views::values) coef = -coef; return std::move(*this); }
	polynomial operator + () && { return std::move(*this); }

	polynomial operator + (const polynomial &rhs) const & { return polynomial(*this) + rhs; }
	polynomial operator - (const polynomial &rhs) const & { return polynomial(*this) - rhs; }
	polynomial operator + (const polynomial &rhs) && { *this += rhs; return std::move(*this); }
	polynomial operator - (const polynomial &rhs) && { *this -= rhs; return std::move(*this); }
	polynomial operator * (const polynomial &rhs) const { return fft_multiply(rhs); }
	polynomial operator / (const polynomial &rhs) const { return divide(rhs, nullptr); }
	polynomial operator % (const polynomial &rhs) const { polynomial ret; divide(rhs, &ret); return ret; }

	polynomial &operator += (const polynomial &rhs) { for (auto &&[mono, coef] : rhs.terms) terms[mono] += coef; return trim(); }
	polynomial &operator -= (const polynomial &rhs) { for (auto &&[mono, coef] : rhs.terms) terms[mono] -= coef; return trim(); }
	polynomial &operator *= (const polynomial &rhs) { return *this = *this * rhs; }
	polynomial &operator /= (const polynomial &rhs) { return *this = *this / rhs; }
	polynomial &operator %= (const polynomial &rhs) { return *this = *this % rhs; }

#ifdef __RESHARPER__
	friend polynomial operator + (polynomial poly, double c) { return poly; } friend polynomial operator + (double c, polynomial poly) { return poly; }
	friend polynomial operator - (polynomial poly, double c) { return poly; } friend polynomial operator - (double c, polynomial poly) { return poly; }
	friend polynomial operator * (polynomial poly, double c) { return poly; } friend polynomial operator * (double c, polynomial poly) { return poly; }
	friend polynomial operator / (polynomial poly, double c) { return poly; } friend polynomial operator / (double c, polynomial poly) { return poly; }
#else
#define POLYNOMIAL_OPERATOR_WITH_DOUBLE(op) \
	friend polynomial operator op (double c, const polynomial &poly) { return poly op c; } \
	friend polynomial operator op (double c, polynomial &&poly) { return std::move(poly) op c; } \
	polynomial operator op (double c) const & { return polynomial(*this) op c; } \
	polynomial operator op (double c) && { *this op= c; return std::move(*this); }
	POLYNOMIAL_OPERATOR_WITH_DOUBLE(+)
	POLYNOMIAL_OPERATOR_WITH_DOUBLE(-)
	POLYNOMIAL_OPERATOR_WITH_DOUBLE(*)
	POLYNOMIAL_OPERATOR_WITH_DOUBLE(/)
#undef POLYNOMIAL_OPERATOR_WITH_DOUBLE
#endif

	polynomial &operator += (double c) { terms[{}] += c; return *this; }
	polynomial &operator -= (double c) { terms[{}] -= c; return *this; }
	polynomial &operator *= (double c) { for (auto &coef : terms | std::views::values) coef *= c; return *this; }
	polynomial &operator /= (double c) { for (auto &coef : terms | std::views::values) coef /= c; return *this; }

	bool operator == (const polynomial &rhs) const = default;
	bool operator != (const polynomial &rhs) const = default;

	[[nodiscard]] polynomial power(unsigned int exp) const {
		assert(exp >= 0);
		if (exp == 0) return identity();
		if (exp == 1) return *this;
		unsigned int m = exp / 2;
		polynomial rt = power(m);
		return exp % 2 ? rt * rt * *this : rt * rt;
	}

	polynomial divide(const polynomial &divisor, polynomial *remainder) const {
		// *this / d = q ... r
		std::map<mono, double, std::greater<>> r(terms.begin(), terms.end()), d(divisor.terms.begin(), divisor.terms.end()), q;
		while (!r.empty() && !d.empty() && r.begin()->first >= d.begin()->first) {
			auto s_mono = r.begin()->first / d.begin()->first;
			auto s_coef = r.begin()->second / d.begin()->second;
			q[s_mono] += s_coef;
			for (auto &&[mono, coef] : d)
				r[mono * s_mono] -= coef * s_coef;
			std::erase_if(r, [](auto &&pair) { return std::abs(pair.second) < eps; });
		}
		if (remainder)
			remainder->terms = unordered_map<mono, double>(r.begin(), r.end());
		polynomial ret;
		ret.terms = unordered_map<mono, double>(q.begin(), q.end());
		return ret;
	}

	int dimension() const {
		int m = 0;
		for (const auto &mono : terms | std::views::keys)
			m = std::max(m, mono.dimension());
		return m;
	}

	polynomial &trim(double eps = ::spaceless::eps) { std::erase_if(terms, [&](auto &&p) { return std::abs(p.second) <= eps; }); return *this; }
	[[nodiscard]] polynomial trimmed() const { polynomial ret = *this; ret.trim(); return ret; }

	unordered_map<mono, double> terms; // mono -> coef

private:

	polynomial naive_multiply(const polynomial &rhs) const {
		polynomial poly;
		for (auto &&[lmono, lcoef] : terms)
			for (auto &&[rmono, rcoef] : rhs.terms)
				poly.terms[lmono * rmono] += lcoef * rcoef;
		return poly;
	}
	
	polynomial fft_multiply(const polynomial &rhs) const {
		int D = std::max(dimension(), rhs.dimension());
		auto get_value = [D](const polynomial &poly, int init, auto &&select) {
			std::vector ret(D, init);
			for (const auto &mono : poly.terms | std::views::keys)
				for (int i = 0; i < D; i++)
					ret[i] = select(ret[i], mono.get_exponent(i));
			return ret;
		};
		auto get_min = [&](const polynomial &poly) { return get_value(poly, std::numeric_limits<int>::max(), OVERLOADED(std::min)); };
		auto get_max = [&](const polynomial &poly) { return get_value(poly, std::numeric_limits<int>::min(), OVERLOADED(std::max)); };

		auto minl = get_min(*this), maxl = get_max(*this), minr = get_min(rhs), maxr = get_max(rhs);

		std::vector<int64_t> shape(D), shape_stride(D);
		size_t ndata = 1;
		for (int i = shape.size() - 1; i >= 0; i--) {
			shape[i] = maxl[i] + maxr[i] - minl[i] - minr[i] + 1;
			shape_stride[i] = ndata;
			ndata *= shape[i];
		}

#ifndef USE_ONEMKL
		// simply because PocketFFT isn't fast enough
		if (terms.size() * rhs.terms.size() <= ndata * std::log2(ndata) / 20)
			return naive_multiply(rhs);
#else
		if (terms.size() == 1 || rhs.terms.size() == 1)
			return naive_multiply(rhs);
#endif

		auto get_flat_index = [&](const mono &mono, const std::vector<int> &min) {
			size_t ret = 0;
			for (int i = 0; i < D; i++)
				ret += (mono.get_exponent(i) - min[i]) * shape_stride[i];
			assert(0 <= ret && ret < ndata);
			return ret;
		};
		auto from_flat_index = [&](size_t flat_index) {
			mono ret;
			for (int i = 0; i < D; i++) {
				ret.set_exponent_unsafe(i, flat_index / shape_stride[i] + minl[i] + minr[i]);
				flat_index %= shape_stride[i];
			}
			return ret.trim();
		};
		auto get_input = [&](const polynomial &poly, const std::vector<int> &min) {
			std::vector<std::complex<double>> ret(ndata);
			for (auto &&[mono, coef] : poly.terms)
				ret[get_flat_index(mono, min)] += coef;
			return ret;
		};

		auto outputl = detail::fft(get_input(*this, minl), shape), outputr = detail::fft(get_input(rhs, minr), shape);
		auto mul = detail::dot(outputl, outputr);
		auto result = detail::ifft(mul, shape);

		polynomial ret;
		auto keep = [](auto &&c) { return std::abs(c.real()) > eps; };
		ret.terms.reserve(std::ranges::count_if(result, keep));
		for (size_t flat_index = 0; flat_index < ndata; flat_index++)
			if (auto v = result[flat_index]; keep(v))
				ret.terms.emplace(from_flat_index(flat_index), v.real());
		return ret;
	}
};

static std::ostream &operator << (std::ostream &os, const polynomial &x) {
	auto to_var = [](int i) { return i < 3 ? char('x' + i) : char('a' + i - 3); };
	bool first = true;
	for (auto &&[mono, coef] : get_ordered(x.terms, std::greater{})) {
		if (coef >= 0 && !first)
			os << '+';
		if (std::abs(coef - 1) > eps)
			os << coef;
		bool is_C = true;
		for (int i = 0; i < mono.dimension(); i++)
			if (int exp = mono[i]) {
				os << to_var(i);
				if (exp != 1)
					os << '^' << exp;
				is_C = false;
			}
		if (is_C && std::abs(coef - 1) <= eps)
			os << coef;
		first = false;
	}
	if (x.terms.empty())
		os << 0;
	return os;
}

}

template<>
struct ankerl::unordered_dense::hash<spaceless::polynomial> {
	using is_avalanching = void;

	auto operator()(const spaceless::polynomial &poly) const noexcept -> uint64_t {
		return get_hash(poly.terms);
	}
};
