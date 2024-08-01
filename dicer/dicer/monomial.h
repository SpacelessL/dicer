#pragma once

#include "utils.h"

namespace spaceless {

class monomial final {
public:
	monomial() noexcept = default;
	monomial(const monomial &) noexcept = default;
	monomial &operator = (const monomial &) noexcept = default;

	template<std::ranges::sized_range Range>
	monomial(const Range &range) noexcept : dim_(int(std::ranges::size(range))) { std::ranges::copy(range, exponents_.begin()); trim(); }
	monomial(const std::initializer_list<int> &list) noexcept : monomial(std::views::all(list)) {}

	void resize(int n) noexcept {
		if (n < dim_) [[unlikely]] std::memset(exponents_.data() + n, 0, (dim_ - n) * sizeof(int));
		dim_ = n;
	}

	int dimension() const noexcept { return dim_; }

	// these 4 functions are unsafe, use with your own risk! trim after you finishing editing
	template<typename Self> auto &&operator[] (this Self &&self, int idx) noexcept { return std::forward<Self>(self).exponents_[idx]; }
	template<typename Self> auto data(this Self &&self) noexcept { return std::forward<Self>(self).exponents_.data(); }
	template<typename Self> auto &&begin(this Self &&self) noexcept { return std::forward<Self>(self).exponents_.begin(); }
	template<typename Self> auto &&end(this Self &&self) noexcept { return std::forward<Self>(self).exponents_.end(); }

	int get_exponent(int idx) const noexcept { return exponents_[idx]; }
	void set_exponent_unsafe(int idx, int exponent) {
		exponents_[idx] = exponent;
		dim_ = std::max(idx + 1, dim_);
	}
	void set_exponent(int idx, int exponent) noexcept {
		if (idx >= dimension()) {
			if (!exponent) return;
			dim_ = idx + 1;
			exponents_[idx] = exponent;
			return;
		}
		exponents_[idx] = exponent;
		trim();
	}

	monomial &invert() noexcept { for (int &x : *this) x = -x; return *this; }
	monomial inverse() const noexcept { auto ret = *this; return ret.invert(); }

	monomial &trim() noexcept {
		while (dim_ && !exponents_[dim_ - 1])
			--dim_;
		return *this;
	}

	monomial &operator *= (const monomial &rhs) noexcept  {
		int d = std::min(dimension(), rhs.dimension());
		if (rhs.dimension() > dimension())
			dim_ = rhs.dimension();
		for (int i = 0; i < d; i++)
			exponents_[i] += rhs.exponents_[i];
		return trim();
	}
	monomial &operator /= (const monomial &rhs) noexcept { return *this *= rhs.inverse(); }
	monomial operator * (const monomial &rhs) const noexcept { auto ret = *this; return ret *= rhs; }
	monomial operator / (const monomial &rhs) const noexcept { auto ret = *this; return ret /= rhs; }
	auto operator <=> (const monomial &rhs) const noexcept = default;
	bool operator == (const monomial &rhs) const noexcept = default;

private:
	int dim_ = 0;
	std::array<int, 16> exponents_{};
};

}

template<>
struct ankerl::unordered_dense::hash<spaceless::monomial> {
	using is_avalanching = void;

	auto operator()(const spaceless::monomial &m) const noexcept -> uint64_t {
		return detail::wyhash::hash(m.data(), m.dimension() * sizeof(int));
	}
};
