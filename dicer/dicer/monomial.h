#pragma once

#include "utils.h"

namespace spaceless {

class monomial final {
public:
	monomial() : resource_(&buffer_, buffer_.size()), exponents_(&resource_) {}
	monomial(const monomial &rhs) : monomial() { *this = rhs; }
	monomial &operator = (const monomial &rhs) { exponents_ = rhs.exponents_; return *this; }

	template<std::ranges::input_range Range>
	monomial(const Range &range) : exponents_(std::ranges::begin(range), std::ranges::end(range)) { trim(); }
	monomial(const std::initializer_list<int> &list) : exponents_(std::ranges::begin(list), std::ranges::end(list)) { trim(); }

	void resize(int n) { exponents_.resize(n); }

	int dimension() const { return int(exponents_.size()); }

	template<typename Self> auto &&operator[] (this Self &&self, int idx) { return std::forward<Self>(self).exponents_[idx]; }
	template<typename Self> auto data(this Self &&self) { return std::forward<Self>(self).exponents_.data(); }
	template<typename Self> auto &&begin(this Self &&self) { return std::forward<Self>(self).exponents_.begin(); }
	template<typename Self> auto &&end(this Self &&self) { return std::forward<Self>(self).exponents_.end(); }

	int get_exponent(int idx) const { return idx < dimension() ? exponents_[idx] : 0; }
	void set_exponent(int idx, int exponent) { if (!exponent) return; if (idx >= dimension()) resize(idx + 1); exponents_[idx] = exponent; }

	monomial &invert() { for (int &x : *this) x = -x; return *this; }
	monomial inverse() const { auto ret = *this; return ret.invert(); }

	monomial &trim() {
		while (!exponents_.empty() && !exponents_.back())
			exponents_.pop_back();
		return *this;
	}

	monomial &operator *= (const monomial &rhs) {
		if (rhs.dimension() > dimension())
			exponents_.resize(rhs.dimension());
		for (int i = 0; i < std::min(dimension(), rhs.dimension()); i++)
			exponents_[i] += rhs.exponents_[i];
		return trim();
	}
	monomial &operator /= (const monomial &rhs) { return *this *= rhs.inverse(); }
	monomial operator * (const monomial &rhs) const { auto ret = *this; return ret *= rhs; }
	monomial operator / (const monomial &rhs) const { auto ret = *this; return ret /= rhs; }
	auto operator <=> (const monomial &rhs) const { return exponents_ <=> rhs.exponents_; }
	bool operator == (const monomial &rhs) const { return exponents_ == rhs.exponents_; }

private:
	std::array<std::byte, sizeof(int) * 16> buffer_{};
	std::pmr::monotonic_buffer_resource resource_;
	std::pmr::vector<int> exponents_;
};

}

template<>
struct ankerl::unordered_dense::hash<spaceless::monomial> {
	using is_avalanching = void;

	auto operator()(const spaceless::monomial &m) const noexcept -> uint64_t {
		return detail::wyhash::hash(m.data(), m.dimension() * sizeof(int));
	}
};
