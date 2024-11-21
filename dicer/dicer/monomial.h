#pragma once

#include "utils.h"

namespace spaceless {

template<size_t N, std::signed_integral T = int8_t>
class alignas(N ? std::bit_ceil(N * sizeof(T)) : 8) monomial final {
public:
	static constexpr size_t num_symbol = N;
	using underlying_type = T;

	constexpr monomial() noexcept = default;
	constexpr monomial(const monomial &) noexcept(N) = default;
	constexpr monomial(monomial &&) noexcept = default;
	constexpr monomial &operator = (const monomial &) noexcept(N) = default;
	constexpr monomial &operator = (monomial &&) noexcept = default;

	constexpr operator monomial<DYNAMIC, T>() const requires(N) { return monomial<DYNAMIC, T>{ *this }; }

	template<std::ranges::sized_range Range>
	explicit monomial(const Range &range) noexcept(N) {
		if constexpr (N)
			std::ranges::copy_n(range, std::min(std::ranges::size(range), N), exponents_.begin());
		else
			exponents_.insert_range(exponents_.end(), range);
	}
	template<typename R>
	explicit constexpr monomial(const std::initializer_list<R> &list) noexcept(N) : monomial(std::views::all(list)) {}

	constexpr monomial operator - () const noexcept(N) { return inverse(); }
	constexpr monomial operator + () const noexcept(N) { return *this; }
	// one
	constexpr static monomial identity() noexcept { return {}; }

	constexpr void resize(int n) noexcept(N) {
		if constexpr (N) {
			if (n < N) [[unlikely]]
				std::ranges::fill(exponents_.begin() + n, exponents_.end(), 0);
		}
		else
			exponents_.resize(n);
	}
	constexpr void resize_if_needed(int n) noexcept(N) {
		if constexpr (!N) {
			if (exponents_.size() < n)
				exponents_.resize(n);
		}
	}

	constexpr int dimension() const noexcept { return exponents_.size(); }

	constexpr auto &&operator[] (this auto &&self, int idx) noexcept { return std::forward<decltype(self)>(self).exponents_[idx]; }
	constexpr auto data(this auto &&self) noexcept { return std::forward<decltype(self)>(self).exponents_.data(); }
	constexpr auto begin(this auto &&self) noexcept { return std::forward<decltype(self)>(self).exponents_.begin(); }
	constexpr auto end(this auto &&self) noexcept { return std::forward<decltype(self)>(self).exponents_.end(); }

	constexpr T get_exponent(int idx) const noexcept { return idx < dimension() ? exponents_[idx] : 0; }
	constexpr void set_exponent(int idx, T exponent) noexcept(N) {
		resize_if_needed(idx + 1);
		exponents_[idx] = exponent;
		if (!exponents_[idx])
			trim();
	}
	constexpr void add_exponent(int idx, T diff) noexcept(N) {
		resize_if_needed(idx + 1);
		exponents_[idx] += diff;
		if (!exponents_[idx])
			trim();
	}

	constexpr monomial &invert() noexcept { for (T &x : exponents_) x = -x; return *this; }
	[[nodiscard]]
	constexpr monomial inverse() const noexcept(N) { auto ret = *this; return ret.invert(); }

	constexpr monomial &trim() noexcept {
		if constexpr (!N)
			while (!exponents_.empty() && !exponents_.back())
				exponents_.pop_back();
		return *this;
	}

	constexpr monomial &operator *= (const monomial &rhs) noexcept(N) {
		int n = N;
		if constexpr (!N) {
			if (rhs.dimension() > dimension())
				resize(rhs.dimension());
			n = std::min(dimension(), rhs.dimension());
		}
		for (int i = 0; i < n; i++)
			exponents_[i] += rhs.exponents_[i];
		if constexpr (!N)
			return trim();
		return *this;
	}
	constexpr monomial &operator /= (const monomial &rhs) noexcept(N) { return *this *= rhs.inverse(); }
	constexpr monomial operator * (const monomial &rhs) const noexcept(N) { auto ret = *this; return ret *= rhs; }
	constexpr monomial operator / (const monomial &rhs) const noexcept(N) { auto ret = *this; return ret /= rhs; }
	constexpr auto operator <=> (const monomial &rhs) const noexcept { return exponents_ <=> rhs.exponents_; }
	constexpr bool operator == (const monomial &rhs) const noexcept { return exponents_ == rhs.exponents_; }

private:
	template<size_t, std::signed_integral>
	friend class polynomial;

	constexpr void set_exponent_unsafe(int idx, int exponent) noexcept {
		exponents_[idx] = exponent;
	}

	std::conditional_t<N != DYNAMIC, std::array<T, N>, std::vector<T>> exponents_{};
};

template<typename T>
struct is_monomial : std::false_type {};
template<size_t N, std::signed_integral T>
struct is_monomial<monomial<N, T>> : std::true_type {};

template<typename T>
concept MonomialType = is_monomial<T>::value;

template<typename F, typename M>
concept MonomialTransformer = MonomialType<M> && requires(F f, M m) {
	{ f(m) } -> MonomialType;
};

}

template<uint8_t N, std::signed_integral T>
struct ankerl::unordered_dense::hash<spaceless::monomial<N, T>> {
	using is_avalanching = void;

	auto operator()(const spaceless::monomial<N, T> &m) const noexcept -> uint64_t {
		return detail::wyhash::hash(m.data(), m.dimension() * sizeof(int));
	}
};
