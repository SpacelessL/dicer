#pragma once

#include <ranges>
#include <cassert>
#include <random>
#include <utility>
#include <iostream>
#include <concepts>
#include <type_traits>

#include "external/ankerl/unordered_dense.h"

namespace spaceless {

// TODO : dynamic sized indexer/Monomial
// Generic reroll logic

static constexpr uint8_t DYNAMIC = 0;

#define OVERLOADED(...) [&](auto &&...args) -> decltype(auto) { return (__VA_ARGS__)(std::forward<decltype(args)>(args)...); }

inline constexpr double eps = 1e-15;

template <typename R, typename V>
concept range_of = std::ranges::range<R> && std::convertible_to<std::ranges::range_value_t<R>, V>;

// get random number in [0.0, 1.0]
[[nodiscard]]
inline double get_random() {
	thread_local std::random_device rd;
	thread_local std::mt19937 gen(rd());
	thread_local std::uniform_real_distribution<> dis(0.0, 1.0);
	return dis(gen);
}

template<typename K, typename V>
using unordered_map = ankerl::unordered_dense::map<K, V>;

template<typename K, typename V, typename Cmp = std::less<>>
inline auto get_ordered(const unordered_map<K, V> &m, const Cmp &cmp = Cmp{}) {
	std::vector<std::pair<K, V>> ret(m.begin(), m.end());
	std::sort(ret.begin(), ret.end(), cmp);
	return ret;
}

template<typename... Args>
inline auto make_reference_range(Args&... args) {
	return std::views::iota(0, int(sizeof...(args))) | std::views::transform([arr = std::vector{ &args... }](int i) -> decltype(auto) { return *arr[i]; });
}

template<std::floating_point T = double>
class accumulator {
public:
	accumulator(T init = 0) : sum(init), comp(0) {}
	accumulator &operator += (T x) noexcept { Add(x); return *this; }
	accumulator &operator -= (T x) noexcept { return *this += (-x); }
	operator T() const noexcept { return GetSum(); }
private:
	void Add(T x) noexcept {
		T tmp = sum + x;
		if (std::fabs(sum) >= std::fabs(x))
			comp += (sum - tmp) + x;
		else
			comp += (x - tmp) + sum;
		sum = tmp;
	}
	T GetSum() const noexcept { return sum + comp; }
	T sum, comp;
};

template<typename T>
concept Hashable = requires(T t) { unordered_map<T, uint8_t>{}; };

static uint64_t hash_combine(uint64_t lhs, uint64_t rhs) { return lhs ^ (rhs + 0x517cc1b727220a95 + (lhs << 6) + (lhs >> 2)); }
template<Hashable T>
static uint64_t get_hash(const T &t) { return ankerl::unordered_dense::hash<T>()(t); }

}

template<typename T>
struct ankerl::unordered_dense::hash<std::optional<T>> {
	using is_avalanching = void;

	[[nodiscard]] auto operator()(const std::optional<T> &opt_t) const noexcept -> uint64_t {
		return opt_t ? spaceless::get_hash(*opt_t) : 0;
	}
};

template<typename K, typename V>
struct ankerl::unordered_dense::hash<std::pair<K, V>> {
	using is_avalanching = void;

	[[nodiscard]] auto operator()(const std::pair<K, V> &pair) const noexcept -> uint64_t {
		return spaceless::hash_combine(spaceless::get_hash(pair.first), spaceless::get_hash(pair.second));
	}
};

template<std::ranges::range Range>
struct ankerl::unordered_dense::hash<Range> {
	using is_avalanching = void;

	[[nodiscard]] auto operator()(const Range &range) const noexcept -> uint64_t {
		uint64_t ret = 0;
		for (auto &&x : range)
			ret = spaceless::hash_combine(ret, spaceless::get_hash(x));
		return ret;
	}
};

template<std::totally_ordered K, std::totally_ordered V>
struct ankerl::unordered_dense::hash<spaceless::unordered_map<K, V>> {
	using is_avalanching = void;

	[[nodiscard]] auto operator()(const spaceless::unordered_map<K, V> &map) const noexcept -> uint64_t {
		auto vec = spaceless::get_ordered(map);
		std::sort(vec.begin(), vec.end());
		return spaceless::get_hash(vec);
	}
};
