#pragma once

#include <ranges>
#include <cassert>
#include <random>
#include <utility>
#include <iostream>
#include <concepts>
#include <type_traits>

#include "hash.h"

namespace spaceless {

// TODO : Generic reroll logic

#define OVERLOADED(...) [&](auto &&...args) -> decltype(auto) { return (__VA_ARGS__)(std::forward<decltype(args)>(args)...); }

inline constexpr double eps = 1e-15;

template <typename R, typename V>
concept range_of = std::ranges::range<R> && std::convertible_to<std::ranges::range_value_t<R>, V>;

// get random number in [0.0, 1.0]
template<std::floating_point T = double>
[[nodiscard]] inline double get_random() {
	thread_local std::random_device rd;
	thread_local std::mt19937 gen(rd());
	thread_local std::uniform_real_distribution<T> dis(0, 1);
	return dis(gen);
}

template<typename... Args>
inline auto make_reference_range(Args&... args) {
	return std::views::iota(0, int(sizeof...(args))) | std::views::transform([arr = std::vector{ &args... }](int i) -> decltype(auto) { return *arr[i]; });
}

template<std::floating_point T = double>
class accumulator final {
public:
	accumulator(T init = 0) noexcept : sum_(init), comp_(0) {}
	accumulator &operator += (T x) noexcept { add(x); return *this; }
	accumulator &operator -= (T x) noexcept { return *this += -x; }
	operator T() const noexcept { return sum_ + comp_; }
private:
	void add(T x) noexcept {
		T tmp = sum_ + x;
		if (std::fabs(sum_) >= std::fabs(x))
			comp_ += (sum_ - tmp) + x;
		else
			comp_ += (x - tmp) + sum_;
		sum_ = tmp;
	}
	T sum_, comp_;
};

}
