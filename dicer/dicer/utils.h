#pragma once

#include <ranges>
#include <cassert>
#include <random>
#include <utility>
#include <iostream>
#include <concepts>
#include <type_traits>
#include <stacktrace>
#include <source_location>

#include "hash.h"

namespace spaceless {

static constexpr size_t DYNAMIC = 0;

namespace detail {
inline auto check_ensure_value(auto &&value, const char *reason = "Unspecified Reason") { return std::make_pair(bool(value), reason); }
[[noreturn]] inline void report_error(const char *assertion, const char *reason, std::stacktrace trace = std::stacktrace::current(), std::source_location loc = std::source_location::current()) {
	std::cerr << std::format(
		"{}: ENSURE({}) failed\n"
		"At {}({}:{}) `{}`\n"
		"Stacktrace:\n{}\n"
		, reason, assertion
		, loc.file_name(), loc.line(), loc.column(), loc.function_name()
		, trace);
	throw std::runtime_error(reason);
}
}

#define ENSURE(...) if (auto [value, reason] = detail::check_ensure_value(__VA_ARGS__); !value) [[unlikely]] detail::report_error(#__VA_ARGS__, reason)

#define OVERLOADED(...) [&](auto &&...args) -> decltype(auto) { return (__VA_ARGS__)(std::forward<decltype(args)>(args)...); }

inline constexpr double eps = 1e-15;

template <typename R, typename V>
concept range_of = std::ranges::range<R> && std::convertible_to<std::ranges::range_value_t<R>, V>;

template<size_t N>
struct string_literal {
	std::array<char, N> value;

	constexpr string_literal(const char (&str)[N + 1]) { std::ranges::copy_n(str, N, value.begin()); }
	constexpr operator std::string_view() const { return { value.data(), value.size() }; }
	constexpr size_t length() const { return N; }
	constexpr char operator[] (size_t idx) const { return value[idx]; }
	constexpr auto operator <=> (const string_literal &) const noexcept = default;
	constexpr bool operator == (const string_literal &) const noexcept = default;
};

template<size_t N>
string_literal(const char(&)[N]) -> string_literal<N - 1>;

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

struct timer {
	timer(std::string_view name, bool log = true) : start(std::chrono::steady_clock::now()), name(name), log(log) {
		if (log)
			std::cout << std::format("[{}] started! \n", name);
	}
	std::chrono::duration<double> elapsed() const {
		return std::chrono::steady_clock::now() - start;
	}
	~timer() {
		if (log)
			std::cout << std::format("[{}] ended, total time : {}s\n", name, elapsed().count());
	}

	std::chrono::time_point<std::chrono::steady_clock> start;
	std::string_view name;
	bool log;
};

template<std::floating_point T = double>
class accumulator final {
public:
	accumulator(T init = 0) noexcept : sum_(init), comp_(0) {}
	template<std::ranges::input_range Range>
	accumulator(Range &&range) noexcept : accumulator() { for (auto &&x : range) add(x); }
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

inline auto squared(const auto &t) { return t * t; }
inline auto cubed(const auto &t) { return t * t * t; }

}
