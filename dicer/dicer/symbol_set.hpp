#pragma once

#include "hash.hpp"
#include "utils.hpp"
#include "sus/misc.hpp"

namespace spaceless {

template<typename T>
concept SymbolSet = requires { typename T::symbol; { T::N } -> std::convertible_to<size_t>; }
	&& requires(const T t, int8_t index, typename T::symbol x) {
		{ t.get_index(x) } -> std::integral;
		{ t.from_index(index) } -> std::convertible_to<typename T::symbol>;
		{ t.contains(x) } -> std::convertible_to<bool>;
		{ t.size() } -> std::integral;
	} && (T::N != DYNAMIC
		|| requires(T t, typename T::symbol x) {
			{ t.add_symbol(x) } -> std::integral;
		}
	);

template<auto ...Symbols>
requires (std::is_same_v<head_of_t<decltype(Symbols)...>, decltype(Symbols)> && ...)
class static_symbol_set final {
public:
	using symbol = head_of_t<decltype(Symbols)...>;
	static constexpr size_t N = sizeof...(Symbols);

	constexpr int8_t get_index(const symbol &s) const {
		int8_t i = 0, ret = -1;
		((s == Symbols ? (ret = i, true) : (++i, false)) || ...);
		return ret;
	}

	constexpr symbol from_index(int8_t index) const {
		int8_t i = 0;
		symbol ret{};
		((i++ == index ? (ret = Symbols) : false) || ...);
		return ret;
	}

	constexpr bool contains(const symbol &s) const {
		return ((s == Symbols) || ...);
	}

	constexpr int size() const { return N; }
};

template<Hashable SymbolType>
class dynamic_symbol_set final {
public:
	using symbol = SymbolType;
	static constexpr size_t N = DYNAMIC;

	dynamic_symbol_set() = default;
	dynamic_symbol_set(const dynamic_symbol_set &) = default;
	dynamic_symbol_set(dynamic_symbol_set &&) noexcept = default;
	dynamic_symbol_set &operator = (const dynamic_symbol_set &) = default;
	dynamic_symbol_set &operator = (dynamic_symbol_set &&) noexcept = default;
	~dynamic_symbol_set() noexcept = default;

	template<std::ranges::input_range Range>
	dynamic_symbol_set(Range &&range) { for (auto &&s : range) add_symbol(s); }
	dynamic_symbol_set(std::initializer_list<symbol> init) : dynamic_symbol_set(std::views::all(init)) {}

	void reserve(size_t n) {
		to_index.reserve(n);
		symbols.reserve(n);
	}

	int8_t add_symbol(const symbol &symbol) {
		if (contains(symbol))
			return get_index(symbol);
		int8_t ret = int8_t(to_index.size());
		to_index.emplace(symbol, ret);
		symbols.emplace_back(symbol);
		return ret;
	}

	int8_t get_index(const symbol &symbol) const { if (auto it = to_index.find(symbol); it != to_index.end()) return it->second; return -1; }

	const symbol &from_index(int8_t index) const { return symbols[index]; }

	bool contains(const symbol &symbol) const { return to_index.contains(symbol); }

	int size() const { return symbols.size(); }

private:
	unordered_map<symbol, int8_t> to_index;
	std::vector<symbol> symbols;
};

namespace detail {
template<string_literal Symbols, std::size_t... I>
auto make_simple_static_symbol_set(std::index_sequence<I...>) {
	static_assert((((Symbols[I] < '0' || Symbols[I] > '9') && Symbols[I] != '-' && Symbols[I] != ',') && ...));
	return std::make_shared<static_symbol_set<Symbols[I]...>>();
}
}

template<string_literal Symbols>
auto make_simple_static_symbol_set() {
	return detail::make_simple_static_symbol_set<Symbols>(std::make_index_sequence<Symbols.length()>{});
}

}
