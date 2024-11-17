#pragma once

#include "utils.h"

namespace spaceless {

template<typename T>
concept SymbolSet = requires { typename T::symbol; { T::dynamic } -> std::convertible_to<bool>; }
	&& requires(const T t, uint8_t index, typename T::symbol x) {
		{ t.get_index(x) } -> std::integral;
		{ t.from_index(index) } -> std::convertible_to<typename T::symbol>;
		{ t.contains(x) } -> std::convertible_to<bool>;
		{ t.size() } -> std::integral;
	} && (!T::dynamic
		|| requires(T t, typename T::symbol x) {
			{ t.add_symbol(x) } -> std::integral;
		}
	);

template<auto FirstSymbol, auto ...RestSymbols>
requires (std::is_same_v<decltype(FirstSymbol), decltype(RestSymbols)> && ...)
class static_symbol_set final {
public:
	using symbol = decltype(FirstSymbol);
	static constexpr bool dynamic = false;
	static_assert(SymbolSet<static_symbol_set>);

	constexpr uint8_t get_index(const symbol &symbol) const {
		if (symbol == FirstSymbol) return 0;
		return 1 + static_symbol_set<RestSymbols...>{}.get_index(symbol);
	}

	constexpr symbol from_index(uint8_t index) const {
		if (index == 0) return FirstSymbol;
		return static_symbol_set<RestSymbols...>{}.from_index(index - 1);
	}

	constexpr bool contains(const symbol &symbol) const {
		return symbol == FirstSymbol || ((symbol == RestSymbols) || ...);
	}

	constexpr int size() const { return 1 + sizeof...(RestSymbols); }
};

template<Hashable SymbolType>
class dynamic_symbol_set final {
public:
	using symbol = SymbolType;
	static constexpr bool dynamic = true;
	static_assert(SymbolSet<dynamic_symbol_set>);

	dynamic_symbol_set() = default;
	dynamic_symbol_set(const dynamic_symbol_set &) = default;
	dynamic_symbol_set(dynamic_symbol_set &&) noexcept = default;
	dynamic_symbol_set &operator = (const dynamic_symbol_set &) = default;
	dynamic_symbol_set &operator = (dynamic_symbol_set &&) noexcept = default;
	~dynamic_symbol_set() noexcept = default;

	template<std::ranges::input_range Range>
	dynamic_symbol_set(Range &&range) { for (auto &&s : range) add_symbol(s); }
	dynamic_symbol_set(std::initializer_list<symbol> init) : dynamic_symbol_set(std::views::all(init)) {}

	uint8_t add_symbol(const symbol &symbol) {
		if (contains(symbol))
			return get_index(symbol);
		uint8_t ret = uint8_t(to_index.size());
		to_index.emplace(symbol, ret);
		symbols.emplace_back(symbol);
		return ret;
	}

	uint8_t get_index(const symbol &symbol) const { return to_index.at(symbol); }

	const symbol &from_index(uint8_t index) const { return symbols[index]; }

	bool contains(const symbol &symbol) const { return to_index.contains(symbol); }

	int size() const { return symbols.size(); }

private:
	unordered_map<symbol, uint8_t> to_index;
	std::vector<symbol> symbols;
};

template<Hashable SymbolType>
class indexer final {
public:
	using symbol = SymbolType;

	indexer() = default;
	indexer(const indexer &) = default;
	indexer(indexer &&) noexcept = default;
	indexer &operator = (const indexer &) = default;
	indexer &operator = (indexer &&) noexcept = default;
	~indexer() = default;

	template<std::input_iterator InputIter>
	indexer(InputIter begin, InputIter end) { while (begin != end) add_symbol(*begin++); }
	indexer(std::initializer_list<symbol> init) : indexer(init.begin(), init.end()) {}
	indexer(const std::vector<symbol> &vec) : indexer(vec.begin(), vec.end()) {}

	uint8_t add_symbol(const symbol &symbol) {
		if (contains(symbol))
			return get_index(symbol);
		uint8_t ret = uint8_t(to_index.size());
		symbols.emplace_back(symbol);
		to_index.emplace(symbol, ret);
		return ret;
	}

	uint8_t get_index(const symbol &symbol) const { return to_index.at(symbol); }

	const symbol &from_index(uint8_t index) const { return symbols[index]; }

	bool contains(const symbol &symbol) const { return to_index.contains(symbol); }

private:
	unordered_map<symbol, uint8_t> to_index;
	std::vector<symbol> symbols;
};

template<Hashable Symbol, Hashable... Args>
auto MakeIndexer(Symbol &&symbol, Args&&... args) {
	indexer<std::remove_cvref_t<Symbol>> ret;
	ret.add_symbol(std::forward<Symbol>(symbol));
	(ret.add_symbol(std::forward<Args>(args)), ...);
	return std::make_shared<decltype(ret)>(ret);
}

template<std::input_iterator Iter>
auto MakeIndexer(Iter begin, Iter end) {
	return std::make_shared<indexer<std::remove_cvref_t<decltype(*begin)>>>(begin, end);
}

}
