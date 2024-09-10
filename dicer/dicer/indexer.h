#pragma once

#include "utils.h"

#include "external/ankerl/unordered_dense.h"

namespace spaceless {

template<auto FirstSymbol, auto ...RestSymbols>
requires (std::is_same_v<decltype(FirstSymbol), decltype(RestSymbols)> && ...)
class static_symbol_set final {
public:
	using symbol = decltype(FirstSymbol);

	template<symbol SymbolValue>
	uint8_t get_index() const noexcept {
		if constexpr (SymbolValue == FirstSymbol) return 0;
		return 1 + static_symbol_set<RestSymbols...>{}.template get_index<SymbolValue>();
	}
	template<uint8_t Index>
	symbol from_index() const noexcept {
		if constexpr (Index == 0) return FirstSymbol;
		return static_symbol_set<RestSymbols...>{}.template from_index<Index - 1>();
	}
	template<symbol SymbolValue>
	bool contains() const noexcept {
		if constexpr (FirstSymbol == SymbolValue) return true;
		return static_symbol_set<RestSymbols...>{}.template contains<SymbolValue>();
	}

	uint8_t get_index(const symbol &symbol) const {
		if (symbol == FirstSymbol) return 0;
		return 1 + static_symbol_set<RestSymbols...>{}.get_index(symbol);
	}

	symbol from_index(uint8_t index) const {
		if (index == 0) return FirstSymbol;
		return static_symbol_set<RestSymbols...>{}.from_index(index - 1);
	}

	bool contains(const symbol &symbol) const { return symbol == FirstSymbol || ((symbol == RestSymbols) || ...); }
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
	indexer(InputIter begin, InputIter end) { while (begin != end) add_symbol(*(begin++)); }
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
