#pragma once

#include "utils.h"

#include "external/ankerl/unordered_dense.h"

namespace spaceless {

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
