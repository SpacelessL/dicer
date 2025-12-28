#pragma once

#include <ranges>
#include <algorithm>

#include "external/ankerl/unordered_dense.h"

namespace spaceless {

template<typename K, typename V>
using unordered_map = ankerl::unordered_dense::map<K, V>;

template<typename T>
concept Hashable = requires(T t) { unordered_map<T, uint8_t>{}; };

static uint64_t hash_combine(uint64_t lhs, uint64_t rhs) { return lhs ^ (rhs + 0x517cc1b727220a95 + (lhs << 6) + (lhs >> 2)); }
template<Hashable T>
static uint64_t get_hash(const T &t) { return ankerl::unordered_dense::hash<T>()(t); }

template<typename K, typename V, typename Cmp = std::less<>>
inline auto get_ordered(const unordered_map<K, V> &m, const Cmp &cmp = Cmp{}) {
	std::vector<std::pair<K, V>> ret(m.begin(), m.end());
	std::ranges::sort(ret, cmp);
	return ret;
}

}

template<typename T>
struct ankerl::unordered_dense::hash<std::optional<T>> {
	using is_avalanching = void;

	[[nodiscard]] auto operator()(const std::optional<T> &opt_t) const noexcept -> uint64_t {
		return opt_t ? spaceless::get_hash(*opt_t) : spaceless::get_hash(0);
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
		std::vector<uint64_t> hashes(map.size());
		for (auto &&[index, p] : std::views::enumerate(map))
			hashes[index] = spaceless::get_hash(p);
		std::ranges::sort(hashes);
		return spaceless::get_hash(hashes);
	}
};
