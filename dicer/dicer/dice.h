#pragma once

#include "polynomial.h"
#include "indexer.h"

#include <concepts>
#include <functional>
#include <ranges>
#include <algorithm>
#include <type_traits>

namespace spaceless {

template<Hashable Symbol>
using dice_face = unordered_map<Symbol, int>;
template<typename Dice>
using combination_dice_face_symbol = std::pair<dice_face<typename Dice::symbol>, const Dice *>;
template<typename Dice>
using combination_dice_face = dice_face<combination_dice_face_symbol<Dice>>;

template<Hashable Symbol>
int count_symbol(const dice_face<Symbol> &face, const Symbol &s = {}) { if (auto it = face.find(s); it != face.end()) return it->second; return 0; }

template<Hashable>
class dice;

struct number_dice_symbol {
	auto operator <=> (const number_dice_symbol &rhs) const = default;
	operator char() const { return '.'; }
};

template<typename Func, typename SymbolType>
concept reroll_func = requires (Func &&func, const dice_face<SymbolType> &face) {
	{ func(face) } -> std::convertible_to<dice_face<SymbolType>>;
};

template<typename Func, typename SymbolType>
concept explode_func = requires (Func && func, const dice_face<SymbolType> &face) {
	{ func(face) } -> std::convertible_to<std::optional<dice<SymbolType>>>;
};

}

template<>
struct ankerl::unordered_dense::hash<spaceless::number_dice_symbol> {
	using is_avalanching = void;
	uint64_t operator()(const spaceless::number_dice_symbol &) const noexcept { return 0; }
};

namespace spaceless {

namespace detail {
template<typename Map>
void add_to_map(Map *map) {}

template<typename Map, typename Symbol, typename... Args>
void add_to_map(Map *map, const Symbol &first_symbol, int first_count, Args &&... args) {
	map->operator[](first_symbol) += first_count;
	add_to_map(map, std::forward<decltype(args)>(args)...);
}

inline auto make_number_dice_indexer() {
	auto ret = std::make_shared<indexer<number_dice_symbol>>();
	ret->add_symbol(number_dice_symbol{});
	return ret;
}
}

inline auto number_dice_indexer() { static auto ret = detail::make_number_dice_indexer(); return ret; }

template<typename Symbol, typename... Args>
auto make_dice_face(const Symbol &first_symbol, int first_count, Args &&... args) {
	unordered_map<Symbol, int> ret;
	detail::add_to_map(&ret, first_symbol, first_count, std::forward<decltype(args)>(args)...);
	return dice_face<Symbol>(ret.begin(), ret.end());
}

template<typename Symbol>
auto make_dice_face(const Symbol &s) { return make_dice_face(s, 1); }

}

namespace spaceless {

using number_dice = dice<number_dice_symbol>;
using number_dice_face = dice_face<number_dice_symbol>;

struct result_statistics {
	void debug_log() const {
		std::cout << "min : " << min << ", max : " << max << std::endl;
		std::cout << "E : " << E << ", V : " << V << ", D : " << D << std::endl;
	}

	int min = std::numeric_limits<int>::max(), max = std::numeric_limits<int>::min();
	accumulator<> E, E2;
	double V = 0, D = 0;
};

template<Hashable SymbolType>
class dice final {
public:
	using mono = monomial;
	using poly = polynomial;
	using symbol = SymbolType;

	dice() = default;
	dice(const dice &) = default;
	dice(dice &&) noexcept = default;
	dice &operator = (const dice &) = default;
	dice &operator = (dice &&) noexcept = default;

	// range<range<symbol, count>, prob>
	template<std::input_iterator Iter>
	dice(Iter faces_begin, Iter faces_end, std::shared_ptr<indexer<symbol>> input_idx) : idx(std::move(input_idx)) {
		for (auto &&[face, coef] : std::ranges::subrange(faces_begin, faces_end)) {
			std::vector<int> arr;
			int dim = 0;
			for (auto &&[symbol, count] : face)
				dim = std::max(dim, idx->get_index(symbol) + 1);
			arr.resize(dim);
			for (auto &&[symbol, count] : face)
				arr[idx->get_index(symbol)] += count;
			pl.terms[arr] += coef;
		}
	}
	// range<range<symbol, count>, prob>
	template<std::ranges::input_range Range>
	dice(const Range &range, std::shared_ptr<indexer<symbol>> idx) : dice(std::ranges::begin(range), std::ranges::end(range), idx) {}

	auto get_indexer() const { return idx; }

	auto to_number_dice(const symbol &s) const {
		int index = idx->get_index(s);
		unordered_map<int, accumulator<>> values;
		for (auto &&[mn, coef] : pl.terms)
			values[mn.get_exponent(index)] += coef;
		number_dice ret;
		for (auto [num, prob] : values)
			ret.pl.terms[std::array<int, 1>{ num }] = prob;
		ret.idx = number_dice_indexer();
		return ret;
	}

	// returns range<dice_face, Prob>
	auto to_faces_with_prob() const { return pl.terms | std::views::transform([this](auto &&p) { return std::make_pair(to_face(p.first), p.second); }); }

	static dice blank_dice(std::shared_ptr<indexer<symbol>> idx) {
		dice ret;
		ret.idx = std::move(idx);
		ret.pl = poly::identity();
		return ret;
	}

	dice &clamp(double eps = ::spaceless::eps) { pl.trim(eps); return *this; }
	[[nodiscard]]
	dice clamped(double eps = ::spaceless::eps) const { dice ret = *this; ret.clamp(eps); return ret; }

	void normalize() {
		accumulator<> sum;
		for (auto &coef : pl.terms | std::views::values)
			sum += coef;
		for (auto &coef : pl.terms | std::views::values)
			coef /= sum;
	}

	// merge range of pair<dice, Prob> into a new dice
	template<std::input_iterator Iter>
	static dice merge_dice(Iter dice_with_prob_begin, Iter dice_with_prob_end) {
		unordered_map<mono, accumulator<>> map;
		for (auto &&[dice, prob] : std::ranges::subrange(dice_with_prob_begin, dice_with_prob_end))
			for (auto &&[mn, coef] : dice.pl.terms)
				map[mn] += coef * prob;
		dice ret;
		ret.idx = dice_with_prob_begin->first.idx;
		ret.pl.terms.reserve(map.size());
		for (auto &&[mn, prob] : map)
			ret.pl.terms.emplace(mn, prob);
		ret.normalize();
		return ret;
	}
	template<std::ranges::input_range Range>
	static dice merge_dice(const Range &dice_with_prob) { return merge_dice(std::ranges::begin(dice_with_prob), std::ranges::end(dice_with_prob)); }

	dice operator + (const dice &rhs) const { dice ret; ret.idx = idx, ret.pl = pl * rhs.pl; ret.normalize(); return ret; }
	dice &operator += (const dice &rhs) { *this = *this + rhs; return *this; }

	dice operator * (int times) const { dice ret; ret.idx = idx, ret.pl = pl.power(times); ret.normalize(); return ret; }
	friend dice operator * (int times, const dice &dice) { return dice * times; }

	// edit the 0
	dice operator + (int num) const {
		dice ret;
		ret.idx = idx;
		ret.pl.terms.reserve(pl.terms.size());
		for (auto [mn, prob] : pl.terms) {
			mn.resize(1);
			mn[0] += num;
			ret.pl.terms.emplace(mn, prob);
		}
		return ret;
	}
	dice operator - (int num) const { return *this + (-num); }
	dice &operator += (int num) { *this = (*this + num); return *this; }
	dice &operator -= (int num) { return *this += (-num); }

	// Func(dice_face) -> dice_face
	template<typename Func, typename NewSymbol>
	[[nodiscard]] auto transformed(Func &&func, std::shared_ptr<indexer<NewSymbol>> new_indexer) const {
		auto current_faces = to_faces_with_prob();
		std::vector<std::pair<dice_face<NewSymbol>, double>> new_faces;
		new_faces.reserve(current_faces.size());
		for (auto &&[face, prob] : current_faces)
			new_faces.emplace_back(func(face), prob);
		return dice<NewSymbol>(new_faces.begin(), new_faces.end(), new_indexer);
	}
	template<typename Func>
	[[nodiscard]] auto transformed(Func &&func) const requires reroll_func<Func, SymbolType> { return transformed(std::forward<decltype(func)>(func), idx); }

	template<range_of<dice> Range>
	static auto get_combination_indexer(const Range &range, bool include_source_dice = false) {
		auto ret = std::make_shared<indexer<combination_dice_face_symbol<dice>>>();
		for (auto &&d : range)
			for (auto &&[face, prob] : d.to_faces_with_prob())
				ret->add_symbol(std::make_pair(face, include_source_dice ? &d : nullptr));
		return ret;
	}

	template<typename CombIndexer>
	auto to_combination_dice(std::shared_ptr<CombIndexer> comb_indexer) const {
		bool include_source_dice = bool(comb_indexer->from_index(0).second);
		return transformed([this, include_source_dice](const dice_face<symbol> &face) {
			return make_dice_face(std::make_pair(face, include_source_dice ? this : nullptr));
		}, comb_indexer);
	}

	auto roll() const {
		double s = get_random();
		accumulator<> sum;
		for (auto &&[mn, coef] : pl.terms)
			if (double(sum += coef) - s >= -eps)
				return to_face(mn);
		// shouldn't be able to reach here
		return to_face(pl.terms.begin()->first);
	}

	auto get_result_statistics(const symbol &symbol) const {
		result_statistics ret;
		int index = idx->get_index(symbol);
		for (auto &&[mn, coef] : pl.terms) {
			int value = mn[index];
			ret.min = std::min(ret.min, value);
			ret.max = std::max(ret.max, value);
			ret.E += value * coef;
			ret.E2 += value * value * coef;
		}
		ret.V = ret.E2 - ret.E * ret.E;
		ret.D = std::sqrt(ret.V);
		return ret;
	}

	// Func(symbol) -> Outputable
	template<typename Func = std::identity>
	void debug_log(Func &&func = {}) const;

	// re-roll Func(dice_face) -> std::optional<dice>
	// returned dice must be the same type and use the same indexer
	template<typename Func>
	[[nodiscard]] dice explode(Func &&func, int recursive_limit = 10) const requires explode_func<Func, SymbolType> {
		return explode_poly([&](const mono &mn) -> std::optional<poly> {
			auto &&opt_dice = func(to_face(mn));
			if (!opt_dice)
				return std::nullopt;
			return opt_dice->pl;
		}, recursive_limit);
	}

private:
	template<Hashable>
	friend class dice;

	dice(const poly &pl, std::shared_ptr<indexer<symbol>> input_idx) : pl(pl), idx(std::move(input_idx)) {}

	auto to_face(const mono &mn) const {
		dice_face<symbol> ret;
		ret.reserve(mn.dimension());
		for (int i = 0; i < mn.dimension(); i++)
			if (mn[i])
				ret[idx->from_index(i)] += mn[i];
		return ret;
	}

	// re-roll Func(mono) -> std::optional<Polynomial>
	template<typename Func>
	[[nodiscard]]
	dice explode_poly(Func &&func, int recursive_limit) const {
		// <to_be_rolled, <cur, prob>>
		unordered_map<std::optional<poly>, unordered_map<poly, accumulator<>>> ret;
		ret[pl][poly::identity()] = 1.0;

		for (int i = 0; i <= recursive_limit; i++) {
			decltype(ret) new_ret;
			for (auto &&[opt_dice, poly_to_prob] : ret)
				if (!opt_dice)
					if (new_ret.contains(opt_dice))
						for (auto &&[pl, prob] : poly_to_prob)
							new_ret[opt_dice][pl] += prob;
					else
						new_ret.emplace(opt_dice, std::move(poly_to_prob));
				else {
					for (auto &&[mn, coef] : opt_dice->terms) {
						const auto &to_be_rolled = func(mn);
						for (auto &&[pl, prob] : poly_to_prob)
							new_ret[to_be_rolled][pl * mn] += prob * coef;
					}
				}
			ret = std::move(new_ret);
		}

		std::vector<std::pair<dice, double>> vec;
		for (auto &poly_to_prob : ret | std::views::values)
			for (auto &&[pl, prob] : poly_to_prob)
				vec.emplace_back(dice(pl, idx), prob);
		return merge_dice(vec);
	}

	std::shared_ptr<indexer<symbol>> idx;
	poly pl;
};

// Func(symbol) -> Outputable
template<Hashable Symbol>
template<typename Func>
void dice<Symbol>::debug_log(Func &&func) const {
	for (auto &&[mn, coef] : get_ordered(pl.terms)) {
		bool first = true;
		std::cout << "{ ";
		for (int i = 0; i < mn.dimension(); i++)
			if (mn[i]) {
				if (!first)
					std::cout << " + ";
				first = false;
				std::cout << mn[i] << '*' << func(idx->from_index(i));
			}
		std::cout << " }: " << coef * 100 << "%" << std::endl;
	}
}

template<typename Symbol, typename... Args>
auto make_dice(Args&&... args) {
	indexer<Symbol> i;
	auto add_to_indexer = [&i](const auto &face) { for (auto &&[symbol, count] : face) i.add_symbol(symbol); };
	(add_to_indexer(args), ...);
	return dice<Symbol>(std::vector{ std::forward<Args>(args)... } | std::views::transform([](auto &&face) { return std::make_pair(face, 1.0 / sizeof...(Args)); }), std::make_shared<indexer<Symbol>>(std::move(i)));
}

template<std::input_iterator Iter>
auto make_number_dice(Iter begin, Iter end) {
	int n = int(std::distance(begin, end));
	auto view = std::ranges::subrange(begin, end) | std::views::transform([n](int idx) {
		return std::make_pair(std::views::single(std::make_pair(number_dice_symbol{}, idx)), 1.0 / n);
	});
	return number_dice(view.begin(), view.end(), number_dice_indexer());
}
auto make_number_dice(const std::ranges::input_range auto &range) { return make_number_dice(std::ranges::begin(range), std::ranges::end(range)); }

template<typename... Args>
auto make_number_dice(Args... numbers) { return make_number_dice(std::vector{ numbers... }); }

inline auto d(int n) { return make_number_dice(std::views::iota(1, n + 1)); };

inline const auto d2 = d(2);
inline const auto d3 = d(3);
inline const auto d4 = d(4);
inline const auto d6 = d(6);
inline const auto d8 = d(8);
inline const auto d10 = d(10);
inline const auto d12 = d(12);
inline const auto d20 = d(20);
inline const auto d30 = d(30);
inline const auto d60 = d(60);
inline const auto d100 = d(100);
inline const auto d120 = d(120);

inline const auto fudge = make_dice<char>(make_dice_face('+'), make_dice_face('+'), make_dice_face(' '), make_dice_face(' '), make_dice_face('-'), make_dice_face('-'));

}
