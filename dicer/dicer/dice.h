#pragma once

#include "distribution.h"
#include "symbol_set.h"

namespace spaceless {

template<SymbolSet S, std::signed_integral T = int8_t>
class dice {
public:
	using symbol_set = S;
	static constexpr auto N = symbol_set::N;
	using mono = monomial<N, T>;
	using poly = polynomial<N, T>;
	using dist = distribution<N, T>;

	dice() = default;
	dice(const dice &) = default;
	dice(dice &&) noexcept = default;
	dice &operator = (const dice &) = default;
	dice &operator = (dice &&) noexcept = default;

	dice(std::shared_ptr<symbol_set> s, poly p) : s(std::move(s)), d(std::move(p)) {}
	dice(std::shared_ptr<symbol_set> s, dist d) : s(std::move(s)), d(std::move(d)) {}

	operator const dist &() const { return d; }
	const dist &distribution() const { return d; }
	const symbol_set &symbols() const { return s; }

private:
	std::shared_ptr<symbol_set> s;
	dist d;
};

template<SymbolSet S, std::signed_integral T = int8_t>
auto make_simple_dice(std::shared_ptr<S> s, std::string_view text) requires std::is_same_v<typename S::symbol, char> {
	polynomial<S::N, T> pl;
	for (auto word : text | std::views::split(',')) {
		monomial<S::N, T> mn;
		auto get_num = [](std::string_view num) -> std::optional<T> {
			T ret = 0;
			auto ed = num.data() + num.length();
			auto [ptr, ec] = std::from_chars(num.data(), ed, ret);
			if (ec == std::errc{} && ptr == ed)
				return ret;
			return{};
		};
		auto process = [&](std::string_view part) {
			if (auto idx = s->get_index(part.back()); idx != -1) {
				if (part.length() == 1)
					mn.add_exponent(idx, 1);
				else {
					auto num = get_num(part.substr(0, part.length() - 1));
					ENSURE(num);
					mn.add_exponent(idx, num.value());
				}
			}
			else {
				ENSURE(S::N == 1);
				auto num = get_num(part);
				ENSURE(num);
				mn.add_exponent(0, num.value());
			}
		};
		for (int i = 0, j = -1; i < word.size(); i++)
			if (s->contains(word[i]) || i == word.size() - 1) {
				process(std::string_view(word).substr(j + 1, i - j));
				j = i;
			}
		pl += mn;
	}
	return dice(s, pl);
}

}
