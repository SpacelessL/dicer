# Dicer

Dicer is a **WIP** C++23 library tailored for computing dice-rolling probability distributions for various board games. Utilizing NDFFT (ONEMKL) for accelerated performance, it also features two classes - monomial and polynomial - that could be beneficial for specific tasks. Not only does it support all the functionalities found in [AnyDice](https://anydice.com/), but it also boasts greater speed and a wider variety of features, such as support for multiple symbol types.

## Usage

### Basic Dice Rolling
```cpp
#include "dicer.hpp"
using namespace spaceless;

// Create a symbol set for a standard d6 (single numeric dimension)
auto symbols = make_simple_static_symbol_set<"N">();  // N for numeric

// Create a d6: each face (1-6) has equal probability
auto d6 = make_simple_dice(symbols, "1,2,3,4,5,6");

// Roll multiple dice: 2d6
auto two_d6 = d6 * 2;  // or: d6 + d6

// Get statistics
auto stats = two_d6.get_statistics(0);
// stats.min = 2, stats.max = 12
// stats.E ≈ 7.0 (expected value)
// stats.V ≈ 5.83 (variance)
// stats.D ≈ 2.42 (standard deviation)
```

### Custom Symbol Dice
```cpp
// Create dice with multiple symbol types
// S = sword, M = magic, D = defense
auto combat_symbols = make_simple_static_symbol_set<"SMD">();

// Combat die: faces can have multiple symbols
// Format: "SS" means 2 swords, "SD" means 1 sword + 1 defense
// Alternative format: "2S,1S1D,1S,2M,1M,1D" (explicit counts)
auto combat_die = make_simple_dice(combat_symbols, "0,0,SS,SD,S,MM,M,D");

// Defense die
auto defense_die = make_simple_dice(combat_symbols, "0,D,D,DD,DD,DDD");

// Roll 3 combat dice
auto three_rolls = combat_die * 3;

// Get statistics for swords (index 0 = 'S')
auto sword_stats = three_rolls.get_statistics(0);  // E[swords]

// Sample a random roll
auto roll = three_rolls.sample();  // returns a monomial with symbol counts

// Supports large numbers of dice with FFT acceleration
auto massive_roll = combat_die * 1000 + defense_die * 500;
```

### Polynomials (Generating Functions)
```cpp
// Create polynomials directly
polynomial<1> p;  // single-variable polynomial
p += monomial<1>({2});   // x^2
p += monomial<1>({1});   // x
p -= 1.0;                // -1
// p = x^2 + x - 1

// Polynomial arithmetic
polynomial<1> a({{{1}, 1}, {{0}, -1}});  // x - 1
polynomial<1> b({{{2}, 1}, {{0}, -1}});  // x^2 - 1

auto product = a * b;      // (x-1)(x^2-1)
auto quotient = b / a;     // (x^2-1)/(x-1) = x + 1
auto remainder = b % a;    // 0

// Multi-variable polynomials
polynomial<2> mv;  // two variables (x, y)
mv += monomial<2>({2, 1});  // x^2 * y
mv += monomial<2>({1, 2});  // x * y^2
// mv = x^2*y + x*y^2

std::cout << mv << '\n';  // prints: x^2y+xy^2
```

### Monomials
```cpp
// Static size monomials (fixed number of variables at compile time)
monomial<3> m1({2, 1, 0});  // x^2 * y^1 * z^0
monomial<3> m2({1, 1, 1});  // x * y * z

auto product = m1 * m2;  // x^3 * y^2 * z
auto inverse = m1.inverse();  // x^(-2) * y^(-1)

// Dynamic size monomials
monomial<DYNAMIC> dm({1, 2, 3});
dm.add_exponent(0, 1);  // increment first exponent
```

### Symbol Sets
```cpp
// Static symbol set (compile-time fixed)
auto static_set = make_simple_static_symbol_set<"ABC">();
static_set->get_index('A');  // returns 0
static_set->get_index('B');  // returns 1
static_set->from_index(2);   // returns 'C'

// Dynamic symbol set (runtime extensible)
dynamic_symbol_set<char> dyn_set({'X', 'Y'});
dyn_set.add_symbol('Z');     // returns 2
dyn_set.contains('X');       // true
dyn_set.size();              // 3
```

### Distribution Transformations
```cpp
auto d6 = make_simple_dice(symbols, "1,2,3,4,5,6");
auto two_d6 = d6 * 2;

// Transform outcomes (e.g., count successes where value >= 4)
auto successes = two_d6.transformed([](const auto& mn) {
    monomial<1> result;
    result.set_exponent(0, mn[0] >= 4 ? 1 : 0);  // 1 if >= 4, else 0
    return result;
});
```

## Todo
1.  Compare performance differences between `std::pmr::vector w/ std::array` and other containers like `std::vector`, or self-implemented array-based vector.
2.  Improve current implementations, especially the helper functions.
3.  Get more helper functions, particularly for reroll and explode mechanics.
4.  Support calculated complex explode probabilities, where dice can explode into arbitrary combinations of symbols and other dice.
5.  Add support for proper graphic output.
6.  Finally, implement DSL support to enable users to easily create their calculation scripts.

# License

This project is licensed under the MIT License - see the LICENSE file for details.
