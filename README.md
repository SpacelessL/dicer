# Dicer

Dicer is a **WIP** C++23 library tailored for computing dice-rolling probability distributions for various board games. Utilizing FFT (pocketfft or ONEMKL) for accelerated performance, it also features two classes - monomial and polynomial - that could be beneficial for specific tasks. Not only does it support all the functionalities found in [AnyDice](https://anydice.com/), but it also boasts greater speed and a wider variety of features, such as support for multiple symbol types.

## Todo
1.  Compare performance differences between `std::pmr::vector w/ std::array` and other containers like `std::vector`, or self-implemented array-based vector.
2.  Improve current implementations, especially the helper functions.
3.  Get more helper functions, particularly for reroll and explode mechanics.
4.  Add support for proper graphic output.
5.  Finally, implement DSL support to enable users to easily create their calculation scripts.

# License

This project is licensed under the MIT License - see the LICENSE file for details.
