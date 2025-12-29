#pragma once

#include <complex>
#include <vector>
#include <span>

namespace spaceless::detail {

using vec_complex = std::vector<std::complex<double>>;
vec_complex fft(const vec_complex &input, std::span<const int64_t> sizes);
vec_complex ifft(const vec_complex &input, std::span<const int64_t> sizes);
vec_complex dot(const vec_complex &lhs, const vec_complex &rhs);

}
