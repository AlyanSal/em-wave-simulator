#pragma once

#include <cmath>
#include <functional>

namespace em::source {
using Source = std::move_only_function<float(float) const>;

inline auto const sin1D = [](float t) -> float { return std::sin(t); };
inline auto const cos1D = [](float t) -> float { return std::cos(t); };
} // namespace em::source
