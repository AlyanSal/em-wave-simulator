#pragma once

#include <cmath>
#include <functional>

namespace em::source {
using Source = std::move_only_function<float(float) const>;

inline auto const sin1D = [](float time) -> float { return std::sin(time); };
inline auto const cos1D = [](float time) -> float { return std::cos(time); };
} // namespace em::source
