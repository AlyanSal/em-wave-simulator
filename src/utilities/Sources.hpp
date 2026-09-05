#pragma once

#include <cmath>
#include <functional>

namespace em::source {

inline constexpr auto sin1D{[](float time) -> float { return std::sin(time); }};

inline constexpr auto cos1D{[](float time) -> float { return std::cos(time); }};

} // namespace em::source
