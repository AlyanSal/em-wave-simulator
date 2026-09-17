#pragma once

#include <cmath>
#include <functional>

namespace em::source {

inline constexpr auto sin1D{[](float time) -> float { return std::sin(time); }};

inline constexpr auto cos1D{[](float time) -> float { return std::cos(time); }};

inline constexpr auto gaussian{[](float omega_t) -> float {
  const float tau{(omega_t - 12.0f) / 3.0f};
  return std::exp(-tau * tau);
}};

} // namespace em::source
