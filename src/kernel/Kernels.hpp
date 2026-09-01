#pragma once

#include <vector>

#include "Constants.hpp"

namespace em::kernel {

/**
 * Calculates the H-Field vectors for a time step into the future from its
 * current state
 */
inline auto calculateFutureEField(std::vector<float>& EField,
                                  std::vector<float>& HField,
                                  std::vector<float>& eps, const float del_t,
                                  const float del_x, const float sigma) noexcept
    -> void {
  const float electric_factor{(del_t * sigma) / (2 * math::constants::eps0)};

  for (auto i{1uz}; i < EField.size() - 1; ++i) {
    const float relative_ef{electric_factor / eps[i]};
    const float damp_coeff{(1.0f - relative_ef) / (1.0f + relative_ef)};

    const float permittivity{math::constants::eps0 * eps[i]};
    const float source_coeff{del_t /
                             (del_x * permittivity * (1.0f + relative_ef))};
    EField[i] =
        (damp_coeff * EField[i]) - (source_coeff * (HField[i] - HField[i - 1]));
  }
}

/**
 * Calculates the E-Field vectors for a time step into the future from its
 * current state
 */
inline auto calculateFutureHField(std::vector<float>& HField,
                                  std::vector<float>& EField, const float del_t,
                                  const float del_x) noexcept -> void {
  const float factor = del_t / (del_x * math::constants::mu0);
  for (auto i{0uz}; i < HField.size() - 1; ++i) {
    HField[i] -= factor * (EField[i + 1] - EField[i]);
  }
}

/**
 * Approximates a boundary condition
 */
inline auto applyBoundaryCondition(std::vector<float>& Efield,
                                   std::pair<float, float>& last_two) noexcept
    -> void {
  Efield[0] = last_two.second;
  last_two.second = last_two.first;
  last_two.first = Efield[1];
}

} // namespace em::kernel
