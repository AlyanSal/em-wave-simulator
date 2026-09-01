#pragma once

#include <functional>
#include <vector>

#include "Constants.hpp"

namespace em::kernel {

/**
 * Calculates the H-Field vectors for a time step into the future from its
 * current state
 */
inline auto calculateFutureEField(std::vector<float>& EField,
                                  std::vector<float>& HField,
                                  std::vector<float>& damp_coeff,
                                  std::vector<float>& source_coeff) noexcept
    -> void {
  for (auto i{1uz}; i < EField.size() - 1; ++i) {
    EField[i] = (damp_coeff[i] * EField[i]) -
                (source_coeff[i] * (HField[i] - HField[i - 1]));
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

inline auto initializeFieldConditions(
    const float dx_, std::function<float(float)>& permittivity_distribution,
    std::function<float(float)>& conductivity_distribution,
    std::vector<float>& permittivity, std::vector<float>& conductance) noexcept
    -> void {
  for (auto i{0uz}; i < permittivity.size(); ++i) {
    permittivity[i] = permittivity_distribution(static_cast<float>(i) * dx_);
    conductance[i] = conductivity_distribution(static_cast<float>(i) * dx_);
  }
}

inline auto precomputeEFieldCalculationCoefficients(
    const std::vector<float>& permittivity,
    const std::vector<float>& conductivity, std::vector<float>& damp_coeff,
    std::vector<float>& source_coeff, const float del_t, const float del_x)
    -> void {
  const float electric_factor{del_t / (2 * math::constants::eps0)};

  for (auto i{0uz}; i < permittivity.size(); ++i) {
    const float rel_factor{electric_factor *
                           (conductivity[i] / permittivity[i])};
    damp_coeff[i] = (1 - rel_factor) / (1 + rel_factor);
    source_coeff[i] = del_t / (del_x * permittivity[i] * math::constants::eps0 *
                               (1.0f + rel_factor));
  }
}

} // namespace em::kernel
