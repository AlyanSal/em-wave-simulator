#pragma once

#include <functional>
#include <vector>

#include "Constants.hpp"

namespace em::kernel {

/**
 * Calculates the H-Field vectors for a time step into the future from its
 * current state
 */
inline auto calculateFutureDEISFields(
    std::vector<float>& DField, std::vector<float>& EField,
    std::vector<float>& IField, std::vector<float>& SField,
    const std::vector<float>& HField, const std::vector<float>& e_den_coeff,
    const std::vector<float>& i_mult_coeff,
    const std::vector<float>& s_mult_coeff,
    const std::vector<float>& s_decay_coeff, const float del_t,
    const float del_x) noexcept -> void {
  const float curl_coeff{del_t / (del_x * math::constants::eps0)};

  for (auto i{1uz}; i < DField.size(); ++i) {
    DField[i] -= curl_coeff * (HField[i] - HField[i - 1]);

    EField[i] = (DField[i] - IField[i] - (s_decay_coeff[i] * SField[i])) /
                e_den_coeff[i];

    IField[i] = IField[i] + (i_mult_coeff[i] * EField[i]);

    SField[i] = (s_decay_coeff[i] * SField[i]) + (s_mult_coeff[i] * EField[i]);
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

/**
 * Initializes the fields permittivity and conductivity distributions
 */
inline auto initializeFieldConditions(
    std::function<float(float)>& permittivity_distribution,
    std::function<float(float)>& conductivity_distribution,
    std::function<float(float)>& chi_1_distribution,
    std::function<float(float)>& t_0_distribution,
    std::vector<float>& permittivity, std::vector<float>& conductivity,
    std::vector<float>& chi_1, std::vector<float>& t_0,
    const float dx_) noexcept -> void {
  for (auto i{0uz}; i < permittivity.size(); ++i) {
    const float pos{static_cast<float>(i) * dx_};
    permittivity[i] = permittivity_distribution(pos);
    conductivity[i] = conductivity_distribution(pos);
    chi_1[i] = chi_1_distribution(pos);

    const float t0_val{t_0_distribution(pos)};
    t_0[i] = (t0_val == 0.0f) ? 1.0f : t0_val;
  }
}

/**
 * Precomputes the coefficients necessary for the D, E, and I
 * ElectroMagneticFields
 */
inline auto precomputeDEISCoefficients(
    const std::vector<float>& permittivity,
    const std::vector<float>& conductivity, const std::vector<float>& chi_1,
    const std::vector<float>& t_0, std::vector<float>& e_den_coeff,
    std::vector<float>& i_mult_coeff, std::vector<float>& s_mult_coeff,
    std::vector<float>& s_decay_coeff, const float del_t) noexcept -> void {
  for (auto i{0uz}; i < permittivity.size(); ++i) {
    const float sig_dt_eps0{(conductivity[i] * del_t) / math::constants::eps0};
    const float dt_over_t0{del_t / t_0[i]};

    i_mult_coeff[i] = sig_dt_eps0;
    s_decay_coeff[i] = std::exp(-dt_over_t0);
    s_mult_coeff[i] = chi_1[i] * dt_over_t0;
    e_den_coeff[i] = permittivity[i] + sig_dt_eps0 + s_mult_coeff[i];
  }
}

inline auto updateFrequencyDomain(std::vector<float>& real_E,
                                  std::vector<float>& imag_E,
                                  const std::vector<float>& EField,
                                  const float current_time,
                                  const float target_frequency) noexcept
    -> void {
  const float omega_t{2.0f * math::constants::pi * target_frequency *
                      current_time};

  const float cos_val{std::cos(omega_t)};
  const float sin_val{std::sin(omega_t)};

  for (auto i{0uz}; i < EField.size(); ++i) {
    real_E[i] += EField[i] * cos_val;
    imag_E[i] -= EField[i] * sin_val;
  }
}

} // namespace em::kernel
