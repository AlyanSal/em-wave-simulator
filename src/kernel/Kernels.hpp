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
    const std::size_t Nx, const std::size_t Ny, std::vector<float>& Dz,
    std::vector<float>& Ez, std::vector<float>& Iz, std::vector<float>& Sz,
    const std::vector<float>& Hx, const std::vector<float>& Hy,
    const std::vector<float>& e_den_coeff,
    const std::vector<float>& i_mult_coeff,
    const std::vector<float>& s_mult_coeff,
    const std::vector<float>& s_decay_coeff, const float del_t,
    const float del_x) noexcept -> void {
  const float curl_coeff{del_t / (del_x * math::constants::eps0)};

  for (auto j{1uz}; j < Ny; ++j) {
    for (auto i{1uz}; i < Nx; ++i) {
      const std::size_t idx{j * Nx + i};
      const std::size_t idx_im1{idx - 1};
      const std::size_t idx_jm1{idx - j};

      Dz[idx] +=
          curl_coeff * ((Hy[idx] - Hy[idx_im1]) - (Hx[idx] - Hx[idx_jm1]));

      Ez[idx] += (Dz[idx] - Iz[idx] - (s_decay_coeff[idx] * Sz[idx])) /
                 e_den_coeff[idx];

      Iz[idx] += (i_mult_coeff[idx] + Ez[idx]);

      Sz[idx] = (s_decay_coeff[idx] * Sz[idx]) + (s_mult_coeff[idx] * Ez[idx]);
    }
  }
}

/**
 * Calculates the E-Field vectors for a time step into the future from its
 * current state
 */
inline auto calculateFutureHField(const std::size_t Nx, const std::size_t Ny,
                                  std::vector<float>& Hx,
                                  std::vector<float>& Hy,
                                  std::vector<float>& Ez, const float del_t,
                                  const float del_x) noexcept -> void {
  const float factor = del_t / (del_x * math::constants::mu0);
  for (auto j{0uz}; j < Ny - 1; ++j) {
    for (auto i{0uz}; i < Nx - 1; ++i) {
      const std::size_t idx{j * Nx + i};
      const std::size_t idx_ip1{idx + 1};
      const std::size_t idx_jp1{idx + j};

      Hx[idx] += factor * (Ez[idx] - Ez[idx_jp1]);

      Hy[idx] += factor * (Ez[idx_ip1] - Ez[idx]);
    }
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
    const std::size_t Nx, const std::size_t Ny, // NOLINT
    std::function<float(float, float)>& permittivity_distribution,
    std::function<float(float, float)>& conductivity_distribution,
    std::function<float(float, float)>& chi_1_distribution,
    std::function<float(float, float)>& t_0_distribution,
    std::vector<float>& permittivity, std::vector<float>& conductivity,
    std::vector<float>& chi_1, std::vector<float>& t_0,
    const float dx_) noexcept -> void {
  for (auto j{0uz}; j < Ny; ++j) {
    const float pos_x{static_cast<float>(j) * dx_};
    for (auto i{0uz}; i < Nx; ++i) {
      const float pos_y{static_cast<float>(i) * dx_};
      permittivity[i] = permittivity_distribution(pos_x, pos_y);
      conductivity[i] = conductivity_distribution(pos_x, pos_y);
      chi_1[i] = chi_1_distribution(pos_x, pos_y);

      const float t0_val{t_0_distribution(pos_x, pos_y)};
      t_0[i] = (t0_val == 0.0f) ? 1.0f : t0_val;
    }
  }
}

/**
 * Precomputes the coefficients necessary for the D, E, and I
 * ElectroMagneticFields
 */
inline auto precomputeDEISCoefficients(
    const std::size_t N, // NOLINT
    const std::vector<float>& permittivity,
    const std::vector<float>& conductivity, const std::vector<float>& chi_1,
    const std::vector<float>& t_0, std::vector<float>& e_den_coeff,
    std::vector<float>& i_mult_coeff, std::vector<float>& s_mult_coeff,
    std::vector<float>& s_decay_coeff, const float del_t) noexcept -> void {
  for (auto i{0uz}; i < N; ++i) {
    const float sig_dt_eps0{(conductivity[i] * del_t) / math::constants::eps0};
    const float dt_over_t0{del_t / t_0[i]};

    i_mult_coeff[i] = sig_dt_eps0;
    s_decay_coeff[i] = std::exp(-dt_over_t0);
    s_mult_coeff[i] = chi_1[i] * dt_over_t0;
    e_den_coeff[i] = permittivity[i] + sig_dt_eps0 + s_mult_coeff[i];
  }
}

inline auto updateFrequencyDomain(const std::size_t N, // NOLINT
                                  std::vector<float>& real_E,
                                  std::vector<float>& imag_E,
                                  const std::vector<float>& EField,
                                  const float current_time,
                                  const float target_frequency) noexcept
    -> void {
  const float omega_t{2.0f * math::constants::pi * target_frequency *
                      current_time};

  const float cos_val{std::cos(omega_t)};
  const float sin_val{std::sin(omega_t)};

  for (auto i{0uz}; i < N; ++i) {
    real_E[i] += EField[i] * cos_val;
    imag_E[i] -= EField[i] * sin_val;
  }
}

} // namespace em::kernel
