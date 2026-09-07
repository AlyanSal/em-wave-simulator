#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>

#include "Constants.hpp"

namespace em::kernel {

/**
 * Calculates the D and E field vectors (including Debye dispersion and CPML)
 * for a time step into the future.
 */
inline auto calculateFutureDEISFields(
    const std::size_t Nx, const std::size_t Ny, std::vector<float>& Dz,
    std::vector<float>& Ez, std::vector<float>& Iz, std::vector<float>& Sz,
    const std::vector<float>& Hx, const std::vector<float>& Hy,
    std::vector<float>& psi_Ezx, std::vector<float>& psi_Ezy,
    const std::vector<float>& be_x, const std::vector<float>& ae_x,
    const std::vector<float>& inv_ke_x, const std::vector<float>& be_y,
    const std::vector<float>& ae_y, const std::vector<float>& inv_ke_y,
    const std::vector<float>& e_den_coeff,
    const std::vector<float>& i_mult_coeff,
    const std::vector<float>& s_mult_coeff,
    const std::vector<float>& s_decay_coeff, const float del_t,
    const float del_x) noexcept -> void {
  const float inv_del_x{1.0f / del_x};
  const float dt_eps0{del_t / math::constants::eps0};

  for (auto j{1uz}; j < Ny - 1; ++j) {
    const float bey{be_y[j]};
    const float aey{ae_y[j]};
    const float inv_key{inv_ke_y[j]};

    for (auto i{1uz}; i < Nx - 1; ++i) {
      const std::size_t idx{(j * Nx) + i};
      const std::size_t idx_im1{idx - 1};
      const std::size_t idx_jm1{idx - Nx};

      const float dHy_x{(Hy[idx] - Hy[idx_im1]) * inv_del_x};
      const float dHx_y{(Hx[idx] - Hx[idx_jm1]) * inv_del_x};

      psi_Ezx[idx] = (be_x[i] * psi_Ezx[idx]) + (ae_x[i] * dHy_x);
      psi_Ezy[idx] = (bey * psi_Ezy[idx]) + (aey * dHx_y);

      const float curl_H{(inv_ke_x[i] * dHy_x + psi_Ezx[idx]) -
                         (inv_key * dHx_y + psi_Ezy[idx])};

      Dz[idx] += dt_eps0 * curl_H;

      Ez[idx] = (Dz[idx] - Iz[idx] - (s_decay_coeff[idx] * Sz[idx])) /
                e_den_coeff[idx];

      Iz[idx] += (i_mult_coeff[idx] * Ez[idx]);

      Sz[idx] = (s_decay_coeff[idx] * Sz[idx]) + (s_mult_coeff[idx] * Ez[idx]);
    }
  }
}

/**
 * Calculates the H-Field vectors with CPML for a time step into the future
 */
inline auto calculateFutureHField(
    const std::size_t Nx, const std::size_t Ny,
    std::vector<float>& Hx, std::vector<float>& Hy,
    const std::vector<float>& Ez,
    std::vector<float>& psi_Hxy, std::vector<float>& psi_Hyx,
    const std::vector<float>& bh_x, const std::vector<float>& ah_x,
    const std::vector<float>& inv_kh_x, const std::vector<float>& bh_y,
    const std::vector<float>& ah_y, const std::vector<float>& inv_kh_y,
    const float del_t, const float del_x) noexcept -> void {
  const float inv_del_x{1.0f / del_x};
  const float dt_mu0{del_t / math::constants::mu0};

  for (auto j{0uz}; j < Ny - 1; ++j) {
    const float bhy{bh_y[j]};
    const float ahy{ah_y[j]};
    const float inv_khy{inv_kh_y[j]};

    for (auto i{0uz}; i < Nx; ++i) {
      const std::size_t idx{(j * Nx) + i};
      const float dEz_y{(Ez[idx + Nx] - Ez[idx]) * inv_del_x};

      psi_Hxy[idx] = (bhy * psi_Hxy[idx]) + (ahy * dEz_y);
      Hx[idx] -= dt_mu0 * ((inv_khy * dEz_y) + psi_Hxy[idx]);
    }
  }

  for (auto j{0uz}; j < Ny; ++j) {
    for (auto i{0uz}; i < Nx - 1; ++i) {
      const std::size_t idx{(j * Nx) + i};
      const float dEz_x{(Ez[idx + 1] - Ez[idx]) * inv_del_x};

      psi_Hyx[idx] = (bh_x[i] * psi_Hyx[idx]) + (ah_x[i] * dEz_x);
      Hy[idx] += dt_mu0 * ((inv_kh_x[i] * dEz_x) + psi_Hyx[idx]);
    }
  }
}

/**
 * Initializes the fields permittivity and conductivity distributions with PML padding offset
 */
inline auto initializeFieldConditions(
    const std::size_t Nx, const std::size_t Ny,
    const std::size_t pml_cells_x, const std::size_t pml_cells_y,
    std::function<float(float, float)>& permittivity_distribution,
    std::function<float(float, float)>& conductivity_distribution,
    std::function<float(float, float)>& chi_1_distribution,
    std::function<float(float, float)>& t_0_distribution,
    std::vector<float>& permittivity, std::vector<float>& conductivity,
    std::vector<float>& chi_1, std::vector<float>& t_0,
    const float dx_) noexcept -> void {
  for (auto j{0uz}; j < Ny; ++j) {
    const bool in_pml_y = (j < pml_cells_y) || (j >= Ny - pml_cells_y);
    const float pos_y = in_pml_y ? 0.0f : static_cast<float>(j - pml_cells_y) * dx_;

    for (auto i{0uz}; i < Nx; ++i) {
      const bool in_pml_x = (i < pml_cells_x) || (i >= Nx - pml_cells_x);
      const float pos_x = in_pml_x ? 0.0f : static_cast<float>(i - pml_cells_x) * dx_;
      const std::size_t idx{(j * Nx) + i};

      if (in_pml_x || in_pml_y) {
        permittivity[idx] = 1.0f;
        conductivity[idx] = 0.0f;
        chi_1[idx] = 0.0f;
        t_0[idx] = 1.0f;
      } else {
        permittivity[idx] = permittivity_distribution(pos_x, pos_y);
        conductivity[idx] = conductivity_distribution(pos_x, pos_y);
        chi_1[idx] = chi_1_distribution(pos_x, pos_y);

        const float t0_val{t_0_distribution(pos_x, pos_y)};
        t_0[idx] = (t0_val == 0.0f) ? 1.0f : t0_val;
      }
    }
  }
}

/**
 * Backwards-compatible overload without PML offsets
 */
inline auto initializeFieldConditions(
    const std::size_t Nx, const std::size_t Ny,
    std::function<float(float, float)>& permittivity_distribution,
    std::function<float(float, float)>& conductivity_distribution,
    std::function<float(float, float)>& chi_1_distribution,
    std::function<float(float, float)>& t_0_distribution,
    std::vector<float>& permittivity, std::vector<float>& conductivity,
    std::vector<float>& chi_1, std::vector<float>& t_0,
    const float dx_) noexcept -> void {
  initializeFieldConditions(Nx, Ny, 0uz, 0uz, permittivity_distribution,
                            conductivity_distribution, chi_1_distribution,
                            t_0_distribution, permittivity, conductivity,
                            chi_1, t_0, dx_);
}

/**
 * Computes 1D CPML parameters (be, ae, inv_ke, bh, ah, inv_kh)
 * for a dimension with total_cells and pml_cells on each boundary.
 */
inline auto initializeCPML1DProfile(
    const std::size_t total_cells, const std::size_t pml_cells,
    std::vector<float>& be, std::vector<float>& ae, std::vector<float>& inv_ke,
    std::vector<float>& bh, std::vector<float>& ah, std::vector<float>& inv_kh,
    const float dt, const float dx) noexcept -> void {
  std::fill(be.begin(), be.end(), 0.0f);
  std::fill(ae.begin(), ae.end(), 0.0f);
  std::fill(inv_ke.begin(), inv_ke.end(), 1.0f);
  std::fill(bh.begin(), bh.end(), 0.0f);
  std::fill(ah.begin(), ah.end(), 0.0f);
  std::fill(inv_kh.begin(), inv_kh.end(), 1.0f);

  if (pml_cells == 0 || (2uz * pml_cells) >= total_cells) {
    return;
  }

  const float d{static_cast<float>(pml_cells) * dx};
  const float m{math::constants::pml_m};
  const float r0{math::constants::pml_r0};
  const float eta0{std::sqrt(math::constants::mu0 / math::constants::eps0)};
  const float sigma_max{-((m + 1.0f) * std::log(r0)) / (2.0f * eta0 * d)};
  const float kappa_max{math::constants::pml_kappa_max};
  const float alpha_max{math::constants::pml_alpha_max};

  auto calc_cpml = [&](const float dist, float& b, float& a, float& inv_k) {
    if (dist <= 0.0f) {
      return;
    }
    const float norm{dist / d};
    const float norm_m{std::pow(norm, m)};
    const float sigma{sigma_max * norm_m};
    const float kappa{1.0f + ((kappa_max - 1.0f) * norm_m)};
    const float alpha{alpha_max * (1.0f - norm)};

    b = std::exp(-((sigma / (kappa * math::constants::eps0)) +
                   (alpha / math::constants::eps0)) * dt);
    a = (sigma / (kappa * (sigma + (kappa * alpha)))) * (b - 1.0f);
    inv_k = 1.0f / kappa;
  };

  for (std::size_t k{0uz}; k < pml_cells; ++k) {
    const float dist_e{static_cast<float>(k + 1uz) * dx};
    const float dist_h{(static_cast<float>(k + 1uz) - 0.5f) * dx};

    const std::size_t left_idx{(pml_cells - 1uz) - k};
    const std::size_t right_idx{(total_cells - pml_cells) + k};

    calc_cpml(dist_e, be[left_idx], ae[left_idx], inv_ke[left_idx]);
    calc_cpml(dist_h, bh[left_idx], ah[left_idx], inv_kh[left_idx]);

    calc_cpml(dist_e, be[right_idx], ae[right_idx], inv_ke[right_idx]);
    calc_cpml(dist_h, bh[right_idx], ah[right_idx], inv_kh[right_idx]);
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

inline auto updateFrequencyDomain(
    const std::size_t N, std::vector<float>& real_E, std::vector<float>& imag_E,
    const std::vector<float>& EField, const float current_time,
    const float target_frequency) noexcept -> void {
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
