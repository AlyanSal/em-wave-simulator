#include "Simulation.hpp"

namespace em::sim {

using namespace mem;

simulation::simulation(const float frequency, const float width,
                       const float length, const float largest_eps,
                       const float pml_factor, const std::size_t divs)
    : width_{width},
      length_{length},
      largest_eps_{largest_eps},
      frequency_{frequency},
      smallest_wavelength_{(math::constants::c0 / std::sqrt(largest_eps)) /
                           frequency},
      pml_factor_{std::clamp(pml_factor, 0.0f, 1.0f)},
      dx_{smallest_wavelength_ / static_cast<float>(divs)},
      dt_{dx_ / (2 * math::constants::c0)},
      interior_rows_{static_cast<std::size_t>(width / dx_)},
      interior_cols_{static_cast<std::size_t>(length / dx_)},
      pml_cells_x_{pml_factor_ > 0.0f
                       ? std::max(1uz, static_cast<std::size_t>(std::round(
                                           static_cast<float>(interior_cols_) *
                                           pml_factor_ * 0.5f)))
                       : 0uz},
      pml_cells_y_{pml_factor_ > 0.0f
                       ? std::max(1uz, static_cast<std::size_t>(std::round(
                                           static_cast<float>(interior_rows_) *
                                           pml_factor_ * 0.5f)))
                       : 0uz},
      rows_{interior_rows_ + (2uz * pml_cells_y_)},
      cols_{interior_cols_ + (2uz * pml_cells_x_)},
      cells_{rows_ * cols_},
      grid_(rows_, cols_) {}

auto simulation::setup_simulation(
    std::function<float(float, float)>& permittivity_function,
    std::function<float(float, float)>& conductivity_function,
    std::function<float(float, float)>& chi1_function,
    std::function<float(float, float)>& t0_function) -> void {
  kernel::initializeFieldConditions(
      cols_, rows_, pml_cells_x_, pml_cells_y_, permittivity_function,
      conductivity_function, chi1_function, t0_function,
      grid_[MainField::Permittivity], grid_[MainField::Conductivity],
      grid_[MainField::Chi_1], grid_[MainField::T_0], dx_);

  kernel::precomputeDEISCoefficients(
      cells_, grid_[MainField::Permittivity], grid_[MainField::Conductivity],
      grid_[MainField::Chi_1], grid_[MainField::T_0],
      grid_[MainField::inv_E_Coeffs], grid_[MainField::I_Coeffs],
      grid_[MainField::SM_Coeffs], grid_[MainField::SD_Coeffs], dt_);

  kernel::initializeCPML1DProfile(cols_, pml_cells_x_, grid_[PMLX::BE],
                                  grid_[PMLX::AE], grid_[PMLX::Inv_KE],
                                  grid_[PMLX::BH], grid_[PMLX::AH],
                                  grid_[PMLX::Inv_KH], dt_, dx_);

  kernel::initializeCPML1DProfile(rows_, pml_cells_y_, grid_[PMLY::BE],
                                  grid_[PMLY::AE], grid_[PMLY::Inv_KE],
                                  grid_[PMLY::BH], grid_[PMLY::AH],
                                  grid_[PMLY::Inv_KH], dt_, dx_);
}

auto simulation::step_simulation() noexcept -> void {
  handle_sources();

  kernel::calculateFutureDEISFields(
      cols_, rows_, grid_[MainField::DzField], grid_[MainField::EzField],
      grid_[MainField::IzField], grid_[MainField::SzField],
      grid_[MainField::HxField], grid_[MainField::HyField],
      grid_[MainField::Psi_Ez_x], grid_[MainField::Psi_Ez_y], grid_[PMLX::BE],
      grid_[PMLX::AE], grid_[PMLX::Inv_KE], grid_[PMLY::BE], grid_[PMLY::AE],
      grid_[PMLY::Inv_KE], grid_[MainField::inv_E_Coeffs],
      grid_[MainField::I_Coeffs], grid_[MainField::SM_Coeffs],
      grid_[MainField::SD_Coeffs], dt_, dx_);

  kernel::calculateFutureHField(
      cols_, rows_, grid_[MainField::HxField], grid_[MainField::HyField],
      grid_[MainField::EzField], grid_[MainField::Psi_Hx_y],
      grid_[MainField::Psi_Hy_x], grid_[PMLX::BH], grid_[PMLX::AH],
      grid_[PMLX::Inv_KH], grid_[PMLY::BH], grid_[PMLY::AH],
      grid_[PMLY::Inv_KH], dt_, dx_);

  // kernel::updateFrequencyDomain(
  //     cells_, grid_[MainField::Real_E], grid_[MainField::Imag_E],
  //     grid_[MainField::EzField], dt_ * static_cast<float>(timestep_),
  //     frequency_);

  ++timestep_;
}

auto simulation::handle_sources() const noexcept -> void {
  const float time = static_cast<float>(timestep_) * dt_;
  const float omega_t = 2.0f * math::constants::pi * frequency_ * time;

  for (auto const& [src, idx] : sources_) {
    grid_[MainField::DzField][idx] += src(omega_t);
  }
}

auto simulation::add_source(std::function<float(float)> source,
                            const float pos_x, const float pos_y) noexcept
    -> void {
  const auto i{pml_cells_x_ + static_cast<std::size_t>(pos_x / dx_)};
  const auto j{pml_cells_y_ + static_cast<std::size_t>(pos_y / dx_)};
  if (i < cols_ && j < rows_) {
    sources_.emplace_back(std::move(source), (j * cols_) + i);
  }
}

} // namespace em::sim
