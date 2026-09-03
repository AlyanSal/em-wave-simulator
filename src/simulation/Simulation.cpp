#include <cmath>

#include "Simulation.hpp"

namespace em::sim {

simulation::simulation(const float frequency, const float width,
                       const float length, const float largest_eps)
    : width_{width},
      // length_{length},
      largest_eps_{largest_eps},
      frequency_{frequency},
      smallest_wavelength_{(math::constants::c0 / std::sqrt(largest_eps)) /
                           frequency},
      rows_{static_cast<std::size_t>(width / smallest_wavelength_) *
            math::constants::divs},
      // cols_{(length / smallest_wavelength_) * 10},
      dx_{width / static_cast<float>(rows_)},
      // dy_{length / cols_},
      dt_{dx_ / (2 * math::constants::c0)},
      cells_{rows_},
      last_two(0.0f, 0.0f),
      grid_(cells_) {}

auto simulation::setup_simulation(
    std::function<float(float)>& permittivity_function,
    std::function<float(float)>& conductivity_function,
    std::function<float(float)>& chi1_function,
    std::function<float(float)>& t0_function) -> void {
  kernel::initializeFieldConditions(
      permittivity_function, conductivity_function, chi1_function, t0_function,
      grid_.Permittivity(), grid_.Conductance(), grid_.Chi_1(), grid_.T_0(),
      dx_);

  kernel::precomputeDEISCoefficients(
      grid_.Permittivity(), grid_.Conductance(), grid_.Chi_1(), grid_.T_0(),
      grid_.ECoeff(), grid_.ICoeff(), grid_.SMCoeff(), grid_.SDCoeff(), dt_);
}

auto simulation::step_simulation() -> void {
  handle_sources();

  kernel::calculateFutureDEISFields(
      grid_.Dfield(), grid_.Efield(), grid_.Ifield(), grid_.Sfield(),
      grid_.Hfield(), grid_.ECoeff(), grid_.ICoeff(), grid_.SMCoeff(),
      grid_.SDCoeff(), dt_, dx_);

  kernel::calculateFutureHField(grid_.Hfield(), grid_.Efield(), dt_, dx_);

  kernel::updateFrequencyDomain(grid_.RealE(), grid_.ImagE(), grid_.Efield(),
                                dt_ * static_cast<float>(timestep_),
                                frequency_);

  kernel::applyBoundaryCondition(grid_.Efield(), last_two);

  ++timestep_;
}

auto simulation::handle_sources() -> void {
  const float time = static_cast<float>(timestep_) * dt_;
  const float omega_t = 2.0f * math::constants::pi * frequency_ * time;

  for (auto const& [src, idx] : sources_) {
    grid_.Dfield()[idx] += src(omega_t);
  }
}

auto simulation::add_source(source::Source source, const float pos_x) -> void {
  sources_.emplace_back(std::move(source),
                        static_cast<std::size_t>(pos_x / dx_));
}

} // namespace em::sim
