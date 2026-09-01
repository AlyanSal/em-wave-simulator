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
    std::function<float(float)>& conductivity_function) -> void {
  kernel::initializeFieldConditions(dx_, permittivity_function,
                                    conductivity_function, grid_.Permittivity(),
                                    grid_.Conductance());

  kernel::precomputeEFieldCalculationCoefficients(
      grid_.Permittivity(), grid_.Conductance(), grid_.DampCoeff(),
      grid_.SourceCoeff(), dt_, dx_);
}

auto simulation::step_simulation() -> void {
  kernel::calculateFutureEField(grid_.Efield(), grid_.Hfield(),
                                grid_.DampCoeff(), grid_.SourceCoeff());
  kernel::calculateFutureHField(grid_.Hfield(), grid_.Efield(), dt_, dx_);

  kernel::applyBoundaryCondition(grid_.Efield(), last_two);

  handle_sources();

  ++timestep_;
}

auto simulation::handle_sources() -> void {
  const float time = static_cast<float>(timestep_) * dt_;
  const float omega_t = 2.0f * std::numbers::pi_v<float> * frequency_ * time;

  for (auto const& [src, idx] : sources_) {
    grid_.Efield()[idx] += src(omega_t);
  }
}

auto simulation::add_source(source::Source source, const float pos_x) -> void {
  sources_.emplace_back(std::move(source),
                        static_cast<std::size_t>(pos_x / dx_));
}

} // namespace em::sim
