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

auto simulation::setup_simulation() -> void {
  auto& Efield{grid_.Efield()};
  auto& eps{grid_.Eps()};

  for (auto i{Efield.size() * 3 / 4}; i < Efield.size(); ++i) {
    eps[i] = 4;
  }
}

auto simulation::step_simulation() -> void {
  kernel::calculateFutureEField(grid_.Efield(), grid_.Hfield(), grid_.Eps(),
                                dt_, dx_, 0.04f);
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

auto simulation::add_source(source::Source source, const float posx) -> void {
  sources_.emplace_back(std::move(source),
                        static_cast<std::size_t>(posx / dx_));
}

} // namespace em::sim
