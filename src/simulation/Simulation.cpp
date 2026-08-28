#include <cmath>

#include "Simulation.hpp"

namespace em::sim {

simulation::simulation(const float frequency, const float width,
                       const float length)
    : width_{width},
      // length_{length},
      frequency_{frequency},
      wavelength_{math::constants::c0 / frequency},
      rows_{static_cast<std::size_t>(width / wavelength_) *
            math::constants::divs},
      // cols_{(length / wavelength_) * 10},
      dx_{width / static_cast<float>(rows_)},
      // dy_{length / cols_},
      dt_{dx_ / (2 * math::constants::c0)},
      cells_{rows_},
      last_two(0.0f, 0.0f),
      timestep_{},
      grid_(cells_) {}

void simulation::setup_simulation() {
  auto const max_amp{25};
  auto& Efield{grid_.Efield()};
  auto& eps{grid_.Eps()};

  auto midpoint{Efield.size() / 2};

  for (auto i{Efield.size() * 3 / 4}; i < Efield.size(); ++i) {
    eps[i] = i - (Efield.size() * 3 / 4) + 1;
  }

  const float sigma = 15.0f;

  for (auto i{0uz}; i < Efield.size(); ++i) {
    auto dist = (static_cast<float>(i) - static_cast<float>(midpoint));
    Efield[i] = max_amp * std::exp(-(dist * dist) / (2.0f * sigma * sigma));
  }

  kernel::calculateFutureHField(grid_.Hfield(), grid_.Efield(), dt_, dx_);
}

void simulation::step_simulation() {
  kernel::calculateFutureEField(grid_.Efield(), grid_.Hfield(), grid_.Eps(),
                                dt_, dx_);
  kernel::calculateFutureHField(grid_.Hfield(), grid_.Efield(), dt_, dx_);

  kernel::applyBoundaryCondition(grid_.Efield(), last_two);

  //   apply_hard_source();
}

void simulation::apply_hard_source() {}

} // namespace em::sim
