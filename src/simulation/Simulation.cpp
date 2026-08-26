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
      grid_(cells_) {}

void simulation::setup_simulation() {
  auto const max_amp{25};
  auto& Efield{grid_.Hfield()};

  auto midpoint{Efield.size() / 2};

  const float sigma = 15.0f;

  for (auto i{0uz}; i < Efield.size(); ++i) {
    float dist = static_cast<float>(static_cast<long long>(i) -
                                    static_cast<long long>(midpoint));
    Efield[i] = max_amp * std::exp(-(dist * dist) / (2.0f * sigma * sigma));
  }

  kernel::calculateFutureEField(grid_.Efield(), grid_.Hfield(), dt_, dx_);
}

void simulation::step_simulation() {
  kernel::calculateFutureEField(grid_.Efield(), grid_.Hfield(), dt_, dx_);
  kernel::calculateFutureHField(grid_.Hfield(), grid_.Efield(), dt_, dx_);

  //   apply_hard_source();
}

void simulation::apply_hard_source() {}

} // namespace em::sim
