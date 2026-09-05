#include "Simulation.hpp"

namespace em::sim {

simulation::simulation(const float frequency, const float width,
                       const float length, const float largest_eps)
    : width_{width},
      length_{length},
      largest_eps_{largest_eps},
      frequency_{frequency},
      smallest_wavelength_{(math::constants::c0 / std::sqrt(largest_eps)) /
                           frequency},
      dx_{static_cast<float>(math::constants::divs) / smallest_wavelength_},
      dt_{dx_ / (2 * math::constants::c0)},
      rows_{static_cast<std::size_t>(width / dx_)},
      cols_{static_cast<std::size_t>(length / dx_)},
      cells_{rows_ * cols_},
      last_two(0.0f, 0.0f),
      grid_(rows_ * cols_) {}

auto simulation::setup_simulation(
    std::function<float(float, float)>& permittivity_function,
    std::function<float(float, float)>& conductivity_function,
    std::function<float(float, float)>& chi1_function,
    std::function<float(float, float)>& t0_function) -> void {
  kernel::initializeFieldConditions(
      cols_, rows_, permittivity_function, conductivity_function, chi1_function,
      t0_function, grid_.Permittivity(), grid_.Conductance(), grid_.Chi_1(),
      grid_.T_0(), dx_);

  kernel::precomputeDEISCoefficients(
      cells_, grid_.Permittivity(), grid_.Conductance(), grid_.Chi_1(),
      grid_.T_0(), grid_.ECoeff(), grid_.ICoeff(), grid_.SMCoeff(),
      grid_.SDCoeff(), dt_);
}

auto simulation::step_simulation() -> void {
  handle_sources();

  kernel::calculateFutureDEISFields(
      cols_, rows_, grid_.Dfield(), grid_.Ezfield(), grid_.Ifield(),
      grid_.Sfield(), grid_.Hxfield(), grid_.Hyfield(), grid_.ECoeff(),
      grid_.ICoeff(), grid_.SMCoeff(), grid_.SDCoeff(), dt_, dx_);

  kernel::calculateFutureHField(cols_, rows_, grid_.Hxfield(), grid_.Hyfield(),
                                grid_.Ezfield(), dt_, dx_);

  kernel::updateFrequencyDomain(
      cells_, grid_.RealE(), grid_.ImagE(), grid_.Ezfield(),
      dt_ * static_cast<float>(timestep_), frequency_);

  kernel::applyBoundaryCondition(grid_.Ezfield(), last_two);

  ++timestep_;
}

auto simulation::handle_sources() -> void {
  const float time = static_cast<float>(timestep_) * dt_;
  const float omega_t = 2.0f * math::constants::pi * frequency_ * time;

  for (auto const& [src, idx] : sources_) {
    grid_.Dfield()[idx] += src(omega_t);
  }
}

auto simulation::add_source(std::function<float(float)> source,
                            const float pos_x) -> void {
  sources_.emplace_back(std::move(source),
                        static_cast<std::size_t>(pos_x / dx_));
}

} // namespace em::sim
