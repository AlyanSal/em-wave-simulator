#pragma once

#include <cmath>
#include <functional>

#include "Constants.hpp"
#include "Grid.hpp"
#include "Kernels.hpp"
#include "Sources.hpp"

namespace em::sim {

class simulation {
public:
  simulation(float frequency, float width, float length, float largest_eps);

  void
  setup_simulation(std::function<float(float, float)>& permittivity_function,
                   std::function<float(float, float)>& conductivity_function,
                   std::function<float(float, float)>& chi1_function,
                   std::function<float(float, float)>& t0_function);

  void step_simulation();

  void add_source(std::function<float(float)> source, float pos_x, float pos_y);

  void handle_sources();

  auto getE() -> std::vector<float>& { return grid_.Ezfield(); }
  auto getH() -> std::vector<float>& { return grid_.Hxfield(); }

  void setEpsilonDist();
  void setConductanceDist();

private:
  const float width_;
  const float length_;
  float largest_eps_;
  const float frequency_;
  const float smallest_wavelength_;

  const float dx_;
  const float dt_;

  const std::size_t rows_;
  const std::size_t cols_;
  const std::size_t cells_;

  std::vector<std::pair<std::function<float(float)>, std::size_t>> sources_;

  std::size_t timestep_{};

  memory::grid grid_;
};

} // namespace em::sim
