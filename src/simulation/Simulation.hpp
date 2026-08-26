#pragma once

#include <cmath>

#include "Constants.hpp"
#include "Grid.hpp"
#include "Kernels.hpp"

namespace em::sim {

class simulation {
public:
  simulation(float frequency, float width, float length);

  void setup_simulation();

  void step_simulation();

  void apply_hard_source();

  std::vector<float>& getE() { return grid_.Efield(); }
  std::vector<float>& getH() { return grid_.Hfield(); }

private:
  const float width_;
  // const float length_;
  const float frequency_;
  const float wavelength_;

  const std::size_t rows_;
  // const std::size_t cols_;

  const float dx_;
  // const float dy_;
  const float dt_;

  const std::size_t cells_;

  memory::grid grid_;
};

} // namespace em::sim
