#pragma once

#include "Grid.hpp"

namespace em::simulation {

class simulation {
public:
  simulation(float frequency, float width, float length);

  void step_simulation();

private:
  const float width_;
  // const float length_;
  const float frequency_;
  const float wavelength_;

  const float dx_;
  // const float dy_;
  const float dt_;

  memory::grid grid_;
};

} // namespace em::simulation
