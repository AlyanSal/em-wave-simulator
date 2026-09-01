#pragma once

#include <iostream>
#include <vector>

#include "Constants.hpp"

namespace em::memory {

/**
 * Memory handler for the simulation's memory
 */
class grid {
public:
  grid(std::size_t cells);

  auto Efield() -> std::vector<float>& { return Exfield_; }
  auto Hfield() -> std::vector<float>& { return Hyfield_; }
  auto Permittivity() -> std::vector<float>& { return permittivity_; }
  auto Conductance() -> std::vector<float>& { return conductance_; }
  auto DampCoeff() -> std::vector<float>& { return damp_coeff_; }
  auto SourceCoeff() -> std::vector<float>& { return source_coeff_; }

private:
  const std::size_t cells_;
  // const std::size_t cols_;

  std::vector<float> Exfield_;
  std::vector<float> Hyfield_;
  std::vector<float> permittivity_;
  std::vector<float> conductance_;
  std::vector<float> damp_coeff_;
  std::vector<float> source_coeff_;
};

} // namespace em::memory