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

  auto Efield() -> std::vector<float>& { return Efield_; }
  auto Hfield() -> std::vector<float>& { return Hfield_; }
  auto Dfield() -> std::vector<float>& { return Dfield_; }
  auto Ifield() -> std::vector<float>& { return Ifield_; }
  auto Permittivity() -> std::vector<float>& { return permittivity_; }
  auto Conductance() -> std::vector<float>& { return conductance_; }
  auto ECoeff() -> std::vector<float>& { return e_den_coeff_; }
  auto ICoeff() -> std::vector<float>& { return i_mult_coeff_; }

private:
  const std::size_t cells_;

  std::vector<float> Efield_;
  std::vector<float> Hfield_;
  std::vector<float> Dfield_;
  std::vector<float> Ifield_;
  std::vector<float> permittivity_;
  std::vector<float> conductance_;
  std::vector<float> e_den_coeff_;
  std::vector<float> i_mult_coeff_;
};

} // namespace em::memory