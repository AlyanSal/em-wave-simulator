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

  std::vector<float>& Efield() { return Exfield_; }
  std::vector<float>& Hfield() { return Hyfield_; }
  std::vector<float>& Eps() { return eps_; }

private:
  const std::size_t cells_;
  // const std::size_t cols_;

  std::vector<float> Exfield_;
  std::vector<float> eps_;
  std::vector<float> Hyfield_;
};

} // namespace em::memory