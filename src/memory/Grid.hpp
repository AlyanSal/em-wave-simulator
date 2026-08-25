#pragma once

#include <iostream>
#include <vector>

#include "Constants.hpp"
#include "Vector.hpp"

namespace em::memory {

/**
 * Memory handler for the simulation's memory
 */
class grid {
public:
  grid(std::size_t cells);

private:
  const std::size_t cells_;
  // const std::size_t cols_;

  std::vector<vector> Exfield_;
  std::vector<vector> Hyfield_;
};

} // namespace em::memory