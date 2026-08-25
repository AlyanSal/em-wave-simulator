#pragma once

#include <iostream>
#include <vector>

#include "Vector.hpp"

namespace em::memory {

/**
 * Memory handler for the simulation's memory
 *
 * TODO: Make allocator ingest dimensions, and grid sizes will be determined
 * based on frequency
 */
class allocator {
public:
  allocator(std::size_t rows, std::size_t cols);

  [[nodiscard]]
  std::pair<std::size_t, std::size_t> dimensions() const;

private:
  const std::size_t rows_;
  const std::size_t cols_;

  std::vector<vector> Efield_;
  std::vector<vector> Hfield_;
};

} // namespace em::memory