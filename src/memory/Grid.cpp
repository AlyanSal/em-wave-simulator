#include "Grid.hpp"

namespace em::memory {

/**
 * Constructor for the allocator
 */
grid::grid(const std::size_t cells)
    : cells_{cells},
      Exfield_(cells, 0.0f),
      Hyfield_(cells, 0.0f),
      permittivity_(cells, 1.0f),
      conductance_(cells, 0.0f),
      damp_coeff_(cells, 0.0f),
      source_coeff_(cells, 0.0f) {}

} // namespace em::memory