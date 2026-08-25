#include "Grid.hpp"

namespace em::memory {

/**
 * Constructor for the allocator
 */
grid::grid(const std::size_t cells)
    : cells_{cells},
      Exfield_(cells, vector()),
      Hyfield_(cells, vector()) {}

} // namespace em::memory