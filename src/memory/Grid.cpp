#include "Grid.hpp"

namespace em::memory {

/**
 * Constructor for the allocator
 */
grid::grid(const std::size_t rows, const std::size_t cols)
    : rows_{rows},
      cols_{cols},
      cells_{rows * cols},
      mainField_(cells_, MainFieldDefaults()),
      pmlx_(cols_, PMLDefaults()),
      pmly_(rows_, PMLDefaults()) {}

} // namespace em::memory