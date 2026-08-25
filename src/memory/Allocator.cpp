#include "Allocator.hpp"

namespace em::memory {

/**
 * Constructor for the allocator
 * @param rows the number of rows for the grid
 * @param cols the number of cols for the grid
 */
allocator::allocator(std::size_t rows, std::size_t cols)
    : Efield_(rows * cols, vector()),
      Hfield_(rows * cols, vector()),
      rows_{rows},
      cols_{cols} {}

[[nodiscard]]
std::pair<std::size_t, std::size_t> allocator::dimensions() const {
  return {rows_, cols_};
}

} // namespace em::memory
