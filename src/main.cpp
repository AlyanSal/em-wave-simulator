#include "Allocator.hpp"

int main() {
  const auto FRQ{2.4e9f};
  const auto ROWS{1024uz};
  const auto COLS{1024uz};

  em::memory::allocator memory(ROWS, COLS);

  return 0;
}