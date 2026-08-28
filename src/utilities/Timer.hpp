#pragma once

#include <chrono>
#include <iostream>

namespace em::util {

inline void timer(auto&& callable) {
  const auto start{std::chrono::steady_clock::now()};
  callable();
  const auto end{std::chrono::steady_clock::now()};

  const auto duration_ns{
      std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
          .count()};

  std::print(std::cout, "Time Elapsed: {}ns\n", duration_ns);
}

} // namespace em::util
