#pragma once

#include <cstdint>

#include "Buffer.hpp"

namespace em::memory {

template <typename T, std::size_t num_arrays> class SoA {
private:
  std::size_t count_;
  Buffer<T> buffer_;

public:
  explicit SoA(const std::size_t count,
               const std::array<T, num_arrays>& defaults = {})
      : count_{count},
        buffer_(count * num_arrays) {
    for (auto i{0uz}; i < num_arrays; ++i) {
      if (defaults[i] != T{0}) {
        std::fill_n((*this)[i], count_, defaults[i]);
      }
    }
  }

  [[nodiscard]]
  auto constexpr operator[](std::size_t array) noexcept -> T* {
    return buffer_.data() + (count_ * array);
  }

  [[nodiscard]]
  auto constexpr operator[](std::size_t array) const noexcept -> T* {
    return buffer_.data() + (count_ * array);
  }

  template <typename Enum>
    requires std::is_enum_v<Enum>
  [[nodiscard]] constexpr auto operator[](Enum array) const noexcept -> T* {
    return (*this)[std::to_underlying(array)];
  }
};

} // namespace em::memory
