#pragma once

#include <cstdint>

#include "Buffer.hpp"

namespace em::mem {

template <typename T, std::size_t num_arrays, std::size_t alignment> class SoA {
  static_assert(alignment % sizeof(T) == 0,
                "Alignment must be multiple of element size");

private:
  std::size_t count_;
  std::size_t stride_;
  Buffer<T, alignment> buffer_;

  [[nodiscard]]
  static constexpr auto calculate_stride(std::size_t count) noexcept
      -> std::size_t {
    constexpr std::size_t align_elements{alignment / sizeof(T)};

    return (count + align_elements - 1) & ~(align_elements - 1);
  }

public:
  explicit SoA(const std::size_t count,
               const std::array<T, num_arrays>& defaults = {})
      : count_{count},
        stride_{calculate_stride(count)},
        buffer_(stride_ * num_arrays) {
    for (auto i{0uz}; i < num_arrays; ++i) {
      if (defaults[i] != T{0}) {
        std::fill_n((*this)[i], count_, defaults[i]);
      }
    }
  }

  [[nodiscard]] auto constexpr count() const noexcept -> std::size_t {
    return count_;
  }

  [[nodiscard]] auto constexpr stride() const noexcept -> std::size_t {
    return stride_;
  }

  [[nodiscard]]
  auto constexpr operator[](std::size_t array) noexcept -> T* {
    T* ptr{buffer_.data() + (count_ * array)};
    return std::assume_aligned<alignment>(ptr);
  }

  [[nodiscard]]
  auto constexpr operator[](std::size_t array) const noexcept -> T* {
    T* ptr{buffer_.data() + (count_ * array)};
    return std::assume_aligned<alignment>(ptr);
  }

  template <typename Enum>
    requires std::is_enum_v<Enum>
  [[nodiscard]] constexpr auto operator[](Enum array) const noexcept -> T* {
    return (*this)[std::to_underlying(array)];
  }
};

} // namespace em::mem
