#pragma once

#include <cstdint>

#ifndef SIMD_BYTES
#  if defined(__AVX512F__)
#    define SIMD_BYTES 64uz
#  elif defined(__AVX2__) || defined(__AVX__)
#    define SIMD_BYTES 32uz
#  else
#    define SIMD_BYTES 16uz
#  endif
#endif

namespace em::mem {

inline constexpr auto simd_bytes{SIMD_BYTES};
inline constexpr auto alignment_bytes{simd_bytes};

template <typename T>
inline constexpr auto default_align{
    (alignment_bytes > alignof(T)) ? alignment_bytes : alignof(T)};

template <typename T>
inline constexpr auto is_padded{sizeof(T) < alignment_bytes};

template <typename T, std::size_t alignment = default_align<T>>
constexpr auto handle_pad(std::size_t unpadded) noexcept -> std::size_t {
  if constexpr (is_padded<T>) {
    constexpr auto lanes{mem::alignment_bytes / sizeof(T)};
    return (unpadded + lanes - 1uz) & ~(lanes - 1uz);
  } else {
    return unpadded;
  }
}

} // namespace em::mem
