#pragma once

#include <cstdlib>
#include <memory>

namespace em::mem {

template <typename T> struct BufferDeleter {
  auto operator()(T* ptr) const noexcept -> void { std::free(ptr); }
};

template <typename T, std::size_t alignment> class Buffer {
private:
  std::size_t count_;
  std::unique_ptr<T[], BufferDeleter<T>> data_;

  [[nodiscard]] static constexpr auto
  calculate_bytes(const std::size_t count) noexcept -> std::size_t {
    const auto raw_bytes{count * sizeof(T)};
    return (raw_bytes + alignment - 1uz) & ~(alignment - 1uz);
  }

public:
  explicit Buffer(std::size_t count)
      : count_{count},
        data_{static_cast<T*>(
            std::aligned_alloc(alignment, calculate_bytes(count)))} {}

  Buffer(Buffer&&) noexcept = default;
  auto operator=(Buffer&&) noexcept -> Buffer& = default;

  Buffer(const Buffer&) = delete;
  auto operator=(const Buffer&) -> Buffer& = delete;

  [[nodiscard]]
  constexpr auto count() const noexcept -> std::size_t {
    return count_;
  }

  [[nodiscard]]
  constexpr auto data() noexcept -> T* {
    return data_.get();
  }

  [[nodiscard]]
  constexpr auto data() const noexcept -> T* {
    return data_.get();
  }
};

} // namespace em::mem
