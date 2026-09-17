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

public:
  explicit Buffer(std::size_t count)
      : count_{count},
        data_{std::make_unique<T*>(std::aligned_alloc(alignment, count))} {}

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
