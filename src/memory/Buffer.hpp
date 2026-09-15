#pragma once

#include <memory>

namespace em::memory {

template <typename T> class Buffer {
private:
  std::size_t count_;
  std::unique_ptr<T[]> data_;

public:
  explicit Buffer(std::size_t count)
      : count_{count},
        data_{std::make_unique<T[]>(count)} {}

  [[nodiscard]]
  constexpr auto count() const noexcept -> std::size_t {
    return count_;
  }

  [[nodiscard]]
  constexpr auto data() noexcept -> T * {
    return data_.get();
  }

  [[nodiscard]]
  constexpr auto data() const noexcept -> T * {
    return data_.get();
  }
};

} // namespace em::memory
