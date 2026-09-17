#pragma once

#include <algorithm>
#include <condition_variable>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <mutex>
#include <print>
#include <string>
#include <thread>
#include <utility>

#include "Buffer.hpp"
#include "Memory.hpp"

namespace em::output {

class BinaryOutput {
public:
  static constexpr uint32_t MAGIC{0x32574D45};
  static constexpr uint32_t VERSION{1};

  BinaryOutput(const std::string& filename, const std::size_t rows,
               const std::size_t cols, const uint32_t channels = 3)
      : rows_{static_cast<uint32_t>(rows)},
        cols_{static_cast<uint32_t>(cols)},
        channels_{channels},
        channel_bytes_{rows_ * cols_ * sizeof(float)},
        frame_bytes_{sizeof(uint32_t) + (channels_ * channel_bytes_)},
        buf_produce_(frame_bytes_),
        buf_write_(frame_bytes_) {
    file_.open(filename, std::ios::binary | std::ios::out | std::ios::trunc);
    if (!file_.is_open()) {
      std::print(std::cerr, "BinaryOutput:: could not open file: {}\n",
                 filename);
      return;
    }

    // 20-byte header
    file_.write(reinterpret_cast<const char*>(&MAGIC), sizeof(MAGIC));
    file_.write(reinterpret_cast<const char*>(&VERSION), sizeof(VERSION));
    file_.write(reinterpret_cast<const char*>(&rows_), sizeof(rows_));
    file_.write(reinterpret_cast<const char*>(&cols_), sizeof(cols_));
    file_.write(reinterpret_cast<const char*>(&channels_), sizeof(channels_));

    worker_ = std::jthread([this]() -> void { writer_loop(); });
  }

  ~BinaryOutput() {
    {
      std::scoped_lock<std::mutex> lock(mutex_);
      done_ = true;
    }
    cv_produce_.notify_one();

    if (worker_.joinable())
      worker_.join();

    if (file_.is_open()) {
      file_.flush();
      file_.close();
    }
  }

  BinaryOutput(const BinaryOutput&) = delete;
  auto operator=(const BinaryOutput&) -> BinaryOutput& = delete;
  BinaryOutput(BinaryOutput&&) noexcept = delete;
  auto operator=(BinaryOutput&&) noexcept -> BinaryOutput& = delete;

  void write_frame(const std::size_t timestep, const std::size_t size,
                   const float* Ez, const float* Hx, const float* Hy) {
    if (!file_.is_open())
      return;

    {
      std::unique_lock<std::mutex> lock(mutex_);
      cv_consume_.wait(lock, [this] -> bool { return !has_data_; });

      char* dst = buf_produce_.data();
      const auto t_step{static_cast<uint32_t>(timestep)};
      std::memcpy(dst, &t_step, sizeof(uint32_t));
      dst += sizeof(uint32_t);

      const std::size_t bytes{size * sizeof(float)};
      if (Ez) {
        std::memcpy(dst, Ez, bytes);
        dst += bytes;
      }
      if (channels_ >= 2 && Hx) {
        std::memcpy(dst, Hx, bytes);
        dst += bytes;
      }
      if (channels_ >= 3 && Hy) {
        std::memcpy(dst, Hy, bytes);
        dst += bytes;
      }

      std::swap(buf_produce_, buf_write_);
      has_data_ = true;
    }

    cv_produce_.notify_one();
  }

  [[nodiscard]] auto is_open() const noexcept -> bool {
    return file_.is_open();
  }
  [[nodiscard]] auto rows() const noexcept -> uint32_t { return rows_; }
  [[nodiscard]] auto cols() const noexcept -> uint32_t { return cols_; }
  [[nodiscard]] auto channels() const noexcept -> uint32_t { return channels_; }

private:
  void writer_loop() {
    while (true) {
      std::unique_lock<std::mutex> lock(mutex_);
      cv_produce_.wait(lock, [this] -> bool { return has_data_ || done_; });

      if (done_ && !has_data_) {
        break;
      }

      lock.unlock();
      file_.write(buf_write_.data(),
                  static_cast<std::streamsize>(frame_bytes_));

      lock.lock();
      has_data_ = false;
      cv_consume_.notify_one();
    }
  }

  std::ofstream file_;
  uint32_t rows_{0};
  uint32_t cols_{0};
  uint32_t channels_{3};
  std::size_t channel_bytes_{0};
  std::size_t frame_bytes_{0};

  // Double Vuffer
  mem::Buffer<char, mem::default_align<char>> buf_produce_;
  mem::Buffer<char, mem::default_align<char>> buf_write_;

  // Thread Synchro
  std::mutex mutex_;
  std::condition_variable cv_produce_;
  std::condition_variable cv_consume_;
  bool has_data_{false};
  bool done_{false};
  std::jthread worker_;
};

} // namespace em::output
