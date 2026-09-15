#pragma once

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <print>
#include <string>
#include <vector>

namespace em::output {

/**
 * Fast binary output writer for simulation frames.
 * Writes a compact 20-byte header followed by contiguous raw float32 memory
 * blocks for Ez, Hx, and Hy at each recorded timestep.
 */
class BinaryOutput {
public:
  // Magic identifier: 'E', 'M', 'W', '2'
  static constexpr uint32_t MAGIC{0x32574D45};
  static constexpr uint32_t VERSION{1};

  BinaryOutput(const std::string& filename, const std::size_t rows,
               const std::size_t cols, const uint32_t channels = 3)
      : rows_{static_cast<uint32_t>(rows)},
        cols_{static_cast<uint32_t>(cols)},
        channels_{channels} {
    file_.open(filename, std::ios::binary | std::ios::out | std::ios::trunc);
    if (!file_.is_open()) {
      std::print(std::cerr, "BinaryOutput:: could not open file: {}\n",
                 filename);
      return;
    }

    // 20-byte header: MAGIC, VERSION, ROWS, COLS, CHANNELS
    file_.write(reinterpret_cast<const char*>(&MAGIC), sizeof(MAGIC));
    file_.write(reinterpret_cast<const char*>(&VERSION), sizeof(VERSION));
    file_.write(reinterpret_cast<const char*>(&rows_), sizeof(rows_));
    file_.write(reinterpret_cast<const char*>(&cols_), sizeof(cols_));
    file_.write(reinterpret_cast<const char*>(&channels_), sizeof(channels_));
  }

  ~BinaryOutput() {
    if (file_.is_open()) {
      file_.flush();
      file_.close();
    }
  }

  BinaryOutput(const BinaryOutput&) = delete;
  auto operator=(const BinaryOutput&) -> BinaryOutput& = delete;
  BinaryOutput(BinaryOutput&&) noexcept = default;
  auto operator=(BinaryOutput&&) noexcept -> BinaryOutput& = default;

  /**
   * Writes one simulation frame containing timestep and raw float arrays for
   * Ez, Hx, Hy. Performs direct memory block transfers without per-cell loops
   * or string formatting.
   */
  void write_frame(const std::size_t timestep, const std::size_t size,
                   const float* Ez, const float* Hx, const float* Hy) {
    if (!file_.is_open()) {
      return;
    }

    const auto t_step{static_cast<uint32_t>(timestep)};
    file_.write(reinterpret_cast<const char*>(&t_step), sizeof(uint32_t));

    const auto e_size{size};
    const auto hx_size{size};
    const auto hy_size{size};

    file_.write(reinterpret_cast<const char*>(Ez),
                static_cast<std::streamsize>(e_size * sizeof(float)));

    if (channels_ >= 2) {
      file_.write(reinterpret_cast<const char*>(Hx),
                  static_cast<std::streamsize>(hx_size * sizeof(float)));
    }

    if (channels_ >= 3) {
      file_.write(reinterpret_cast<const char*>(Hy),
                  static_cast<std::streamsize>(hy_size * sizeof(float)));
    }
  }

  [[nodiscard]] auto is_open() const noexcept -> bool {
    return file_.is_open();
  }
  [[nodiscard]] auto rows() const noexcept -> uint32_t { return rows_; }
  [[nodiscard]] auto cols() const noexcept -> uint32_t { return cols_; }
  [[nodiscard]] auto channels() const noexcept -> uint32_t { return channels_; }

private:
  std::ofstream file_;
  uint32_t rows_{0};
  uint32_t cols_{0};
  uint32_t channels_{3};
};

} // namespace em::output
