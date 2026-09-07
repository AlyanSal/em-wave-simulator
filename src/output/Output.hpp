#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

class Output {
private:
  std::fstream file;

public:
  Output(const std::string& filename) {
    file = std::fstream(filename, std::ios::out);
    if (!file.is_open()) {
      std::print(std::cout, "EM:: could not open file: {}", filename);
    }

    file << "timestep,row,col,Ez,Hx,Hy" << '\n';
  }

  ~Output() { file.close(); };

  static inline auto float_to_int8(const float val, const float max_val) noexcept -> int8_t {
    if (max_val <= 0.0f) return 0;
    const float clamped = std::clamp(val / max_val, -1.0f, 1.0f);
    return static_cast<int8_t>(std::round(clamped * 127.0f));
  }

  void write_frame(const std::size_t timestep, const std::size_t Nx,
                   const std::size_t Ny, const std::vector<float>& Ez,
                   const std::vector<float>& Hx, const std::vector<float>& Hy) {
    for (auto j{0uz}; j < Ny; ++j) {
      for (auto i{0uz}; i < Nx; ++i) {
        const std::size_t idx{(j * Nx) + i};
        file << timestep << ',' << j << ',' << i << ',' << Ez[idx] << ','
             << Hx[idx] << ',' << Hy[idx] << '\n';
      }
    }
  }

  void write_byte_frame(const std::size_t timestep, const std::size_t Nx,
                        const std::size_t Ny, const std::vector<float>& Ez,
                        const std::vector<float>& Hx, const std::vector<float>& Hy,
                        const float e_max = 1.0f, const float h_max = 1.0f / 376.73f) {
    for (auto j{0uz}; j < Ny; ++j) {
      for (auto i{0uz}; i < Nx; ++i) {
        const std::size_t idx{(j * Nx) + i};
        file << timestep << ',' << j << ',' << i << ','
             << static_cast<int>(float_to_int8(Ez[idx], e_max)) << ','
             << static_cast<int>(float_to_int8(Hx[idx], h_max)) << ','
             << static_cast<int>(float_to_int8(Hy[idx], h_max)) << '\n';
      }
    }
  }
};