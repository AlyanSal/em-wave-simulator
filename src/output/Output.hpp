#pragma once

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

  void write_frame(const std::size_t timestep, const std::size_t Nx,
                   const std::size_t Ny, const std::vector<float>& Ez,
                   const std::vector<float>& Hx, const std::vector<float>& Hy) {
    for (auto j{0uz}; j < Ny; ++j) {
      for (auto i{0uz}; i < Nx; ++i) {
        const std::size_t idx{(j * Nx) + i};
        file << timestep << ',' << i << ',' << j << ',' << Ez[idx] << ','
             << Hx[idx] << ',' << Hy[idx] << '\n';
      }
    }
  }
};