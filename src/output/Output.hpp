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

    file << "timestep,cell,Ex,Hy" << '\n';
  }

  ~Output() { file.close(); };

  void write_frame(std::size_t timestep, std::vector<float>& Efield,
                   std::vector<float>& Hfield) {
    for (auto i{0uz}; i < Efield.size(); ++i) {
      file << timestep << ',' << i << ',' << Efield[i] << ',' << Hfield[i]
           << '\n';
    }
  }
};