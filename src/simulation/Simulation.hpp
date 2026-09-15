#pragma once

#include <cmath>
#include <functional>

#include "Constants.hpp"
#include "Grid.hpp"
#include "Kernels.hpp"
#include "Sources.hpp"

namespace em::sim {

class simulation {
public:
  simulation(float frequency, float width, float length, float largest_eps,
             float pml_factor);

  void
  setup_simulation(std::function<float(float, float)>& permittivity_function,
                   std::function<float(float, float)>& conductivity_function,
                   std::function<float(float, float)>& chi1_function,
                   std::function<float(float, float)>& t0_function);

  void step_simulation() noexcept;

  void add_source(std::function<float(float)> source, float pos_x,
                  float pos_y) noexcept;

  void handle_sources() const noexcept;

  [[nodiscard]] auto constexpr getEmData() const noexcept
      -> std::array<float*, 3> {
    return {grid_[memory::MainField::EzField],
            grid_[memory::MainField::HxField],
            grid_[memory::MainField::HyField]};
  }

  [[nodiscard]] auto constexpr Dims() const noexcept
      -> std::array<std::size_t, 6> { // NOLINT
    return {rows_,          cols_,        interior_rows_,
            interior_cols_, pml_cells_x_, pml_cells_y_};
  }

  [[nodiscard]] auto PmlFactor() const noexcept -> float { return pml_factor_; }

  void setEpsilonDist();
  void setConductanceDist();

private:
  const float width_;
  const float length_;
  float largest_eps_;
  const float frequency_;
  const float smallest_wavelength_;
  const float pml_factor_;

  const float dx_;
  const float dt_;

  const std::size_t interior_rows_;
  const std::size_t interior_cols_;
  const std::size_t pml_cells_x_;
  const std::size_t pml_cells_y_;

  const std::size_t rows_;
  const std::size_t cols_;
  const std::size_t cells_;

  std::vector<std::pair<std::function<float(float)>, std::size_t>> sources_;

  std::size_t timestep_{};

  memory::grid grid_;
};

} // namespace em::sim
