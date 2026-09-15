#include "BinaryOutput.hpp"
#include "Output.hpp"
#include "Simulation.hpp"
#include "Timer.hpp"

auto main() -> int {
  const auto FRQ{1e9f};
  const auto WIDTH{2.5f};
  const auto LENGTH{3.5f};
  const auto LARGEST_EPSILON{1.0f};
  const auto PML_FACTOR{0.2f};
  const auto SLIT_WIDTH{0.45f};
  const auto WALL_X{1.0f};
  const auto SIM_TIME{2000uz};

  std::function<float(float, float)> perm_func{
      []([[maybe_unused]] const float pos_x,
         [[maybe_unused]] const float pos_y) -> float { return 1.0f; }};

  std::function<float(float, float)> cond_func{
      [WIDTH, WALL_X, SLIT_WIDTH]([[maybe_unused]] const float pos_x,
                                  [[maybe_unused]] const float pos_y) -> float {
        if (pos_x >= WALL_X && pos_x <= WALL_X + 0.05f) {
          if (std::abs(pos_y - WIDTH / 2.0f) > (SLIT_WIDTH / 2.0f)) {
            return 1e4f;
          }
        }
        return 0.0f;
      }};

  std::function<float(float, float)> chi1_func{
      []([[maybe_unused]] const float pos_x,
         [[maybe_unused]] const float pos_y) -> float { return 0.0f; }};

  std::function<float(float, float)> t0_func{
      []([[maybe_unused]] const float pos_x,
         [[maybe_unused]] const float pos_y) -> float { return 1.0f; }};

  em::sim::simulation simulation(FRQ, WIDTH, LENGTH, LARGEST_EPSILON,
                                 PML_FACTOR);

  const auto [ROWS, COLS, _, _, _, _] = simulation.Dims();
  em::output::BinaryOutput writer("output.bin", ROWS, COLS);

  auto setup{[&]() -> void {
    simulation.setup_simulation(perm_func, cond_func, chi1_func, t0_func);
  }};

  for (float y{}; y < WIDTH; y += 0.1f) {
    simulation.add_source(em::source::sin1D, LENGTH / 8, y);
  }

  em::util::timer(setup);

  const auto [EZ, HX, HY] = simulation.getEmData();

  auto run{[&]() -> void {
    for (auto i{0uz}; i < SIM_TIME; ++i) {
      simulation.step_simulation();
      if (i % 5 == 0)
        writer.write_frame(i, ROWS * COLS, EZ, HX, HY);
    }
  }};

  em::util::timer(run);

  return 0;
}