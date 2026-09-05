#include "Output.hpp"
#include "Simulation.hpp"
#include "Timer.hpp"

auto main() -> int {
  const auto FRQ{2.4e9f};
  const auto WIDTH{2.5f};
  const auto LENGTH{2.5f};
  const auto LARGEST_EPSILON{4.0f};

  std::function<float(float, float)> perm_func{
      [WIDTH](const float pos_x, const float pos_y) -> float { return 1.0f; }};

  std::function<float(float, float)> cond_func{
      [WIDTH](const float pos_x, const float pos_y) -> float { return 0.0f; }};

  std::function<float(float, float)> chi1_func{
      [WIDTH](const float pos_x, const float pos_y) -> float { return 0.0f; }};

  std::function<float(float, float)> t0_func{
      [WIDTH](const float pos_x, const float pos_y) -> float { return 1.0f; }};

  em::sim::simulation simulation(FRQ, WIDTH, LENGTH, LARGEST_EPSILON);
  Output writer("output.csv");

  auto setup{[&]() -> void {
    simulation.setup_simulation(perm_func, cond_func, chi1_func, t0_func);
  }};

  simulation.add_source(em::source::sin1D, 0.2f);

  em::util::timer(setup);

  auto run{[&]() -> void {
    for (auto i{0uz}; i < 1500; ++i) {
      simulation.step_simulation();
      writer.write_frame(i, simulation.getE(), simulation.getH());
    }
  }};

  em::util::timer(run);

  return 0;
}