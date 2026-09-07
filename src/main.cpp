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

  std::function<float(float, float)> perm_func{
      []([[maybe_unused]] const float pos_x,
         [[maybe_unused]] const float pos_y) -> float { return 1.0f; }};

  std::function<float(float, float)> cond_func{
      []([[maybe_unused]] const float pos_x,
         [[maybe_unused]] const float pos_y) -> float { return 0.0f; }};

  std::function<float(float, float)> chi1_func{
      []([[maybe_unused]] const float pos_x,
         [[maybe_unused]] const float pos_y) -> float { return 0.0f; }};

  std::function<float(float, float)> t0_func{
      []([[maybe_unused]] const float pos_x,
         [[maybe_unused]] const float pos_y) -> float { return 1.0f; }};

  em::sim::simulation simulation(FRQ, WIDTH, LENGTH, LARGEST_EPSILON,
                                 PML_FACTOR);
  BinaryOutput writer("output.bin", simulation.Rows(), simulation.Cols());

  auto setup{[&]() -> void {
    simulation.setup_simulation(perm_func, cond_func, chi1_func, t0_func);
  }};

  simulation.add_source(em::source::sin1D, LENGTH / 4, WIDTH / 2);
  simulation.add_source(em::source::cos1D, LENGTH / 4 * 3, WIDTH / 2);

  em::util::timer(setup);

  auto run{[&]() -> void {
    for (auto i{0uz}; i < 1500; ++i) {
      simulation.step_simulation();
      if (i % 5 == 0)
        writer.write_frame(i, simulation.getEz(), simulation.getHx(),
                           simulation.getHy());
    }
  }};

  em::util::timer(run);

  return 0;
}