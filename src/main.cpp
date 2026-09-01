#include "Output.hpp"
#include "Simulation.hpp"
#include "Timer.hpp"

auto main() -> int {
  const auto FRQ{2.4e9f};
  const auto WIDTH{2.5f};
  const auto LENGTH{2.5f};
  const auto LARGEST_EPSILON{4.0f};

  std::function<float(float)> perm_func{[WIDTH](const float pos_x) -> float {
    return 1.0f + (3.0f * (pos_x > WIDTH / 2));
  }};

  std::function<float(float)> cond_func{[WIDTH](const float pos_x) -> float {
    return 0.04f * (pos_x > WIDTH / 2);
  }};

  em::sim::simulation simulation(FRQ, WIDTH, LENGTH, LARGEST_EPSILON);
  Output writer("output.csv");

  auto setup{
      [&]() -> void { simulation.setup_simulation(perm_func, cond_func); }};

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