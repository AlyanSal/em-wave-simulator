#include "Output.hpp"
#include "Simulation.hpp"

int main() {
  const auto FRQ{2.4e9f};
  const auto WIDTH{2.5f};
  const auto LENGTH{2.5f};

  em::sim::simulation simulation(FRQ, WIDTH, LENGTH);
  Output writer("output.csv");

  simulation.setup_simulation();

  for (auto i{0uz}; i < 1500; ++i) {
    simulation.step_simulation();
    writer.write_frame(i, simulation.getE(), simulation.getH());
  }

  return 0;
}