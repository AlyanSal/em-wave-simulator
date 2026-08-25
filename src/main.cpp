#include "Simulation.hpp"

int main() {
  const auto FRQ{2.4e9f};
  const auto WIDTH{2.5f};
  const auto LENGTH{2.5f};

  em::simulation::simulation simulation(FRQ, WIDTH, LENGTH);

  return 0;
}