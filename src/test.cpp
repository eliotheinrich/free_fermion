#include "ChainSimulator.hpp"
#include <iostream>

int main() {
  size_t L = 50;

  dataframe::Params p;
  p["dt"] = 0.2;
  p["L"] = static_cast<double>(L);
  p["initial_state"] = "checkerboard";
  p["do_measurements"] = 1.0;
  p["coupling_type"] = 1.0;
  p["p1"] = 0.5;
  p["p2"] = 0.7;

  std::unique_ptr<ChainSimulator> f = std::make_unique<ChainSimulator>(p, 1);
  SimulatorDisplay display(std::move(f), 1, 30);
  display.animate();
}
