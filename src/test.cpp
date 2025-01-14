#include "ChainSimulator.hpp"
#include "AdaptiveFermionSimulator.hpp"
#include <iostream>
#include <Frame.h>
#include <TimeSamplingDriver.hpp>

using namespace dataframe;

void test_simulator() {
  size_t L = 2;

  dataframe::Params p;
  p["seed"] = 31.0;
  p["p"] = 0.5;
  p["r"] = 0.5;
  p["system_size"] = static_cast<double>(L);
  p["sample_correlations"] = 1.0;
  p["renyi_indices"] = "1, 2";


  p["equilibration_timesteps"] = 10000.0;
  p["sampling_timesteps"] = 10000.0;
  p["sample_entropy"] = 1.0;

  AdaptiveFermionSimulator simulator(p, 1);
  Eigen::MatrixXcd amplitudes(4, 2);
  amplitudes << 1, 0,
                0, 0,
                0, 0.707,
                0, -0.707;

  simulator.state->amplitudes = amplitudes;
  simulator.state->projective_measurement(1, simulator.randf());
}

int main() {
  test_simulator();
  throw std::runtime_error("FINISHED");

  size_t L = 3;

  dataframe::Params p;
  p["seed"] = 31.0;
  p["p"] = 0.5;
  p["r"] = 0.5;
  p["system_size"] = static_cast<double>(L);
  p["sample_correlations"] = 1.0;
  p["renyi_indices"] = "1, 2";


  p["equilibration_timesteps"] = 10000.0;
  p["sampling_timesteps"] = 10000.0;
  p["sample_entropy"] = 1.0;

  //std::unique_ptr<ChainSimulator> f = std::make_unique<ChainSimulator>(p, 1);
  //SimulatorDisplay display(std::move(f), 1, 30);
  //display.animate();

  TimeSamplingDriver<AdaptiveFermionSimulator> driver(p);
  driver.init_simulator(1);
  auto slide = driver.generate_dataslide();
  //std::cout << slide.to_json() << std::endl;
}
