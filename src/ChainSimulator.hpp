#pragma once

#include <FreeFermion.h>

#define CS_SIMPLE_COUPLING 0
#define CS_RANDOM_COUPLING 1

class ChainSimulator : public FreeFermionSimulator {
  private:
    double p1;
    double p2;
    double beta;

    bool do_measurements;

    int coupling_type;

  public:
    ChainSimulator(dataframe::Params& params, uint32_t num_threads) : FreeFermionSimulator(params, num_threads) {
      p1 = dataframe::utils::get<double>(params, "p1", 0.0);
      p2 = dataframe::utils::get<double>(params, "p2", 0.0);

      beta = dataframe::utils::get<double>(params, "beta", 1.0);

      coupling_type = dataframe::utils::get<int>(params, "coupling_type");

      do_measurements = dataframe::utils::get<int>(params, "do_measurements", 0);

      init_fermion_state(true);

      std::string initial_state = dataframe::utils::get<std::string>(params, "initial_state");
      if (initial_state == "single_particle") {
        state->single_particle();
      } else if (initial_state == "all_particles") {
        state->all_particles();
      } else if (initial_state == "checkerboard") {
        state->checkerboard_particles();
      }
    }

    Eigen::MatrixXcd simple_coupling() const {
      Eigen::MatrixXcd H = Eigen::MatrixXcd::Zero(L, L);
      for (size_t i = 0; i < L; i++) {
        H(i, (i+1)%L) = 1.0;
        H((i+1)%L, i) = 1.0;
      }

      return H;
    }

    Eigen::MatrixXcd random_coupling() {
      Eigen::MatrixXcd H = Eigen::MatrixXcd::Zero(L, L);
      for (size_t i = 0; i < L; i++) {
        double K = 1.0;
        if (randf() < p1) {
          K = -K;
        }

        H(i, (i+1)%L) = K;
        H((i+1)%L, i) = K;
      }
      
      return H;
    }

    Eigen::MatrixXcd measurement_coupling() {
      Eigen::MatrixXcd H = Eigen::MatrixXcd::Zero(L, L);
      for (size_t i = 0; i < L; i++) {
        double K = 1.0;
        if (randf() < p2) {
          K = 0.0;
        }

        H(i, (i+1)%L) = K;
        H((i+1)%L, i) = K;
      }
      
      return H;
    }

    virtual void timesteps(uint32_t num_steps) override {
      Eigen::MatrixXcd H;
      if (coupling_type == CS_SIMPLE_COUPLING) {
        H = simple_coupling();
      } else if (coupling_type == CS_RANDOM_COUPLING) {
        H = random_coupling();
      }

      state->evolve(H);

      if (do_measurements) {
        H = measurement_coupling();
        state->weak_measurement(H, beta);
      }
    }
};


