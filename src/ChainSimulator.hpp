#pragma once

#include <FreeFermion.h>

#define CS_SIMPLE_COUPLING 0
#define CS_RANDOM_COUPLING 1

class ChainSimulator : public FreeFermionSimulator {
  private:
    double dt;
    double p1;
    double p2;

    bool do_measurements;

    int coupling_type;

  public:
    ChainSimulator(dataframe::Params& params, uint32_t num_threads) : FreeFermionSimulator(params) {
      dt = dataframe::utils::get<double>(params, "dt");
      p1 = dataframe::utils::get<double>(params, "p1", 0.0);
      p2 = dataframe::utils::get<double>(params, "p2", 0.0);

      coupling_type = dataframe::utils::get<int>(params, "coupling_type");

      do_measurements = dataframe::utils::get<int>(params, "do_measurement", 0);

      std::string initial_state = dataframe::utils::get<std::string>(params, "initial_state");
      if (initial_state == "single_particle") {
        single_particle();
      } else if (initial_state == "all_particles") {
        all_particles();
      } else if (initial_state == "checkerboard") {
        checkerboard_particles();
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

      std::cout << fmt::format("Orthogonal: {}\n", check_orthogonality());

      Eigen::MatrixXcd A = Eigen::MatrixXcd::Zero(L, L);

      evolve(A, H);
      std::cout << fmt::format("num particles: {}\n", num_particles());

      if (do_measurements) {
        H = measurement_coupling();
        weak_measurement(A, H);
      }
    }
};


