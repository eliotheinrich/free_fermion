#pragma once

#include <numbers>
#include <FreeFermion.h>

class AdaptiveFermionSimulator : public FreeFermionSimulator {
  private:
    double p;
    double r;

    std::vector<bool> active;

  public:
    AdaptiveFermionSimulator(dataframe::Params& params, uint32_t num_threads) : FreeFermionSimulator(params, num_threads) {
      p = dataframe::utils::get<double>(params, "p");
      r = dataframe::utils::get<double>(params, "r");

      active = std::vector<bool>(L, true);
      init_fermion_state(false);
      state->checkerboard_particles();
    }

    Eigen::MatrixXcd hamiltonian(size_t i, size_t j) const {
      Eigen::MatrixXcd H = Eigen::MatrixXcd::Zero(L, L);
      H(i, j) = std::numbers::pi/4;
      H(j, i) = std::numbers::pi/4;

      return H;
    }

    void unitary_timestep(size_t i, size_t j) {
      if (active[i] || active[j]) {
        auto H = hamiltonian(i, j);
        auto A= Eigen::MatrixXcd::Zero(L, L);

        if (rand() % 2) {
          state->evolve(H, A);
        } else {
          state->evolve(A, H);
        }

        active[i] = true;
        active[j] = true;
      }
    }

    void adaptive_timestep(size_t i, size_t j, bool b) {
      bool b1 = state->projective_measurement(i, randf());
      bool b2 = state->projective_measurement(j, randf());
      std::cout << "calling adaptive timestep\n";
      if (b1 == b2) {
        return;
      }

      std::cout << fmt::format("b1, b2, b = {}, {}, {}\n", b1, b2, b);

      if (b1 == b) {
        active[i] = false;
        active[j] = false;
      } else {
        if (randf() < r) {
          active[i] = false;
          active[j] = false;
          std::cout << "calling swap\n";
          state->swap(i, j);
        }
      }
    }

    virtual void timesteps(uint32_t num_steps) override {
      std::cout << "STARTING AMPLITUDES = \n" << state->amplitudes << "\n";
      for (size_t j = 0; j < L/2; j++)  {
        unitary_timestep(2*j, 2*j + 1);
      }

      for (size_t j = 0; j < L/2; j++) {
        adaptive_timestep(2*j, 2*j + 1, true);
      }

      for (size_t j = 0; j < L/2; j++) {
        unitary_timestep(2*j + 1, (2*j + 2) % L);
      }

      for (size_t j = 0; j < L/2; j++) {
        adaptive_timestep(2*j + 1, (2*j + 2) % L, false);
      }
    }
};


