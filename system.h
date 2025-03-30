#ifndef SYSTEM_H
#define SYSTEM_H

// This is the Header file -> for interface
//for the spins, energy calculations, and monte carlo updates

#include <vector>

class IsingSystem {
    private:
        int N; // grid size (NxN)
        double beta; // beta is the inverse temperature
        //std::vector<int> spins; // represents the spin states (either +1 or -1) for 1D
        std::vector<std::vector<int>> spins; // this for 2D grid that has 100 rows and each contains 100 spins

    
    public:
        IsingSystem(int num_atoms, double beta);

        void initialize_spins();    // Randomly initialize spins
        double compute_energy() const; // Compute system energy
        double compute_magnetization() const; // Compute magnetization
        void monte_carlo_step();    // Perform Monte Carlo update

        void simulate(int mc_steps); // Run Monte Carlo simulation
        void print_results(double avg_energy, double avg_magnetization) const; // Print final energy and magnetization

};

#endif