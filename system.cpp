// This is the cpp file -> for implementation
// to implement the ising model
#include "system.h"
#include <fstream> // Include for file output
#include <iostream>
#include <cstdlib>  // For rand() -> to return random number between 0 and 1
#include <cmath>    // For exp()
#include <ctime>    // For seeding random numbers

// Constructor: Initializes system size and temperature
IsingSystem::IsingSystem(int num_atoms, double beta) {
    N = num_atoms;
    this->beta =beta; // to store beta
    //spins.resize(N); // this is for 1d

      for (int i = 0; i < N; i++) {
        spins[i].resize(N);  
    }

    initialize_spins();
}

// Randomly initialize spins as +1 or -1
void IsingSystem::initialize_spins() {
    srand(time(0)); // Seed random number generator
    for (int i = 0; i < N; i++) {
        for (int j=0; j< N; j++){
            spins[i][j] = (rand() % 2 == 0) ? 1 : -1; // If rand() % 2 == 0 is true, assigns 1 to spins[i].
                                              // Otherwise, assigns -1 to spins[i].
        }
    }
}

// Compute total energy of system (the change of energy)
double IsingSystem::compute_energy() const {
    double energy = 0.0; // stores the total energy of the system
    for (int i = 0; i < N - 1; i++) {
        for (int j = 0; j < N - 1; j++){
            int spin = spins[i][j];

            // Periodic boundary conditions (PBC)
            int right = spins[i][(j + 1) % N];
            int down  = spins[(i + 1) % N][j];

            // Interaction energy with right and bottom neighbors
            energy += -spin * (right + down);

        }
    }
    return energy / (N*N); // normalising by the number of spins
}

// Compute total magnetization
double IsingSystem::compute_magnetization() const {
    double magnetization = 0.0; // stores the total magnetisation of the system
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++){
            magnetization += spins[i][j];
        }
    }
    return magnetization/ (N*N); // normalising by the number of spins;
}

// Perform Monte Carlo update using Metropolis Algorithm -> Probability for a system being in energy state Ek -> (P(Ek))
void IsingSystem::monte_carlo_step() {
    for (int i = 0; i < N; i++) { // Loop N*N times for a full Monte Carlo sweep

        int row = rand() % N; // Pick a random row
        int column = rand() % N; // Pick a random column

        int spin_old = spins[row][column];
        int spin_new = -spin_old; // Flip spin

        // Compute energy difference ΔE = 2 * spin_old * sum(neighbors)
        
        // This is for 1D
        // double deltaE = 0;
        // if (index > 0) deltaE += 2 * spin_old * spins[index - 1]; // Left neighbor
        // if (index < N - 1) deltaE += 2 * spin_old * spins[index + 1]; // Right neighbor

        //This is for 2D
        int deltaE = 2 * spin_old * (
            spins[(row + 1) % N][column] + spins[(row - 1 + N) % N][column] + // Top & Bottom
            spins[row][(column + 1) % N] + spins[row][(column - 1 + N) % N]   // Left & Right
        );

        // Metropolis acceptance criterion
        if (deltaE < 0 || exp(-deltaE * beta) > (double)rand() / RAND_MAX) {
            spins[row][column] = spin_new; // Accept flip
        }
    }
}

void IsingSystem::save_spins(const std::string& filename) {
    std::ofstream spinFile(filename, std::ios::app); // open file in append mode

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            spinFile << spins[i][j] << " ";
        }
        spinFile << "\n"; // new line for each row
    }
    spinFile << "END\n"; // mark the end of configuration
    spinFile.close();
}



// To Run Monte Carlo Simulation
void IsingSystem::simulate(int mc_steps){
    std::ofstream spinFile("ising_spins.txt", std::ios::trunc); // to create file
    //spinFile.close(); // Close file

    for(int step=0; step<mc_steps; step++){
        monte_carlo_step(); // performing MC update

        // Save lattice configuration every 100 steps
        if (step % 100 == 0) {
            // std::ofstream spinFile("ising_spins.txt", std::ios::app); // Append mode
            // if (!spinFile) {
            //     std::cerr << "Error: Could not open ising_spins.txt" << std::endl;
            //     return;
            // }

            // for (int i = 0; i < N; i++) {
            //     for (int j = 0; j < N; j++) {
            //         spinFile << spins[i][j] << " ";
            //     }
            //     spinFile << "\n"; // New row
            // }
            // spinFile << "END\n"; // Mark end of one configuration
            std::cout << "Step " << step << ": Saving spin configuration" << std::endl; // Debug print
            save_spins("ising_spins.txt");
        }
    }
    spinFile.close(); // Close file
}

// Print final energy and magnetization
void IsingSystem::print_results(double avg_energy, double avg_magnetization) const{

    std::cout << "Beta: " << beta 
              << " | Final Energy: " << avg_energy 
              << " | Final Magnetization: " << avg_magnetization << std::endl;
}

