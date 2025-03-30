#include "system.h"
#include <iostream>
#include <fstream> // to create files, write to files, and read data from files



int main() {
    int num_atoms =100; //number of atoms spinning in 1D ising model
    double temperature = 1.0; // Temperature
    int mc_steps = 100000;  // Number of Monte Carlo steps
    int N =100; // 100x100 grids

    double betaMin =0.1;
    double betaMax =1.5;
    double betaStep =0.1;

    std::ofstream dataFile("ising_results.csv"); // Output file for results
    dataFile << "Beta,Average Energy,Average Magnetization\n"; // CSV header


    for (double beta =betaMin; beta<= betaMax; beta += betaStep){
        // Create Ising system object
        IsingSystem ising(num_atoms, beta);

        // Run Monte Carlo Simulation
        ising.simulate(mc_steps);

        // Compute energy and magnetization
        double avg_energy = ising.compute_energy() ;
        double avg_magnetization = abs(ising.compute_magnetization()); // Take absolute value

        // Print results
        ising.print_results(avg_energy,avg_magnetization);

        // Write to CSV file
        dataFile << beta << "," << avg_energy << "," << avg_magnetization << "\n";

    }

    dataFile.close(); // Close the output file
    std::cout << "Simulation results saved to 'ising_results.csv'" << std::endl;

    return 0;
}


// to compile and run the program -> g++ main.cpp system.cpp -o result
//                                  ./result

// to run python file (for visualisation) -> python plot_ising.py 