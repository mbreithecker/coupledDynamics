#include <iostream>
#include <chrono>
#include "InputOutput.h"
#include "Lattice2D.cuh"

/**
 * Allocate all needed arrays on the device
 * @param parameters Parameters struct defined in global.h containing all needed parameters for the whole simulation
 * @param lattice Lattice struct defined in global.h containing the device pointers for the lattices
 */
void prepareLatticeForDevice(Parameters &parameters, Lattice &lattice) {

    //Allocate memory for lattice on device:
    const size_t bytes = sizeof(TFloat) * parameters.volume;

    //Lattices
    cudaMalloc(&lattice.d_lattice_sigma, bytes);
    cudaMalloc(&lattice.d_lattice_pi, bytes);
    cudaMalloc(&lattice.d_lattice_n, bytes);
    cudaMalloc(&lattice.d_lattice_nu, bytes);

    cudaMalloc(&lattice.d_next_lattice_sigma, bytes);
    cudaMalloc(&lattice.d_next_lattice_pi, bytes);
    cudaMalloc(&lattice.d_next_lattice_n, bytes);
    cudaMalloc(&lattice.d_next_lattice_nu, bytes);

    //Random States
    cudaMalloc(&lattice.d_curandStateArray, sizeof(curandStateType) * parameters.volume);

    //vector random noise field
    cudaMalloc(&lattice.d_random_noise, bytes * parameters.dimension);
    cudaMalloc(&lattice.d_next_random_noise, bytes * parameters.dimension);

    //Observables arrays
    cudaMalloc(&lattice.d_obs_sigma_avg, bytes);
    cudaMalloc(&lattice.d_obs_sigma_sq_avg, bytes);
    cudaMalloc(&lattice.d_obs_pi_avg, bytes);
    cudaMalloc(&lattice.d_obs_pi_sq_avg, bytes);
    cudaMalloc(&lattice.d_obs_n_avg, bytes);
    cudaMalloc(&lattice.d_obs_n_sq_avg, bytes);
    cudaMalloc(&lattice.d_obs_nu_avg, bytes);
    cudaMalloc(&lattice.d_obs_nu_sq_avg, bytes);
}

/**
 * Release all allocated arrays from device.
 * @param lattice Device pointers for the lattices
 */
void finishLattice(Lattice &lattice) {
    cudaFree(lattice.d_lattice_sigma);
    cudaFree(lattice.d_lattice_pi);
    cudaFree(lattice.d_lattice_n);
    cudaFree(lattice.d_lattice_nu);

    cudaFree(lattice.d_next_lattice_sigma);
    cudaFree(lattice.d_next_lattice_pi);
    cudaFree(lattice.d_next_lattice_n);
    cudaFree(lattice.d_next_lattice_nu);

    cudaFree(lattice.d_curandStateArray);

    cudaFree(lattice.d_random_noise);
    cudaFree(lattice.d_next_random_noise);

    cudaFree(lattice.d_obs_sigma_avg);
    cudaFree(lattice.d_obs_sigma_sq_avg);
    cudaFree(lattice.d_obs_pi_avg);
    cudaFree(lattice.d_obs_pi_sq_avg);
    cudaFree(lattice.d_obs_n_avg);
    cudaFree(lattice.d_obs_n_sq_avg);
    cudaFree(lattice.d_obs_nu_avg);
    cudaFree(lattice.d_obs_nu_sq_avg);
}


int main(int argv, char **args) {

    std::cout << "============== Coupled Model A/B - Lattice Simulation ==============" << std::endl;
    std::cout << "                                               Version 1.1.5-RELEASE" << std::endl;

    //Read parameters from command line
    auto parameters = getParametersFromProgramOptions(args, argv);

    //Stats
    auto start_time = std::chrono::high_resolution_clock::now();
    printParameters(parameters);

    //-----Simulation------

    //Allocate Cuda Memory
    Lattice lattice;
    prepareLatticeForDevice(parameters, lattice);

    if (parameters.dimension == 2) {

        std::cout << "\n=========================== Observables ============================" << std::endl;

        std::vector<Observables> observablesList;

        for (int i = 0; i < parameters.number_of_simulations; i++) {

            //Update seed for next round
            parameters.seed += 3*i + i*i*i  + 2*parameters.volume + parameters.measure_time;

            //Run simulation and wait for the measurement
            Observables observables = Lattice2D::startSimulation(parameters, lattice);
            printObservables(observables);

            observablesList.push_back(observables);
        }

        std::cout << "====================================================================\n" << std::endl;

        writeToDatabase(parameters, observablesList);

        finishLattice(lattice);
    }else {
        std::cerr << "Dimension not supported" << std::endl;
        return 1;
    }

    auto milliSeconds = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start_time).count();
    printRuntime(milliSeconds);

    std::cout << "\n====================================================================\n" << std::endl;

    return 0;
}
