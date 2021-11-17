//
// Created by maximilian on 30.08.21.
//

#ifndef MODELAB_CUDA_GLOBAL_H
#define MODELAB_CUDA_GLOBAL_H

#include <string>

typedef float TFloat;

#define POWER_TWO(a) ((a)*(a))
#define POWER_THREE(a) ((a)*(a)*(a))

#define DYNAMIC_CHANGE_MODE_RE_INIT 0
#define DYNAMIC_CHANGE_MODE_AVG_SUBTRACT 1
#define DYNAMIC_CHANGE_MODE_ZERO 2

#define curandStateType curandStateMRG32k3a

//Cuda threads for each block.
const int NUM_THREADS = 128;

/**
 * Lattice - struct
 * Contains all pointers (DEVCIE) for the field-arrays
 */
typedef struct {

    //Main fields
    TFloat *d_lattice_sigma;
    TFloat *d_lattice_pi;
    TFloat *d_lattice_n;
    TFloat *d_lattice_nu;

    //Buffer field for leapfrog integration
    TFloat *d_next_lattice_sigma;
    TFloat *d_next_lattice_pi;
    TFloat *d_next_lattice_n;
    TFloat *d_next_lattice_nu;

    //Random states
    void * d_curandStateArray; //type: curandStateMRG32k3a_t, can't use this type because .h-file is no .cuh file

    //vector random noise
    TFloat * d_random_noise; // [ dim1,dim2,(dim3)]
    TFloat * d_next_random_noise;

    //Observables arrays
    TFloat * d_obs_sigma_avg;
    TFloat * d_obs_sigma_sq_avg;
    TFloat * d_obs_pi_avg;
    TFloat * d_obs_pi_sq_avg;
    TFloat * d_obs_n_avg;
    TFloat * d_obs_n_sq_avg;
    TFloat * d_obs_nu_avg;
    TFloat * d_obs_nu_sq_avg;
} Lattice;

/**
 * Parameters - struct
 * Contains physical constants, program Options and calculated values
 * Is passed by copy-by-value
 */
typedef struct {
    //lattice props
    unsigned int dimension;
    unsigned int sites;
    unsigned int volume; // = sites^dimension

    //LatticeEquation
    TFloat msq_A;
    TFloat lambda_A;
    TFloat gamma_A;
    TFloat msq_B;
    TFloat gamma_B;
    TFloat J;
    TFloat C;
    TFloat mu;
    TFloat temperature;

    //Time Integration / Simulation
    TFloat integrator_timeDelta;
    bool fastThermalization;
    unsigned int thermalization_time_1;
    unsigned int thermalization_time_2;
    unsigned int thermalization_time_3;
    unsigned int measure_time;
    unsigned long long seed;
    unsigned int number_of_simulations;
    unsigned int dynamicChangeMode;

    //Cached variables to avoid sqrt computation in kernel
    TFloat integrator_timeDelta_SQRT;
    TFloat standard_deviation_B;
    TFloat standard_deviation_A;

    //External field for thermalization
    TFloat symmetryBreakingJ;

    //InputOutput vars/flags
    bool flag_printEveryStep;
    std::string outFileName;
} Parameters;


/**
 * Observalbes - struct
 * Used on the the host-system. The calculated observables form the lattices are stored here
 */
typedef struct {
    TFloat sigma_avg;
    TFloat sigma_sq_avg;
    TFloat pi_avg;
    TFloat pi_sq_avg;

    TFloat n_avg;
    TFloat n_sq_avg;
    TFloat nu_avg;
    TFloat nu_sq_avg;
} Observables;


#endif //MODELAB_CUDA_GLOBAL_H
