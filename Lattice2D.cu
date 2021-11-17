//
// Created by maximilian on 30.08.21.
//

#include <random>
#include "Lattice2D.cuh"
#include "InputOutput.h"
#include "global.h"

/**
 * Contains all kernels for lattice-calculation on a 2D Lattice
 * as well as one host method to perform kernel calls
 */
namespace Lattice2D {


    /**
     * Kernels and methods for Langevin-evolution (physical part)
     */
    namespace langevin {

        /**
         * Init lattice with random distribution. Set observables to zero.
         * @param lat
         * @param parameters
         */
        __global__ void device_initLattice(Lattice lat, Parameters parameters) {
            //Check array range
            unsigned int tid = (blockIdx.x * blockDim.x) + threadIdx.x;
            if(tid >= parameters.volume) return;

            //Init random states once
            curandStateType randomState{};
            curand_init(parameters.seed, tid, 0, &randomState);
            ((curandStateType *)lat.d_curandStateArray)[tid] = randomState;

            //Reset observables
            lat.d_obs_sigma_avg[tid] = 0;
            lat.d_obs_sigma_sq_avg[tid] = 0;
            lat.d_obs_n_avg[tid] = 0;
            lat.d_obs_n_sq_avg[tid] = 0;

            //Init sigma with random values
            lat.d_lattice_sigma[tid] = curand_normal(&randomState);
            lat.d_next_lattice_sigma[tid] = curand_normal(&randomState);
            lat.d_next_lattice_pi[tid] = curand_normal(&randomState);
            lat.d_lattice_pi[tid] = curand_normal(&randomState);

            //Create random-vector-field for diffusive random term
            for(int dim = 0; dim < parameters.dimension; dim++) {
                lat.d_random_noise[dim*parameters.volume + tid] = curand_normal(&randomState);
                lat.d_next_random_noise[dim*parameters.volume + tid] = curand_normal(&randomState);
            }

            //Initialize n with modelA random values
            lat.d_lattice_n[tid] = curand_normal(&randomState);
            lat.d_next_lattice_n[tid] = curand_normal(&randomState);
            lat.d_lattice_nu[tid] = 0;
            lat.d_next_lattice_nu[tid] = 0;
        }


        /**
         * Use the langevin-equation with the coupled Hamiltonian to calculate the next time-step.
         * A leapfrog-integrator scheme is used.
         * @param time
         * @param lat
         * @param parameters
         */
        __global__ void device_evolveLatticePoint(unsigned int time, Lattice lat, Parameters parameters) {
            //Check array range
            unsigned int tid = (blockIdx.x * blockDim.x) + threadIdx.x;
            if(tid >= parameters.volume) return;

            //Calculate indices
            unsigned int x,y,xf1,xf2,xb1,xb2,yf1,yf2,yb1,yb2,xf1yf1,xb1yf1,xf1yb1,xb1yb1;
            {
                y = tid % parameters.sites;
                x = tid / parameters.sites;

                xf1 = parameters.sites * ( (x+1) % parameters.sites) + y;
                xf2 = parameters.sites * ( (x+2) % parameters.sites) + y;
                xb1 = parameters.sites * ( (x-1+parameters.sites) % parameters.sites) + y;
                xb2 = parameters.sites * ( (x-2+parameters.sites) % parameters.sites) + y;

                yf1 = parameters.sites * x + ( (y+1) % parameters.sites);
                yf2 = parameters.sites * x + ( (y+2) % parameters.sites);
                yb1 = parameters.sites * x + ( (y-1+parameters.sites) % parameters.sites);
                yb2 = parameters.sites * x + ( (y-2+parameters.sites) % parameters.sites);

                xf1yf1 = parameters.sites * ( (x+1) % parameters.sites) + ( (y+1) % parameters.sites);
                xf1yb1 = parameters.sites * ( (x+1) % parameters.sites) + ( (y-1+parameters.sites) % parameters.sites);
                xb1yf1 = parameters.sites * ( (x-1+parameters.sites) % parameters.sites) + ( (y+1) % parameters.sites);
                xb1yb1 = parameters.sites * ( (x-1+parameters.sites) % parameters.sites) + ( (y-1+parameters.sites) % parameters.sites);
            }

            //Check for time to apply different parameters for the different stages of thermalization
            if(time == 0) {
                //Leapfrog half step
                parameters.integrator_timeDelta /= 2;
                parameters.integrator_timeDelta_SQRT /= 1.41421f;
            }else if(time < parameters.thermalization_time_1) {
                //External field for thermalization
                parameters.J = parameters.symmetryBreakingJ;
            }else if(parameters.fastThermalization ) {
                if(time-parameters.thermalization_time_1 < parameters.thermalization_time_2/10) {
                    parameters.integrator_timeDelta = 0.1;
                    parameters.integrator_timeDelta_SQRT = 0.3162278;
                }else if(time-parameters.thermalization_time_1 < parameters.thermalization_time_2/5) {
                    parameters.integrator_timeDelta = 0.1;
                    parameters.integrator_timeDelta_SQRT = 0.3162278;
                }
            }


            //Step sigma field (Model A dynamics)
            {
                //step pi
                TFloat next_pi = 0;

                //m^2 * phi
                next_pi -= parameters.msq_A * lat.d_lattice_sigma[tid];

                //nabla^2 phi
                next_pi +=
                        +(lat.d_lattice_sigma[xf1] - 4 * lat.d_lattice_sigma[tid] + lat.d_lattice_sigma[xb1])
                        +(lat.d_lattice_sigma[yf1]  + lat.d_lattice_sigma[yb1]);

                //lambda * phi^3
                next_pi -= parameters.lambda_A/6.0f * POWER_THREE(lat.d_lattice_sigma[tid]);

                //external field: J
                next_pi += parameters.J;

                //Langevin damping
                next_pi -= parameters.gamma_A * lat.d_lattice_pi[tid];

                //Coupling
                next_pi -= parameters.C * lat.d_lattice_n[tid];

                next_pi *= parameters.integrator_timeDelta;

                //Random term
                next_pi += parameters.integrator_timeDelta_SQRT * parameters.standard_deviation_A
                           * curand_normal(&((curandStateType *)lat.d_curandStateArray)[tid]);

                next_pi += lat.d_lattice_pi[tid];

                //Update pi,phi values for field
                lat.d_next_lattice_pi[tid] = next_pi;
                lat.d_next_lattice_sigma[tid] = lat.d_lattice_sigma[tid] + parameters.integrator_timeDelta * next_pi;
            }

            //Step n field with model A dynamics
            if(parameters.thermalization_time_3 == 0 || time < int(parameters.thermalization_time_1 + parameters.thermalization_time_2)) {
                //step pi (= nu)
                TFloat next_pi = 0;

                //m^2 * phi
                next_pi -= parameters.msq_B * lat.d_lattice_n[tid];

                //nabla^2 phi
                next_pi +=
                        +(lat.d_lattice_n[xf1] - 4 * lat.d_lattice_n[tid] + lat.d_lattice_n[xb1])
                        +(lat.d_lattice_n[yf1]  + lat.d_lattice_n[yb1]);

                //external field: J
                next_pi += parameters.J;

                next_pi -= parameters.gamma_B * lat.d_lattice_nu[tid];

                next_pi -= parameters.C * lat.d_lattice_sigma[tid];

                next_pi *= parameters.integrator_timeDelta;

                next_pi += parameters.integrator_timeDelta_SQRT * parameters.standard_deviation_A *
                           (curand_normal(&((curandStateType *)lat.d_curandStateArray)[tid]));

                next_pi += lat.d_lattice_nu[tid];

                //Update pi,phi values for field
                lat.d_next_lattice_nu[tid] = next_pi;
                lat.d_next_lattice_n[tid] = lat.d_lattice_n[tid] + parameters.integrator_timeDelta * next_pi;
            }

            else
                //Step n field with model B dynamics
            {
                //step nu
                TFloat next_nu = 0;

                //m^2 * nabla^2 n
                next_nu +=
                        parameters.msq_B * (
                                //nabla
                                lat.d_lattice_n[xf1] + lat.d_lattice_n[xb1] + lat.d_lattice_n[yf1] + lat.d_lattice_n[yb1] - 4 * lat.d_lattice_n[tid]
                        );

                next_nu -= 0
                           +       lat.d_lattice_n[xb2] + lat.d_lattice_n[xf2] + lat.d_lattice_n[yb2] + lat.d_lattice_n[yf2]
                           +  2 * (lat.d_lattice_n[xb1yf1]+lat.d_lattice_n[xf1yf1]+lat.d_lattice_n[xb1yb1]+lat.d_lattice_n[xf1yb1])
                           -  8 * (lat.d_lattice_n[xb1]+lat.d_lattice_n[xf1]+lat.d_lattice_n[yb1]+lat.d_lattice_n[yf1])
                           + 20 *  lat.d_lattice_n[tid];

                //Coupling
                next_nu +=
                        parameters.C * ( //nabla
                                +(lat.d_lattice_sigma[xf1] - 4 * lat.d_lattice_sigma[tid] + lat.d_lattice_sigma[xb1])
                                +(lat.d_lattice_sigma[yf1]  + lat.d_lattice_sigma[yb1])
                        );

                // mu * (...)
                next_nu *= parameters.mu;

                //Langevin damping
                next_nu -= parameters.gamma_B * lat.d_lattice_nu[tid];

                //leapfrog
                next_nu *= parameters.integrator_timeDelta;

                //randomNoise with sqrt(integrator_timeDelta)
                next_nu += parameters.integrator_timeDelta_SQRT * parameters.standard_deviation_B * (
                        lat.d_random_noise[tid] + lat.d_random_noise[parameters.volume + tid]
                        - lat.d_random_noise[xb1] - lat.d_random_noise[parameters.volume + yb1]
                );

                //leapfrog
                next_nu += lat.d_lattice_nu[tid];

                //Update pi,phi values for field (leapfrog)
                lat.d_next_lattice_nu[tid] = next_nu;
                lat.d_next_lattice_n[tid] = lat.d_lattice_n[tid] + parameters.integrator_timeDelta * next_nu;
            }

            //Fill next random vector field (for x and y dimension)
            lat.d_next_random_noise[tid] = curand_normal(&((curandStateType *)lat.d_curandStateArray)[tid]);
            lat.d_next_random_noise[parameters.volume + tid] = curand_normal(&((curandStateType *)lat.d_curandStateArray)[tid]);

            if(parameters.flag_printEveryStep) {
                //Measure observables (overwrite mode)
                lat.d_obs_sigma_avg[tid] = (lat.d_lattice_sigma[tid]);
                lat.d_obs_sigma_sq_avg[tid] = (POWER_TWO(lat.d_lattice_sigma[tid]));
                lat.d_obs_pi_avg[tid] = (lat.d_lattice_pi[tid]);
                lat.d_obs_pi_sq_avg[tid] = (POWER_TWO(lat.d_lattice_pi[tid]));
                lat.d_obs_n_avg[tid] = (lat.d_lattice_n[tid]);
                lat.d_obs_n_sq_avg[tid] = (POWER_TWO(lat.d_lattice_n[tid]));
                lat.d_obs_nu_sq_avg[tid] = (POWER_TWO(lat.d_lattice_nu[tid]));
                lat.d_obs_nu_avg[tid] = (lat.d_lattice_nu[tid]);
            }
            else
            if(time >= parameters.thermalization_time_3+parameters.thermalization_time_2+parameters.thermalization_time_1) {
                //Measure observables (add mode)
                lat.d_obs_sigma_avg[tid] += (lat.d_lattice_sigma[tid]);
                lat.d_obs_n_avg[tid] += (lat.d_lattice_n[tid]);
                lat.d_obs_sigma_sq_avg[tid] += (POWER_TWO(lat.d_lattice_sigma[tid]));
                lat.d_obs_n_sq_avg[tid] += (POWER_TWO(lat.d_lattice_n[tid]));
            }
        }


        /**
         * If the model dynamics are changed, the nu field needs to be reinitialized
         * @param lat
         * @param parameters
         */
        __global__ void device_reinitializeNuField(Lattice lat, Parameters parameters) {
            //Check array range
            unsigned int tid = (blockIdx.x * blockDim.x) + threadIdx.x;
            if(tid >= parameters.volume) return;

            //Calculate indexes
            unsigned int y = tid % parameters.sites;
            unsigned int x = tid / parameters.sites;
            unsigned int xb1 = parameters.sites * ( (x-1+parameters.sites) % parameters.sites) + y;

            //Divergence of the random-vector-field -> Diffusive random term
            lat.d_lattice_nu[tid] = sqrt(parameters.temperature) * (
                    lat.d_random_noise[tid] + lat.d_random_noise[parameters.volume + tid]
                    - lat.d_random_noise[xb1] - lat.d_random_noise[parameters.volume + xb1]
            );

            //Fill next random vector field (for x and y dimension)
            lat.d_next_random_noise[tid] = curand_normal(&((curandStateType *)lat.d_curandStateArray)[tid]);
            lat.d_next_random_noise[parameters.volume + tid] = curand_normal(&((curandStateType *)lat.d_curandStateArray)[tid]);
        }

        /**
         * Subtract the nu-average from every nu-lattice site.
         * @param lattice
         * @param parameters
         */
        void modelChange_subtractNuAverageFromNuField(const Lattice& lattice, const Parameters& parameters) {
            auto * buffer = new TFloat [parameters.volume];
            cudaMemcpy(buffer, lattice.d_lattice_nu, parameters.volume * sizeof(TFloat), cudaMemcpyDeviceToHost);

            //Use double precision for average calculation
            double nu_avg = 0;
            for(size_t i = 0; i < parameters.volume; i++) {
                nu_avg += (double) buffer[i];
            }
            nu_avg /= (double)parameters.volume;
            for(size_t i = 0; i < parameters.volume; i++) {
                buffer[i] -= (TFloat) nu_avg;
            }

            cudaMemcpy(lattice.d_lattice_nu, buffer, parameters.volume * sizeof(TFloat), cudaMemcpyHostToDevice);
            cudaDeviceSynchronize();
            delete[] buffer;
        }

        /**
         * Set the nu-field to zero
         * @param lattice
         * @param parameters
         */
        void modelChange_setNuFieldZero(const Lattice& lattice, const Parameters& parameters) {
            auto * buffer = new TFloat [parameters.volume];
            //Set local array to zero
            for(size_t i = 0; i < parameters.volume; i++) {buffer[i] = 0;}
            //Push local array on device
            cudaMemcpy(lattice.d_lattice_nu, buffer, parameters.volume * sizeof(TFloat), cudaMemcpyHostToDevice);
            cudaMemcpy(lattice.d_next_lattice_nu, buffer, parameters.volume * sizeof(TFloat), cudaMemcpyHostToDevice);
            cudaDeviceSynchronize();
            delete[] buffer;
            std::cout << "ZERO SET" << std::endl;
        }

    }

    namespace observables {
        //====== Observables
        __global__ void reduceArray(float * d_array, size_t length, size_t depth) {
            size_t tid = (blockIdx.x * blockDim.x) + threadIdx.x;

            //get indices for summation from tid and depth
            size_t first = (1 << depth) * tid;
            size_t second = (1 << depth) * (tid) + (1 << (depth-1));

            if(second >= length){
                //if only second index exceeds array, average firstIndex with zero
                if(first < length) { d_array[first] /= 2; }
                return;
            }
            //sum array
            d_array[first] += d_array[second];
            //average sum
            d_array[first] /= 2;
        }

        TFloat getFieldAverage(float * d_array, size_t length) {
            //smallest power-2 number greater than length; needed for correcting the value if length is not a power-2 number
            unsigned long long powerLength = 1;
            //loop as long as 2**depth < length
            for(int depth = 1; (1 << depth) < 2*length; depth++) {
                dim3 grid = (length/depth) / NUM_THREADS + 1;
                reduceArray<<<grid,NUM_THREADS>>>(d_array, length,depth);
                cudaDeviceSynchronize();
                powerLength <<= 1;
            }
            TFloat arraySum = 0;
            cudaMemcpy(&arraySum, d_array, sizeof(TFloat), cudaMemcpyDeviceToHost);
            //Correct value if length is not a power of two and return
            return arraySum * (float)((double)powerLength/(double)length);
        }

        /**
         * Get all observables form device using ArrayReduction
         * Caution: This modifies the observables Arrays!
         * TODO copy to buffer array in first step to avoid modification of the org array
         * @param lattice
         * @param parameters
         * @return
         */
        Observables getObservables(const Lattice& lattice, const Parameters& parameters) {
            Observables observables;

            observables.sigma_avg = getFieldAverage(lattice.d_obs_sigma_avg, parameters.volume);
            observables.sigma_sq_avg = getFieldAverage(lattice.d_obs_sigma_sq_avg, parameters.volume);
            observables.pi_avg = getFieldAverage(lattice.d_obs_pi_avg, parameters.volume);
            observables.pi_sq_avg = getFieldAverage(lattice.d_obs_pi_sq_avg, parameters.volume);
            observables.n_avg = getFieldAverage(lattice.d_obs_n_avg, parameters.volume);
            observables.n_sq_avg = getFieldAverage(lattice.d_obs_n_sq_avg, parameters.volume);
            observables.nu_avg = getFieldAverage(lattice.d_obs_nu_avg, parameters.volume);
            observables.nu_sq_avg = getFieldAverage(lattice.d_obs_nu_sq_avg, parameters.volume);

            cudaError_t cuda_err = cudaDeviceSynchronize();
            if (cuda_err != cudaSuccess) printf("cudaError: \"%s\".\n", cudaGetErrorString(cuda_err));

            return observables;
        }

    }


    /**
     * Starts a single langevin simulation with the given parameters and lattice
     * The lattice needs to be allocated. But not initialized
     * The device-lattices are changed
     * @param parameters Parameters for simulation
     * @param lattice Allocated lattice on the device
     * @return
     */
    Observables startSimulation(Parameters parameters, Lattice lattice) {

        //For storing all time steps if parameters.flag_printEveryStep is enabled
        std::map<std::string, std::vector<TFloat>> stepMap;
        {
            stepMap["sigma"] = std::vector<TFloat>();
            stepMap["sigma_sq"] = std::vector<TFloat>();
            stepMap["pi"] = std::vector<TFloat>();
            stepMap["pi_sq"] = std::vector<TFloat>();
            stepMap["n"] = std::vector<TFloat>();
            stepMap["n_sq"] = std::vector<TFloat>();
            stepMap["nu"] = std::vector<TFloat>();
            stepMap["nu_sq"] = std::vector<TFloat>();
        }

        //Number of Cuda-cores to invoke
        dim3 grid = parameters.volume / NUM_THREADS + 1;

        //Init lattice
        langevin::device_initLattice<<<grid,NUM_THREADS>>>(lattice, parameters);
        cudaDeviceSynchronize();

        //Set external field for thermalization
        //If programOption is zero choose randomly between +0.1 and - 0.1
        if(parameters.symmetryBreakingJ == 0) {
            std::default_random_engine generator(parameters.seed);
            std::uniform_real_distribution<TFloat> distribution(0,1);
            parameters.symmetryBreakingJ = distribution(generator) < 0.5f ? -1 : 1;
        }

        //Start time evolution
        for(int time = 0 ; time < parameters.thermalization_time_1 + parameters.thermalization_time_2
                                + parameters.thermalization_time_3 + parameters.measure_time; time++) {

            //Thermalization ModelChange correction
            if(parameters.thermalization_time_3 != 0 &&
                time == parameters.thermalization_time_1 + parameters.thermalization_time_2) {
                if(parameters.dynamicChangeMode == DYNAMIC_CHANGE_MODE_RE_INIT) {
                    langevin::device_reinitializeNuField<<<grid, NUM_THREADS>>>(lattice, parameters);
                    cudaDeviceSynchronize();
                }else if(parameters.dynamicChangeMode == DYNAMIC_CHANGE_MODE_AVG_SUBTRACT) {
                    langevin::modelChange_subtractNuAverageFromNuField(lattice, parameters);
                }else if(parameters.dynamicChangeMode == DYNAMIC_CHANGE_MODE_ZERO) {
                    langevin::modelChange_setNuFieldZero(lattice,parameters);
                }
            }

            //Update every lattice-site
            langevin::device_evolveLatticePoint<<<grid, NUM_THREADS>>>(time, lattice,parameters);
            cudaError_t cuda_err = cudaDeviceSynchronize();
            if (cuda_err != cudaSuccess) {
                printf("cudaError: \"%s\".\n", cudaGetErrorString(cuda_err));
                std::exit(1);
            }

            //Swap field-pointers for next calculation
            std::swap(lattice.d_lattice_sigma, lattice.d_next_lattice_sigma);
            std::swap(lattice.d_lattice_pi, lattice.d_next_lattice_pi);
            std::swap(lattice.d_lattice_n, lattice.d_next_lattice_n);
            std::swap(lattice.d_lattice_nu, lattice.d_next_lattice_nu);
            std::swap(lattice.d_random_noise, lattice.d_next_random_noise);

            if(parameters.flag_printEveryStep) {
                Observables stepObs = observables::getObservables(lattice,parameters);

                std::cout << "t=" << time <<"\t";
                printObservables(stepObs);

                //Record sigma and n
                stepMap["sigma"].push_back(stepObs.sigma_avg);
                stepMap["sigma_sq"].push_back(stepObs.sigma_sq_avg);
                stepMap["pi"].push_back(stepObs.pi_avg);
                stepMap["pi_sq"].push_back(stepObs.pi_sq_avg);
                stepMap["n"].push_back(stepObs.n_avg);
                stepMap["n_sq"].push_back(stepObs.n_sq_avg);
                stepMap["nu"].push_back(stepObs.nu_avg);
                stepMap["nu_sq"].push_back(stepObs.nu_sq_avg);
            }
        }

        //At this point the simulation is finished

        //Get observables from lattice after measurement
        Observables observables = observables::getObservables(lattice, parameters);

        //Normalize observables
        if(!parameters.flag_printEveryStep) {
            observables.sigma_avg /= (TFloat) parameters.measure_time;
            observables.sigma_sq_avg /= (TFloat) parameters.measure_time;
            observables.n_avg /= (TFloat) parameters.measure_time;
            observables.n_sq_avg /= (TFloat) parameters.measure_time;
        } else {
            writeAllSteps(parameters, stepMap);
        }

        return observables;
    }

}
