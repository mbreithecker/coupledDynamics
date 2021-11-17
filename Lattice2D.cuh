//
// Created by maximilian on 30.08.21.
//

#ifndef MODELAB_CUDA_LATTICE2D_CUH
#define MODELAB_CUDA_LATTICE2D_CUH


#include "global.h"
#include <curand_kernel.h>
#include <iostream>
#include <utility>
#include "InputOutput.h"

namespace Lattice2D {
    Observables startSimulation(Parameters parameters, Lattice lattice);
}


#endif //MODELAB_CUDA_LATTICE2D_CUH
