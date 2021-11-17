//
// Created by maximilian on 30.08.21.
//

#ifndef MODELAB_CUDA_INPUTOUTPUT_H
#define MODELAB_CUDA_INPUTOUTPUT_H

#include "global.h"
#include <boost/program_options.hpp>
#include <iostream>
#include <chrono>
#include <fstream>
#include <iomanip>

Parameters getParametersFromProgramOptions(char** argc, int argv);

void printParameters(const Parameters& parameters);

void printObservables(Observables observables);

void printRuntime(unsigned long long milliSeconds);

void writeToDatabase(const Parameters parameters, const std::vector<Observables>& observablesList);

void writeAllSteps(const Parameters& parameters, const std::map<std::string, std::vector<TFloat>>& floatList);

#endif //MODELAB_CUDA_INPUTOUTPUT_H