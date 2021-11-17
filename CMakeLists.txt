cmake_minimum_required(VERSION 3.20)

## not recommended
set(CMAKE_C_COMPILER "gcc-8")
set(CMAKE_CXX_COMPILER "g++-8")


set(CMAKE_CUDA_COMPILER  /usr/local/cuda/bin/nvcc)

project(coupledDynamics LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 14)

find_package(Boost 1.76 REQUIRED COMPONENTS log program_options)
include_directories(${Boost_INCLUDE_DIRS})
set(CMAKE_CXX_BOOST_FLAGS "-DBOOST_LOG_DYN_LINK")

add_executable(coupledDynamics main.cu InputOutput.cpp InputOutput.h Lattice2D.cu Lattice2D.cuh global.h)

target_link_libraries(coupledDynamics ${Boost_LIBRARIES})

set_target_properties(coupledDynamics PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)




