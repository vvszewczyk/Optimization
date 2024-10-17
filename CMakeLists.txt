cmake_minimum_required(VERSION 3.29)
project(Optimization_s)

set(CMAKE_CXX_STANDARD 14)

include_directories(optimization/optimization)

add_executable(Optimization_s
    optimization/optimization/main.cpp
    optimization/optimization/matrix.cpp
    optimization/optimization/matrix.h
    optimization/optimization/ode_solver.cpp
    optimization/optimization/ode_solver.h
    optimization/optimization/opt_alg.cpp
    optimization/optimization/opt_alg.h
    optimization/optimization/solution.cpp
    optimization/optimization/solution.h
    optimization/optimization/user_funs.cpp
    optimization/optimization/user_funs.h)
