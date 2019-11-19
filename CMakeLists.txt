cmake_minimum_required(VERSION 3.5)
project(gauss)

set(CMAKE_CXX_STANDARD 14)

set(HEADERS gauss.h iterSystemSolver.h SubFunctions.h Classic_Newton.h RK4.h residual.h)
set(SOURCES main.cpp gauss.cpp iterSystemSolver.cpp SubFunctions.cpp Classic_Newton.cpp RK4.cpp residual.cpp)
add_executable(gauss ${HEADERS} ${SOURCES})