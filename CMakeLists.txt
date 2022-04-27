cmake_minimum_required(VERSION 3.20)
project(affine C)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0") # Profiling
set(CMAKE_C_STANDARD 11) # Profiling

add_executable(affine src/main.c src/structures.c src/structures.h src/affine.c src/affine.h src/orthoderivative.c src/orthoderivative.h src/adjoint.c src/adjoint.h)
