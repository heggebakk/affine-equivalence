cmake_minimum_required(VERSION 3.20)
project(affine C)

set(CMAKE_C_STANDARD 11)

add_executable(affine main.c structures.c structures.h affine.c affine.h)
