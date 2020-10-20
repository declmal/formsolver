cmake_minimum_required(VERSION 3.2)

# specify include directories
include_directories("include")

# compile executives
find_package(CUDA REQUIRED)
cuda_add_executable(test_bin tests/test_fem.cc)
