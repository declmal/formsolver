cmake_minimum_required(VERSION 3.2)

project(formsolver)

# set compilation flags
set(CMAKE_CXX_FLAGS "-std=c++14")

# specify include directories
include_directories("include")

# # find extern packages
# find_package(glog REQUIRED)

# compile shared library
set(LIB_FORM form)
file(GLOB_RECURSE SRCS_FORM src/*.cc)
add_library(${LIB_FORM} SHARED ${SRCS_FORM})

# # link libraries
# target_link_libraries(${LIB_FORM} PRIVATE glog::glog)
