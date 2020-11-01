cmake_minimum_required(VERSION 3.2)

project(formsolver)

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/config.cmake)
  include(${CMAKE_CURRENT_SOURCE_DIR}/config.cmake)
elseif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/cmake/config.cmake)
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/config.cmake)
endif()

# set compilation flags
set(CMAKE_CXX_FLAGS "-std=c++17")

# specify include directories
include_directories("include")

# find cblas
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindBLAS.cmake)

# find glog
find_package(glog REQUIRED)

# compile shared library
set(LIB_FORM form)
file(GLOB_RECURSE SRCS_FORM src/*.cc)
add_library(${LIB_FORM} SHARED ${SRCS_FORM})

# link libraries
target_link_libraries(${LIB_FORM} 
                      PRIVATE glog::glog
                      PRIVATE ${BLAS_LIBRARIES})
