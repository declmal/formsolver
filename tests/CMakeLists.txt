cmake_minimum_required(VERSION 3.2)

project(formsolver_tests)

if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/config.cmake)
  include(${CMAKE_CURRENT_BINARY_DIR}/config.cmake)
elseif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/config.cmake)
  include(${CMAKE_CURRENT_SOURCE_DIR}/config.cmake)
endif()

# set compilation flags
set(CMAKE_CXX_FLAGS "-std=c++11")

# specify include directories
include_directories("../include")

# find extern packages
find_package(glog REQUIRED)
if(USE_BLAS STREQUAL "openblas")
  find_package(BLAS REQUIRED)
elseif(USE_BLAS STREQUAL "apple")
  find_library(BLAS_LIBRARIES Accelerate)
  include_directories(SYSTEM ${BLAS_LIBRARIES}/Versions/Current/Frameworks/vecLib.framework/Versions/Current/Headers/)
endif()
message(STATUS "Using BLAS library " ${BLAS_LIBRARIES})

# locate libform
find_file(LIB_FORM_PATH NAMES libform.so libform.dylib PATHS ../lib)
get_filename_component(LIB_FORM ${LIB_FORM_PATH} ABSOLUTE)

# compile executives and link libraries
file(GLOB_RECURSE SRCS_FORM_TEST *.cc)
foreach (SRC_FILE ${SRCS_FORM_TEST})
  get_filename_component(EXE_FILE ${SRC_FILE} NAME_WE)
  message("compling test case: " ${SRC_FILE} " to executive: " ${EXE_FILE})
  add_executable(${EXE_FILE} ${SRC_FILE})
  target_link_libraries(${EXE_FILE} 
                        PRIVATE ${LIB_FORM} 
                        PRIVATE glog::glog
                        PRIVATE ${BLAS_LIBRARIES})
endforeach()
