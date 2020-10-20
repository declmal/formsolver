cmake_minimum_required(VERSION 3.2)

project(formsolver_tests)

# set compilation flags
set(CMAKE_CXX_FLAGS "-std=c++11")

# specify include directories
include_directories("../include")

# find extern packages
find_package(glog REQUIRED)

# locate libform
find_file(LIB_FORM_PATH NAMES libform.so libform.dylib PATHS ../lib)
get_filename_component(LIB_FORM ${LIB_FORM_PATH} ABSOLUTE)

# compile executives and link libraries
file(GLOB_RECURSE SRCS_FORM_TEST *.cc)
foreach (SRC_FILE ${SRCS_FORM_TEST})
  find_file(SRC_FILE_PATH NAMES ${SRC_FILE} PATHS .)
  get_filename_component(EXE_FILE ${SRC_FILE_PATH} ABSOLUTE)
  add_executable(${EXE_FILE} ${SRC_FILE_PATH})
  target_link_libraries(${EXES_FILE} ${LIB_FORM})
  target_link_libraries(${EXES_FILE} PRIVATE glog::glog)
endforeach()
