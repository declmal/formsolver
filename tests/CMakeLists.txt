cmake_minimum_required(VERSION 3.2)

project(formsolver_tests)

get_filename_component(PARENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/config.cmake)
  include(${CMAKE_CURRENT_SOURCE_DIR}/config.cmake)
elseif(EXISTS ${PARENT_DIR}/cmake/config.cmake)
  include(${PARENT_DIR}/cmake/config.cmake)
endif()

# set compilation flags
set(CMAKE_CXX_FLAGS "-std=c++14")

# specify include directories
include_directories("../include")

# find cblas
include(${PARENT_DIR}/cmake/FindBLAS.cmake)

# find glog
find_package(glog REQUIRED)

# locate libform
find_file(LIB_FORM_PATH NAMES libform.so libform.dylib PATHS ../lib)
get_filename_component(LIB_FORM ${LIB_FORM_PATH} ABSOLUTE)

# locate libccx
find_file(LIB_CCX_PATH NAMES libccx.so libccx.dylib PATHS ../lib)
get_filename_component(LIB_CCX ${LIB_CCX_PATH} ABSOLUTE)

# compile executives and link libraries
file(GLOB_RECURSE SRCS_FORM_TEST *.cc)
foreach (SRC_FILE ${SRCS_FORM_TEST})
  get_filename_component(EXE_FILE ${SRC_FILE} NAME_WE)
  message("compling test case: " ${SRC_FILE} " to executive: " ${EXE_FILE})
  add_executable(${EXE_FILE} ${SRC_FILE})
  target_link_libraries(${EXE_FILE} 
                        PRIVATE ${LIB_FORM} 
                        PRIVATE ${LIB_CCX}
                        PRIVATE glog::glog
                        PRIVATE ${BLAS_LIBRARIES})
endforeach()
