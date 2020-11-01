cmake_minimum_required(VERSION 3.2)

project(formsolver_ccx)

# enbale fortran
enable_language(Fortran)

# compile shared library
set(LIB_CCX ccx)
file(GLOB_RECURSE SRCS_CCX form_*.f)
add_library(${LIB_CCX} SHARED ${SRCS_CCX})
