if(USE_BLAS STREQUAL "openblas")
  find_package(BLAS REQUIRED)
elseif(USE_BLAS STREQUAL "apple")
  find_library(BLAS_LIBRARIES Accelerate)
  include_directories(SYSTEM ${BLAS_LIBRARIES}/Versions/Current/Frameworks/vecLib.framework/Versions/Current/Headers/)
endif()
message(STATUS "Using BLAS library " ${BLAS_LIBRARIES})