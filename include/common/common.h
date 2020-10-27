#ifndef COMMON_COMMON_H_
#define COMMON_COMMON_H_

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <glog/logging.h>

#define DP_ZERO 1e-6

template <typename T>
void init_rand(T* const data, const unsigned int size) {
  for (unsigned int i = 0; i < size; ++i) {
    data[i] = (T)(rand() & 0xFF) / 10.0;
  }
}

template <typename T>
void init_arange(T* const data, const unsigned int size) {
  for (unsigned int i = 0; i < size; ++i) {
    data[i] = (T)(2*i + 1);
  }
}

template <typename T>
void init_zero(T* const data, const unsigned int size) {
  for (unsigned int i = 0; i < size; ++i) {
    data[i] = (T)0;
  }
}

template <typename T>
void init_diag_unit(T* const data, const unsigned int ndim) {
  auto d_ = data;
  for (unsigned int i = 0; i < ndim; ++i) {
    for (unsigned int j = 0; j < ndim; ++j) {
      d_[j] = (T)(i==j ? 1 : 0);
    }
    d_ += ndim;
  }
}

template <typename T>
T rand_gen(const T low, const T high) {
  T rnd = static_cast<T>(rand())/static_cast<T>(RAND_MAX);
  return low + rnd*(high-low);
}

template <typename T>
void print_mat(
  const T* const a, const unsigned int nrow, const unsigned int ncol) {
  auto a_ = a;
  for (unsigned int i = 0; i < nrow; ++i) {
    for (unsigned int j = 0; j < ncol; ++j) {
      std::cout << a_[j] << " ";
    }
    std::cout << std::endl;
    a_ += ncol;
  }
  std::cout << std::endl;
}

template <typename T>
bool validate(
  const T* const a, const T* const b,
  const unsigned int size, double tol=1e-6) {
  for (unsigned int i = 0; i < size; ++i) {
    if ((a[i]<DP_ZERO) && (b[i]<DP_ZERO)) {
      continue;
    }
    double error = abs(a[i]/b[i]-1);
    if (error > tol) {
      LOG(WARNING) << "out of tolerance, a[i]: " << a[i] << ", b[i]: " 
        << b[i] << ", relative error: " << error << ", tol: " << tol;
      return false;
    }
  }
  return true;
}

#endif // COMMON_COMMON_H_
