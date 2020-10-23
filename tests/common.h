#ifndef TESTS_COMMON_H_
#define TESTS_COMMON_H_

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

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

void print_mat_dp(const double* const a, const unsigned int nrow, const unsigned int ncol) {
  auto a_ = a;
  for (unsigned int i = 0; i < nrow; ++i) {
    for (unsigned int j = 0; j < ncol; ++j) {
      printf("%10lf ", a_[j]);
    }
    printf("\n");
    a_ += ncol;
  }
  printf("\n");
}

bool validate_dp(
  const double* const a, const double* const b,
  const unsigned int size, double tol=1e-6) {
  for (unsigned int i = 0; i < size; ++i) {
    if (abs(a[i]-b[i]) > tol) {
      return false;
    }
  }
  return true;
}

#endif // TESTS_COMMON_H_
