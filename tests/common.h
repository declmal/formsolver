#ifndef TESTS_COMMON_H_
#define TESTS_COMMON_H_

#include <stdlib.h>
#include <stdio.h>

template <typename T>
void rand_init(T* const data, const unsigned int size) {
  for (unsigned int i = 0; i < size; ++i) {
    data[i] = (T)(rand() & 0xFF) / 10.0;
  }
}

template <typename T>
void print_mat(const T* const a, const unsigned int nrow, const unsigned int ncol) {
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

#endif // TESTS_COMMON_H_
