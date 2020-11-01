#ifndef TESTS_FEM_BRICK_COMMON_BRICK_H_
#define TESTS_FEM_BRICK_COMMON_BRICK_H_

#include <common/common.h>

template <typename T>
void init_x(T* const X, unsigned int NR, unsigned int NC) {
  if ((NR == 3) && (NC == 8)) {
    T X0[24] = {
      1,0,0,1,1,0,0,1,
      1,1,0,0,1,1,0,0,
      1,1,1,1,0,0,0,0
    };
    for (unsigned int i = 0; i < 24; ++i) {
      X[i] = X0[i];
    }
  } else if ((NR == 3) && (NC == 20)) {
    T X0[60] = {
      0, 0, 0,
      1, 0, 0,
      1, 1, 0,
      0, 1, 0,
      0, 0, 1,
      1, 0, 1,
      1, 1, 1,
      0, 1, 1,
      0.5, 0, 0,
      1, 0.5, 0,
      0.5, 1, 0,
      0, 0.5, 0,
      0.5, 0, 1,
      1, 0.5, 1,
      0.5, 1, 1,
      0, 0.5, 1,
      0, 0, 0.5,
      1, 0, 0.5,
      1, 1, 0.5,
      0, 1, 0.5,
    };
    T X0t[60];
    transpose<T>(X0, 20, 3, X0t);
    for (unsigned int i = 0; i < 60; ++i) {
      X[i] = X0t[i];
    }
  }
}

#endif  // TESTS_FEM_BRICK_COMMON_BRICK_H_
