#ifndef FEM_MAT_ELASTIC_H_
#define FEM_MAT_ELASTIC_H_

#include "material.h"

namespace fem {
template <typename T>
struct Ela3DInitializer : MatInitializer<T> {
  static inline void initialize(T* const C, T* const p) {
    auto E = p[0];
    auto nu = p[1];
    for (unsigned int i = 0; i < 36; ++i) {
      C[i] = (T)0;
    }
    T lambda = (T)(E*nu/((1+nu)*(1-2*nu)));
    T mu = (T)(E/(2*(1+nu)));
    C[0] = C[7] = C[14] = (T)(lambda+2*mu);
    C[1] = C[2] = C[6] = C[8] = C[12] = C[13] = lambda;
    C[21] = C[28] = C[35] = mu;
  }
};

template <typename T>
using Ela3D = Mat<T,6,2,Ela3DInitializer>;

} // namespace fem

#endif // FEM_MAT_ELASTIC_H_
