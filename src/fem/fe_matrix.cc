#include <fem/fe_matrix.h>

namespace fem {
template <typename T>
void matmul_3nn3(
  const T* const a, const T* const b, const unsigned int N, T* const c) {
  auto a0_ = a;
  auto b_ = b;
  auto a1_ = a0_ + N;
  auto a2_ = a1_ + N;
  for (unsigned int k = 0; k < N; ++k) {
    c[0] += a0_[0] * b_[0];
    c[1] += a0_[0] * b_[1];
    c[2] += a0_[0] * b_[2];
    a0_ += 1;
    c[3] += a1_[0] * b_[0];
    c[4] += a1_[0] * b_[1];
    c[5] += a1_[0] * b_[2];
    a1_ += 1;
    c[6] += a2_[0] * b_[0];
    c[7] += a2_[0] * b_[1];
    c[8] += a2_[0] * b_[2];
    a2_ += 1;
    b_ += 3;
  }
}
template void matmul_3nn3(
  const double* const a, const double* const b,
  const unsigned int N, double* const c);

template <typename T>
void matmul_n333(
  const T* const a, const T* const b, const unsigned int N, T* const c) {
  auto a_ = a;
  auto c_ = c;
  for (unsigned int k = 0; k < N; ++k) {
    c_[0] = a_[0]*b[0] + a_[1]*b[3] + a_[2]*b[6];
    c_[1] = a_[0]*b[1] + a_[1]*b[4] + a_[2]*b[7];
    c_[2] = a_[0]*b[2] + a_[1]*b[5] + a_[2]*b[8];
    a_ += 3;
    c_ += 3;
  }
}
template void matmul_n333(
  const double* const a, const double* const b,
  const unsigned int N, double* const c);

template <typename T>
void det_33(const T* const a, T* const det) {
  det[0] = 
    (a[3]*a[7] - a[4]*a[6]) * a[2] +
    (a[1]*a[6] - a[0]*a[7]) * a[5] +
    (a[0]*a[4] - a[1]*a[3]) * a[8];
}
template void det_33(const double* const a, double* const det);

template <typename T>
void inv_33(const T* const a, const T* const det, T* const inv) {
  inv[0] = a[4]*a[8] - a[5]*a[7];
  inv[1] = a[2]*a[7] - a[1]*a[8];
  inv[2] = a[1]*a[5] - a[2]*a[4];
  inv[3] = a[5]*a[6] - a[3]*a[8];
  inv[4] = a[0]*a[8] - a[2]*a[6];
  inv[5] = a[2]*a[3] - a[0]*a[5];
  inv[6] = a[3]*a[7] - a[4]*a[6];
  inv[7] = a[1]*a[6] - a[0]*a[7];
  inv[8] = a[0]*a[4] - a[1]*a[3];
  for (unsigned int i = 0; i < 8; ++i) {
    inv[i] /= det[0];
  }
}
template void inv_33(
  const double* const a, const double* const det, double* const inv);
} // namespace fem
