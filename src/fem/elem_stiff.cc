#include <fem/elem_stiff.h>

namespace fem {
template <typename T>
void lin_trans_mat_tl(
  const T* const h0, const T* const u0t, const unsigned int N, T* const B0t_L) {
  auto _3N = 3 * N;
  auto B0 = B0t_L;
  auto B1 = B0 + _3N;
  auto B2 = B1 + _3N;
  auto B3 = B2 + _3N;
  auto B4 = B3 + _3N;
  auto B5 = B4 + _3N;
  auto H = h0;
  for (unsigned int i = 0; i < N; ++i) {
    B0[0] = H[0] * ((T)1 + u0t[0]);
    B0[1] = H[0] * u0t[3];
    B0[2] = H[0] * u0t[6];
    B0 += 3;
    B1[0] = H[1] * u0t[1];
    B1[1] = H[1] * ((T)1 + u0t[4]);
    B1[2] = H[1] * u0t[7];
    B1 += 3;
    B2[0] = H[2] * u0t[2];
    B2[1] = H[2] * u0t[5];
    B2[2] = H[2] * ((T)1 + u0t[8]);
    B2 += 3;
    B3[0] = H[0]*u0t[1] + H[1]*((T)1+u0t[0]);
    B3[1] = H[0]*((T)1+u0t[4]) + H[1]*u0t[3];
    B3[2] = H[0]*u0t[7] + H[1]*u0t[6];
    B3 += 3;
    B4[0] = H[1]*u0t[2] + H[2]*u0t[1];
    B4[1] = H[1]*u0t[5] + H[2]*((T)1+u0t[4]);
    B4[2] = H[1]*((T)1+u0t[8]) + H[2]*u0t[7];
    B4 += 3;
    B5[0] = H[0]*u0t[2] + H[2]*((T)1+u0t[0]);
    B5[1] = H[0]*u0t[5] + H[2]*u0t[3];
    B5[2] = H[0]*((T)1+u0t[8]) + H[2]*u0t[6];
    B5 += 3;
    H += 3;
  }
}
template void lin_trans_mat_tl(
  const double* const h0, const double* const u0t,
  const unsigned int N, double* const B0t_L);
} // namespace fem
