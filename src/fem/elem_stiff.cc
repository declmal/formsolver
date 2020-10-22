#include <fem/elem_stiff.h>

namespace fem {
template <typename T>
void jacobi(const T* const h, const T* const X, const unsigned int N, T* const J) {
  for (unsigned int i = 0; i < N; ++i) {
  }
}

template <typename T>
void det_jacobi(const T* const J, T* const detJ) {
  detJ[0] = 
    (J[3]*J[7] - J[4]*J[6]) * J[2] +
    (J[1]*J[6] - J[0]*J[7]) * J[5] +
    (J[0]*J[4] - J[1]*J[3]) * J[8];
}

template <typename T>
void inv_jacobi(const T* const J, const T* const detJ, T* const Jinv) {
  Jinv[0] = J[4]*J[8] - J[5]*J[7];
  Jinv[1] = J[2]*J[7] - J[1]*J[8];
  Jinv[2] = J[1]*J[5] - J[2]*J[4];
  Jinv[3] = J[5]*J[6] - J[3]*J[8];
  Jinv[4] = J[0]*J[8] - J[2]*J[6];
  Jinv[5] = J[2]*J[3] - J[0]*J[5];
  Jinv[6] = J[3]*J[7] - J[4]*J[6];
  Jinv[7] = J[1]*J[6] - J[0]*J[7];
  Jinv[8] = J[0]*J[4] - J[1]*J[3];
  for (unsigned int i = 0; i < 8; ++i) {
    Jinv /= detJ[0];
  }
}

template <typename T>
void intr_deriv(
  const T* const Jinv, const T* const h, const unsigned int N, T* const htau) {
  auto Htau = htau;
  auto H = h;
  for (unsigned int i = 0; i < N; ++i) {
    Htau[0] = Jinv[0]*H[0] + Jinv[1]*H[1] + Jinv[2]*H[2];
    Htau[1] = Jinv[3]*H[0] + Jinv[4]*H[1] + Jinv[5]*H[2];
    Htau[2] = Jinv[6]*H[0] + Jinv[7]*H[1] + Jinv[8]*H[2];
    H += 3;
    Htau += 3;
  }
}

template <typename T>
void disp_deriv_tl(
  const T* const h0, const T* const Ut, const unsigned int N, T* const u0t) {
  auto u = u0t;
  auto U0 = Ut;
  auto U1 = U0 + N;
  auto U2 = U1 + N;
  auto H = h0;
  for (unsigned int i = 0; i < 9; ++i) {
    u0t[i] = (T)0;
  }
  for (unsigned int i = 0; i < N; ++i) {
    u0t[0] += H[0] * U0[0];
    u0t[1] += H[1] * U0[0];
    u0t[2] += H[2] * U0[0];
    U0 += 1;
    u0t[3] += H[0] * U1[0];
    u0t[4] += H[1] * U1[0];
    u0t[5] += H[2] * U1[0];
    U1 += 1;
    u0t[6] += H[0] * U2[0];
    u0t[7] += H[1] * U2[0];
    u0t[8] += H[2] * U2[0];
    U2 += 1;
    H += 3;
  }
}

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
} // namespace fem
