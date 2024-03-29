#include <common/matrix.h>

template <typename T>
void matmul_dnnd_3d(
  const T* const a, const T* const b, const unsigned int N, T* const c) {
  auto a0_ = a;
  auto a1_ = a0_ + N;
  auto a2_ = a1_ + N;
  auto b_ = b;
  for (unsigned int k = 0; k < 9; ++k) {
    c[k] = 0; 
  }
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

template <typename T>
void matmul_nddd_3d(
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

template <typename T>
void matmul_ndddt_3d(
  const T* const a, const T* const b, const unsigned int N, T* const c) {
  auto a_ = a;
  auto c_ = c;
  for (unsigned int k = 0; k < N; ++k) {
    c_[0] = a_[0]*b[0] + a_[1]*b[1] + a_[2]*b[2];
    c_[1] = a_[0]*b[3] + a_[1]*b[4] + a_[2]*b[5];
    c_[2] = a_[0]*b[6] + a_[1]*b[7] + a_[2]*b[8];
    a_ += 3;
    c_ += 3;
  }
}

template <typename T>
T det_dd_3d(const T* const a) {
  return
    (a[3]*a[7] - a[4]*a[6]) * a[2] +
    (a[1]*a[6] - a[0]*a[7]) * a[5] +
    (a[0]*a[4] - a[1]*a[3]) * a[8];
}

template <typename T>
void inv_dd_3d(const T* const a, const T det, T* const inv) {
  inv[0] = a[4]*a[8] - a[5]*a[7];
  inv[1] = a[2]*a[7] - a[1]*a[8];
  inv[2] = a[1]*a[5] - a[2]*a[4];
  inv[3] = a[5]*a[6] - a[3]*a[8];
  inv[4] = a[0]*a[8] - a[2]*a[6];
  inv[5] = a[2]*a[3] - a[0]*a[5];
  inv[6] = a[3]*a[7] - a[4]*a[6];
  inv[7] = a[1]*a[6] - a[0]*a[7];
  inv[8] = a[0]*a[4] - a[1]*a[3];
  for (unsigned int i = 0; i < 9; ++i) {
    inv[i] /= det;
  }
}
template <typename T>
void mattile_diag_dd_3d(const T* const a, T* const tile) {
  for (unsigned int i = 0; i < 81; ++i) {
    tile[i] = (T)0;
  }
  tile[0] = tile[30] = tile[60] = a[0];
  tile[1] = tile[31] = tile[61] = a[1];
  tile[2] = tile[32] = tile[62] = a[2];
  tile[9] = tile[39] = tile[69] = a[3];
  tile[10] = tile[40] = tile[70] = a[4];
  tile[11] = tile[41] = tile[71] = a[5];
  tile[18] = tile[48] = tile[78] = a[6];
  tile[19] = tile[49] = tile[79] = a[7];
  tile[20] = tile[50] = tile[80] = a[8];
}
template <typename T>
void matmul2_dne_ee_edn_3d(
  const T* const a, const T* const b, const unsigned int N,
  T* const buffer, T* const c) {
  auto _3N = 3 * N;
  auto a0 = a;
  auto a1 = a0 + _3N;
  auto a2 = a1 + _3N;
  auto a3 = a2 + _3N;
  auto a4 = a3 + _3N;
  auto a5 = a4 + _3N;
  auto c_ = c;
  for (unsigned int i = 0; i < _3N; ++i) {
    buffer[0] = a0[i]*b[0] + a1[i]*b[6] + a2[i]*b[12] +
                a3[i]*b[18] + a4[i]*b[24] + a5[i]*b[30];
    buffer[1] = a0[i]*b[1] + a1[i]*b[7] + a2[i]*b[13] +
                a3[i]*b[19] + a4[i]*b[25] + a5[i]*b[31];
    buffer[2] = a0[i]*b[2] + a1[i]*b[8] + a2[i]*b[14] +
                a3[i]*b[20] + a4[i]*b[26] + a5[i]*b[32];
    buffer[3] = a0[i]*b[3] + a1[i]*b[9] + a2[i]*b[15] +
                a3[i]*b[21] + a4[i]*b[27] + a5[i]*b[33];
    buffer[4] = a0[i]*b[4] + a1[i]*b[10] + a2[i]*b[16] +
                a3[i]*b[22] + a4[i]*b[28] + a5[i]*b[34];
    buffer[5] = a0[i]*b[5] + a1[i]*b[11] + a2[i]*b[17] +
                a3[i]*b[23] + a4[i]*b[29] + a5[i]*b[35];
    for (unsigned int j = 0; j < _3N; ++j) {
      c_[j] = buffer[0]*a0[j] + buffer[1]*a1[j] + buffer[2]*a2[j] +
              buffer[3]*a3[j] + buffer[4]*a4[j] + buffer[5]*a5[j];
    }
    c_ += _3N;
  }
}

template <typename T>
void matmul2_dnf_ff_fdn_3d(
  const T* const a, const T* const b, const unsigned int N,
  T* const buffer, T* const c) {
  auto _3N = 3 * N;
  auto a0 = a;
  auto a1 = a0 + _3N;
  auto a2 = a1 + _3N;
  auto a3 = a2 + _3N;
  auto a4 = a3 + _3N;
  auto a5 = a4 + _3N;
  auto a6 = a5 + _3N;
  auto a7 = a6 + _3N;
  auto a8 = a7 + _3N;
  auto c_ = c;
  for (unsigned int i = 0; i < _3N; ++i) {
    buffer[0] = a0[i]*b[0] + a1[i]*b[9] + a2[i]*b[18] +
                a3[i]*b[27] + a4[i]*b[36] + a5[i]*b[45] +
                a6[i]*b[54] + a7[i]*b[63] + a8[i]*b[72];
    buffer[1] = a0[i]*b[1] + a1[i]*b[10] + a2[i]*b[19] +
                a3[i]*b[28] + a4[i]*b[37] + a5[i]*b[46] +
                a6[i]*b[55] + a7[i]*b[64] + a8[i]*b[73];
    buffer[2] = a0[i]*b[2] + a1[i]*b[11] + a2[i]*b[20] +
                a3[i]*b[29] + a4[i]*b[38] + a5[i]*b[47] +
                a6[i]*b[56] + a7[i]*b[65] + a8[i]*b[74];
    buffer[3] = a0[i]*b[3] + a1[i]*b[12] + a2[i]*b[21] +
                a3[i]*b[30] + a4[i]*b[39] + a5[i]*b[48] +
                a6[i]*b[57] + a7[i]*b[66] + a8[i]*b[75];
    buffer[4] = a0[i]*b[4] + a1[i]*b[13] + a2[i]*b[22] +
                a3[i]*b[31] + a4[i]*b[40] + a5[i]*b[49] +
                a6[i]*b[58] + a7[i]*b[67] + a8[i]*b[76];
    buffer[5] = a0[i]*b[5] + a1[i]*b[14] + a2[i]*b[23] +
                a3[i]*b[32] + a4[i]*b[41] + a5[i]*b[50] +
                a6[i]*b[59] + a7[i]*b[68] + a8[i]*b[77];
    buffer[6] = a0[i]*b[6] + a1[i]*b[15] + a2[i]*b[24] +
                a3[i]*b[33] + a4[i]*b[42] + a5[i]*b[51] +
                a6[i]*b[60] + a7[i]*b[69] + a8[i]*b[78];
    buffer[7] = a0[i]*b[7] + a1[i]*b[16] + a2[i]*b[25] +
                a3[i]*b[34] + a4[i]*b[43] + a5[i]*b[52] +
                a6[i]*b[61] + a7[i]*b[70] + a8[i]*b[79];
    buffer[8] = a0[i]*b[8] + a1[i]*b[17] + a2[i]*b[26] +
                a3[i]*b[35] + a4[i]*b[44] + a5[i]*b[53] +
                a6[i]*b[62] + a7[i]*b[71] + a8[i]*b[80];
    for (unsigned int j = 0; j < _3N; ++j) {
      c_[j] = buffer[0]*a0[j] + buffer[1]*a1[j] + buffer[2]*a2[j] +
              buffer[3]*a3[j] + buffer[4]*a4[j] + buffer[5]*a5[j] +
              buffer[6]*a6[j] + buffer[7]*a7[j] + buffer[8]*a8[j];
    }
    c_ += _3N;
  }
}

template <typename T>
void matadd(
  const T* const a, const T* const b, const unsigned int size, T* const c) {
  for (unsigned int i = 0; i < size; ++i) {
    c[i] = a[i] + b[i];
  }
}

template <typename T>
void matinc_mul(const T* const a, const unsigned int size, T* const b, T k) {
  for (unsigned int i = 0; i < size; ++i) {
    b[i] += a[i] * k;
  }
}

template <typename T>
void matinc(const T* const a, const unsigned int size, T* const b) {
  for (unsigned int i = 0; i < size; ++i) {
    b[i] += a[i];
  }
}

template <typename T>
void matdec(const T* const a, const unsigned int size, T* const b) {
  for (unsigned int i = 0; i < size; ++i) {
    b[i] -= a[i];
  }
}

FORM_REGISTER_MATRIX_OP(float)
FORM_REGISTER_MATRIX_OP(double)
