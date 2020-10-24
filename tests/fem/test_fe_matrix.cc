#include <math.h>
#include <cblas.h>
#include <glog/logging.h>
#include <fem/fe_matrix.h>
#include "../common/common.h"

template <typename T>
void test_matmul_cblas(bool layout=true) {
  // init a
  unsigned int ra = 3;
  unsigned int ca = 2;
  unsigned int nEntryA = ra * ca;
  auto a = (T*)malloc(nEntryA*sizeof(T));
  init_arange<T>(a, nEntryA);
  // init b
  unsigned int rb = 2;
  unsigned int cb = 2;
  unsigned int nEntryB = rb * cb;
  auto b = (T*)malloc(nEntryB*sizeof(T));
  init_arange<T>(b, nEntryB);
  // init c
  unsigned int rc = ra;
  unsigned int cc = cb;
  unsigned int nEntryC = rc * cc;
  auto c = (T*)malloc(nEntryC*sizeof(T));
  // execute
  Matmul<T>::impl(
      CblasRowMajor, CblasNoTrans, CblasNoTrans,
      rc, cc, ca, (T)1.0, a, ca, b, cb, (T)0.0, c, cc);
  if (layout) {
    LOG(INFO) << "matrix a layout";
    print_mat<T>(a, ra, ca);
    LOG(INFO) << "matrix b layout";
    print_mat<T>(b, rb, cb);
    LOG(INFO) << "matrix c layout";
    print_mat<T>(c, rc, cc);
  }
  // free
  free(a);
  free(b);
  free(c);
}

template <typename T>
void test_matmul_n333_cpu(bool layout=false, double tol=1e-6) {
  // init a
  unsigned int N = 8;
  unsigned int nEntryA = N * 3;
  unsigned int nBytesA = nEntryA * sizeof(T);
  auto a = (T*)malloc(nBytesA);
  init_rand<T>(a, nEntryA);
  // init b
  unsigned int nEntryB = 3 * 3;
  unsigned int nBytesB = nEntryB * sizeof(T);
  auto b = (T*)malloc(nBytesB);
  init_rand<T>(b, nEntryB);
  // init c
  unsigned int nEntryC = N * 3;
  unsigned int nBytesC = nEntryC * sizeof(T);
  auto c = (T*)malloc(nBytesC);
  // init d
  unsigned int nEntryD = N * 3;
  unsigned int nBytesD = nEntryD * sizeof(T);
  auto d = (T*)malloc(nBytesD);
  // execute
  fem::matmul_n333<T>(a, b, N, c);
  // validate
  Matmul<T>::impl(
    CblasRowMajor, CblasNoTrans, CblasNoTrans,
    N, 3, 3, (T)1.0, a, 3, b, 3, (T)0.0, d, 3);
  if (layout) {
    LOG(INFO) << "matrix a layout";
    print_mat<T>(a, N, 3);
    LOG(INFO) << "matrix b layout";
    print_mat<T>(b, 3, 3);
    LOG(INFO) << "matrix c layout";
    print_mat<T>(c, N, 3);
    LOG(INFO) << "matrix d layout";
    print_mat<T>(d, N, 3);
  }
  bool flag = validate<T>(c, d, nEntryC, tol);
  // free
  free(a);
  free(b);
  free(c);
  free(d);
  if (flag) {
    LOG(INFO) << "test_matmul_n333_cpu succeed";
  } else {
    LOG(FATAL) << "test_matmul_n333_cpu failed";
  }
}

template <typename T>
void test_matmul_3nn3_cpu(bool layout=false, double tol=1e-6) {
  // init a
  unsigned int N = 8;
  unsigned int nEntryA = 3 * N;
  unsigned int nBytesA = nEntryA * sizeof(T);
  auto a = (T*)malloc(nBytesA);
  init_rand<T>(a, nEntryA);
  // init b
  unsigned int nEntryB = N * 3;
  unsigned int nBytesB = nEntryB * sizeof(T);
  auto b = (T*)malloc(nBytesB);
  init_rand<T>(b, nEntryB);
  // init c
  unsigned int nEntryC = 3 * 3;
  unsigned int nBytesC = nEntryC * sizeof(T);
  auto c = (T*)malloc(nBytesC);
  // init d
  unsigned int nEntryD = 3 * 3;
  unsigned int nBytesD = nEntryD * sizeof(T);
  auto d = (T*)malloc(nBytesD);
  // execute
  fem::matmul_3nn3<T>(a, b, N, c);
  // validate
  Matmul<T>::impl(
    CblasRowMajor, CblasNoTrans, CblasNoTrans,
    3, 3, N, (T)1.0, a, N, b, 3, (T)0.0, d, 3);
  if (layout) {
    LOG(INFO) << "matrix a layout";
    print_mat<T>(a, 3, N);
    LOG(INFO) << "matrix b layout";
    print_mat<T>(b, N, 3);
    LOG(INFO) << "matrix c layout";
    print_mat<T>(c, 3, 3);
    LOG(INFO) << "matrix d layout";
    print_mat<T>(d, 3, 3);
  }
  bool flag = validate<T>(c, d, nEntryC, tol);
  // free
  free(a);
  free(b);
  free(c);
  free(d);
  if (flag) {
    LOG(INFO) << "test_matmul_3nn3_cpu succeed";
  } else {
    LOG(FATAL) << "test_matmul_3nn3_cpu failed";
  }
}

template <typename T>
void test_inv_33_cpu(bool layout=false, double tol=1e-6) {
  // init a
  unsigned int nEntryA = 3 * 3;
  unsigned int nBytesA = nEntryA * sizeof(T);
  auto a = (T*)malloc(nBytesA);
  init_rand<T>(a, nEntryA);
  // init det
  unsigned int nBytesDet = 1 * sizeof(T);
  auto det = (T*)malloc(nBytesDet);
  // init inv
  unsigned int nEntryInv = 3 * 3;
  unsigned int nBytesInv = nEntryInv * sizeof(T);
  auto inv = (T*)malloc(nBytesInv);
  // init mul
  unsigned int nEntryMul = 3 * 3;
  unsigned int nBytesMul = nEntryMul * sizeof(T);
  auto mul = (T*)malloc(nBytesMul);
  // init unit
  unsigned int nEntryUnit = 3 * 3;
  unsigned int nBytesUnit = nEntryUnit * sizeof(T);
  auto unit = (T*)malloc(nBytesUnit);
  init_diag_unit<T>(unit, 3);
  // execute
  fem::det_33(a, det);
  if (abs(det[0]) < 1e-6) {
    free(a);
    free(det);
    free(inv);
    free(mul);
    free(unit);
    LOG(FATAL) << "singular matrix encountered, det: " << det[0];
  }
  fem::inv_33(a, det, inv);
  // validate
  Matmul<T>::impl(
    CblasRowMajor, CblasNoTrans, CblasNoTrans,
    3, 3, 3, (T)1.0, a, 3, inv, 3, (T)0.0, mul, 3);
  if (layout) {
    LOG(INFO) << "matrix a layout";
    print_mat<T>(a, 3, 3);
    LOG(INFO) << "matrix det layout";
    print_mat<T>(det, 1, 1);
    LOG(INFO) << "matrix inv layout";
    print_mat<T>(inv, 3, 3);
    LOG(INFO) << "matrix mul layout";
    print_mat<T>(mul, 3, 3);
  }
  bool flag = validate<T>(mul, unit, 9, tol);
  // free
  free(a);
  free(det);
  free(inv);
  free(mul);
  free(unit);
  if (flag) {
    LOG(INFO) << "test_inv_33_cpu succeed";
  } else {
    LOG(FATAL) << "test_inv_33_cpu failed";
  }
}

template <typename T>
void test_mattile_diag_33_cpu(bool layout=false, double tol=1e-6) {
  // init a
  unsigned int nEntryA = 3 * 3;
  unsigned int nBytesA = nEntryA * sizeof(T);
  auto a = (T*)malloc(nBytesA);
  init_rand<T>(a, nEntryA);
  // init tile
  unsigned int nEntryTile = 9 * 9;
  unsigned int nBytesTile = nEntryTile * sizeof(T);
  auto tile = (T*)malloc(nBytesTile);
  // init I matrices
  auto I0 = (T*)malloc(27*sizeof(T));
  init_zero<T>(I0, 27);
  I0[0] = I0[10] = I0[20] = 1.0;
  auto I1 = (T*)malloc(27*sizeof(T));
  init_zero<T>(I1, 27);
  I1[3] = I1[13] = I1[23] = 1.0;
  auto I2 = (T*)malloc(27*sizeof(T));
  init_zero<T>(I2, 27);
  I2[6] = I2[16] = I2[26] = 1.0;
  // init tmps
  auto tmp = (T*)malloc(27*sizeof(T));
  auto mul0 = (T*)malloc(81*sizeof(T));
  auto mul1 = (T*)malloc(81*sizeof(T));
  auto mul2 = (T*)malloc(81*sizeof(T));
  // init add
  auto add0 = (T*)malloc(81*sizeof(T));
  auto add1 = (T*)malloc(81*sizeof(T));
  // execute
  fem::mattile_diag_33<T>(a, tile);
  // validate
  Matmul<T>::impl(
    CblasRowMajor, CblasTrans, CblasNoTrans,
    9, 3, 3, (T)1.0, I0, 9, a, 3, (T)0.0, tmp, 3);
  Matmul<T>::impl(
    CblasRowMajor, CblasNoTrans, CblasNoTrans,
    9, 9, 3, (T)1.0, tmp, 3, I0, 9, (T)0.0, mul0, 9);
  Matmul<T>::impl(
    CblasRowMajor, CblasTrans, CblasNoTrans,
    9, 3, 3, (T)1.0, I1, 9, a, 3, (T)0.0, tmp, 3);
  Matmul<T>::impl(
    CblasRowMajor, CblasNoTrans, CblasNoTrans,
    9, 9, 3, (T)1.0, tmp, 3, I1, 9, (T)0.0, mul1, 9);
  Matmul<T>::impl(
    CblasRowMajor, CblasTrans, CblasNoTrans,
    9, 3, 3, (T)1.0, I2, 9, a, 3, (T)0.0, tmp, 3);
  Matmul<T>::impl(
    CblasRowMajor, CblasNoTrans, CblasNoTrans,
    9, 9, 3, (T)1.0, tmp, 3, I2, 9, (T)0.0, mul2, 9);
  matadd<T>(mul0, mul1, 81, add0);
  matadd<T>(add0, mul2, 81, add1);
  if (layout) {
    LOG(INFO) << "matrix a layout";
    print_mat<T>(a, 3, 3);
    LOG(INFO) << "matrix tile layout";
    print_mat<T>(tile, 9, 9);
    LOG(INFO) << "matrix add1 layout";
    print_mat<T>(add1, 9, 9);
  }
  bool flag = validate<T>(tile, add1, 81, tol);
  // free
  free(a);
  free(tile);
  free(I0);
  free(I1);
  free(I2);
  free(tmp);
  free(mul0);
  free(mul1);
  free(mul2);
  free(add0);
  free(add1);
  if (flag) {
    LOG(INFO) << "test_mattile_diag33_cpu succeed";
  } else {
    LOG(FATAL) << "test_mattile_diag33_cpu failed";
  }
}

template <typename T>
void test_matmul2_3n6_66_63n_cpu(bool layout=false, double tol=1e-6) {
  // init a
  unsigned int N = 8;
  unsigned int _3N = 3 * N;
  unsigned int nEntryA = 6 * _3N;
  unsigned int nBytesA = nEntryA * sizeof(T);
  auto a = (T*)malloc(nBytesA);
  init_rand<T>(a, nEntryA);
  // init b
  unsigned int nEntryB = 6 * 6;
  unsigned int nBytesB = nEntryB * sizeof(T);
  auto b = (T*)malloc(nBytesB);
  init_rand<T>(b, nEntryB);
  // init buffer
  unsigned int nEntryBuffer = 6;
  unsigned int nBytesBuffer = nEntryBuffer * sizeof(T);
  auto buffer = (T*)malloc(nBytesBuffer);
  // init c
  unsigned int nEntryC = _3N * _3N;
  unsigned int nBytesC = nEntryC * sizeof(T);
  auto c = (T*)malloc(nBytesC);
  // init d
  unsigned int nEntryD = _3N * 6;
  unsigned int nBytesD = nEntryD * sizeof(T);
  auto d = (T*)malloc(nBytesD);
  // init e
  unsigned int nEntryE = _3N * _3N;
  unsigned int nBytesE = nEntryE * sizeof(T);
  auto e = (T*)malloc(nBytesE);
  // execute
  fem::matmul2_3n6_66_63n(a, b, N, buffer, c);
  // validate
  Matmul<T>::impl(
    CblasRowMajor, CblasTrans, CblasNoTrans,
    _3N, 6, 6, (T)1.0, a, _3N, b, 6, (T)0.0, d, 6);
  Matmul<T>::impl(
    CblasRowMajor, CblasNoTrans, CblasNoTrans,
    _3N, _3N, 6, (T)1.0, d, 6, a, _3N, (T)0.0, e, _3N);
  if (layout) {
    LOG(INFO) << "matrix a layout";
    print_mat<T>(a, 6, _3N);
    LOG(INFO) << "matrix b layout";
    print_mat<T>(b, 6, 6);
    LOG(INFO) << "matrix c layout";
    print_mat<T>(c, _3N, _3N);
    LOG(INFO) << "matrix e layout";
    print_mat<T>(e, _3N, _3N);
  }
  bool flag = validate<T>(c, e, _3N*_3N, tol);
  // free
  free(a);
  free(b);
  free(buffer);
  free(c);
  free(d);
  free(e);
  if (flag) {
    LOG(INFO) << "test_matmul2_3n6_66_63n_cpu succeed";
  } else {
    LOG(FATAL) << "test_matmul2_3n6_66_63n_cpu failed";
  }
}

template <typename T>
void test_matmul2_3n9_99_93n_cpu(bool layout=false, double tol=1e-6) {
  // init a
  unsigned int N = 8;
  unsigned int _3N = 3 * N;
  unsigned int nEntryA = 9 * _3N;
  unsigned int nBytesA = nEntryA * sizeof(T);
  auto a = (T*)malloc(nBytesA);
  init_rand<T>(a, nEntryA);
  // init b
  unsigned int nEntryB = 9 * 9;
  unsigned int nBytesB = nEntryB * sizeof(T);
  auto b = (T*)malloc(nBytesB);
  init_rand<T>(b, nEntryB);
  // init buffer
  unsigned int nEntryBuffer = 9;
  unsigned int nBytesBuffer = nEntryBuffer * sizeof(T);
  auto buffer = (T*)malloc(nBytesBuffer);
  // init c
  unsigned int nEntryC = _3N * _3N;
  unsigned int nBytesC = nEntryC * sizeof(T);
  auto c = (T*)malloc(nBytesC);
  // init d
  unsigned int nEntryD = _3N * 9;
  unsigned int nBytesD = nEntryD * sizeof(T);
  auto d = (T*)malloc(nBytesD);
  // init e
  unsigned int nEntryE = _3N * _3N;
  unsigned int nBytesE = nEntryE * sizeof(T);
  auto e = (T*)malloc(nBytesE);
  // execute
  fem::matmul2_3n9_99_93n(a, b, N, buffer, c);
  // validate
  Matmul<T>::impl(
    CblasRowMajor, CblasTrans, CblasNoTrans,
    _3N, 9, 9, (T)1.0, a, _3N, b, 9, (T)0.0, d, 9);
  Matmul<T>::impl(
    CblasRowMajor, CblasNoTrans, CblasNoTrans,
    _3N, _3N, 9, (T)1.0, d, 9, a, _3N, (T)0.0, e, _3N);
  if (layout) {
    LOG(INFO) << "matrix a layout";
    print_mat<T>(a, 9, _3N);
    LOG(INFO) << "matrix b layout";
    print_mat<T>(b, 9, 9);
    LOG(INFO) << "matrix c layout";
    print_mat<T>(c, _3N, _3N);
    LOG(INFO) << "matrix e layout";
    print_mat<T>(e, _3N, _3N);
  }
  bool flag = validate<T>(c, e, _3N*_3N, tol);
  // free
  free(a);
  free(b);
  free(buffer);
  free(c);
  free(d);
  free(e);
  if (flag) {
    LOG(INFO) << "test_matmul2_3n9_99_93n_cpu succeed";
  } else {
    LOG(FATAL) << "test_matmul2_3n9_99_93n_cpu failed";
  }
}

int main(int argc, char* argv[]) {
  // log init
  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr = 1;
  // double precison tests
  test_matmul_cblas<double>(false);
  test_matmul_n333_cpu<double>();
  test_matmul_3nn3_cpu<double>();
  test_inv_33_cpu<double>();
  test_mattile_diag_33_cpu<double>();
  test_matmul2_3n6_66_63n_cpu<double>();
  test_matmul2_3n9_99_93n_cpu<double>();
  LOG(INFO) << "double precision test passed";
  // single precision tests
  test_matmul_cblas<float>(false);
  test_matmul_n333_cpu<float>();
  test_matmul_3nn3_cpu<float>();
  test_inv_33_cpu<float>();
  test_mattile_diag_33_cpu<float>();
  test_matmul2_3n6_66_63n_cpu<float>();
  test_matmul2_3n9_99_93n_cpu<float>();
  return 0;
}
