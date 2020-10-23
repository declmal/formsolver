#include <math.h>
#include <cblas.h>
#include <lapack.h>
#include <glog/logging.h>
#include <fem/fe_matrix.h>
#include "../common.h"

void test_matmul_cblas(bool layout=true) {
  // init a
  unsigned int ra = 3;
  unsigned int ca = 2;
  unsigned int nEntryA = ra * ca;
  auto a = (double*)malloc(nEntryA*sizeof(double));
  init_arange<double>(a, nEntryA);
  // init b
  unsigned int rb = 2;
  unsigned int cb = 2;
  unsigned int nEntryB = rb * cb;
  auto b = (double*)malloc(nEntryB*sizeof(double));
  init_arange<double>(b, nEntryB);
  // init c
  unsigned int rc = ra;
  unsigned int cc = cb;
  unsigned int nEntryC = rc * cc;
  auto c = (double*)malloc(nEntryC*sizeof(double));
  // execute
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rc, cc, ca, 1.0, a, ca, b, cb, 0.0, c, cc);
  if (layout) {
    LOG(INFO) << "matrix a layout";
    print_mat_dp(a, ra, ca);
    LOG(INFO) << "matrix b layout";
    print_mat_dp(b, rb, cb);
    LOG(INFO) << "matrix c layout";
    print_mat_dp(c, rc, cc);
  }
  // free
  free(a);
  free(b);
  free(c);
}

void test_inv_cblas(bool layout=true) {
  // init a
  auto a = (double*)malloc(9*sizeof(double));
  init_arange<double>(a, 9);
  // execute
  if (layout) {
    LOG(INFO) << "matrix a layout";
    print_mat_dp(a, 3, 3);
  }
  
  // free
  free(a);
}

void test_matmul_n333_dp_cpu(bool layout=false, double tol=1e-6) {
  // init a
  unsigned int N = 8;
  unsigned int nEntryA = N * 3;
  unsigned int nBytesA = nEntryA * sizeof(double);
  auto a = (double*)malloc(nBytesA);
  init_rand<double>(a, nEntryA);
  // init b
  unsigned int nEntryB = 3 * 3;
  unsigned int nBytesB = nEntryB * sizeof(double);
  auto b = (double*)malloc(nBytesB);
  init_rand<double>(b, nEntryB);
  // init c
  unsigned int nEntryC = N * 3;
  unsigned int nBytesC = nEntryC * sizeof(double);
  auto c = (double*)malloc(nBytesC);
  // init d
  unsigned int nEntryD = N * 3;
  unsigned int nBytesD = nEntryD * sizeof(double);
  auto d = (double*)malloc(nBytesD);
  // execute
  fem::matmul_n333<double>(a, b, N, c);
  // validate
  cblas_dgemm(
    CblasRowMajor, CblasNoTrans, CblasNoTrans,
    N, 3, 3, 1.0, a, 3, b, 3, 0.0, d, 3);
  if (layout) {
    LOG(INFO) << "matrix a layout";
    print_mat_dp(a, N, 3);
    LOG(INFO) << "matrix b layout";
    print_mat_dp(b, 3, 3);
    LOG(INFO) << "matrix c layout";
    print_mat_dp(c, N, 3);
    LOG(INFO) << "matrix d layout";
    print_mat_dp(d, N, 3);
  }
  if (validate_dp(c, d, nEntryC, tol)) {
    LOG(INFO) << "test_matmul_n333_dp_cpu succeed";
  } else {
    LOG(WARNING) << "test_matmul_n333_dp_cpu failed";
  }
  // free
  free(a);
  free(b);
  free(c);
  free(d);
}

void test_matmul_3nn3_dp_cpu(bool layout=false, double tol=1e-6) {
  // init a
  unsigned int N = 8;
  unsigned int nEntryA = 3 * N;
  unsigned int nBytesA = nEntryA * sizeof(double);
  auto a = (double*)malloc(nBytesA);
  init_rand<double>(a, nEntryA);
  // init b
  unsigned int nEntryB = N * 3;
  unsigned int nBytesB = nEntryB * sizeof(double);
  auto b = (double*)malloc(nBytesB);
  init_rand<double>(b, nEntryB);
  // init c
  unsigned int nEntryC = 3 * 3;
  unsigned int nBytesC = nEntryC * sizeof(double);
  auto c = (double*)malloc(nBytesC);
  // init d
  unsigned int nEntryD = 3 * 3;
  unsigned int nBytesD = nEntryD * sizeof(double);
  auto d = (double*)malloc(nBytesD);
  // execute
  fem::matmul_3nn3<double>(a, b, N, c);
  // validate
  cblas_dgemm(
    CblasRowMajor, CblasNoTrans, CblasNoTrans,
    3, 3, N, 1.0, a, N, b, 3, 0.0, d, 3);
  if (layout) {
    LOG(INFO) << "matrix a layout";
    print_mat_dp(a, 3, N);
    LOG(INFO) << "matrix b layout";
    print_mat_dp(b, N, 3);
    LOG(INFO) << "matrix c layout";
    print_mat_dp(c, 3, 3);
    LOG(INFO) << "matrix d layout";
    print_mat_dp(d, 3, 3);
  }
  if (validate_dp(c, d, nEntryC)) {
    LOG(INFO) << "test_matmul_3nn3_dp_cpu succeed";
  } else {
    LOG(WARNING) << "test_matmul_3nn3_dp_cpu failed";
  }
  // free
  free(a);
  free(b);
  free(c);
  free(d);
}

void test_inv_33_dp_cpu(bool layout=false, double tol=1e-6) {
  // init a
  unsigned int nEntryA = 3 * 3;
  unsigned int nBytesA = nEntryA * sizeof(double);
  auto a = (double*)malloc(nBytesA);
  init_rand<double>(a, nEntryA);
  // init det
  unsigned int nBytesDet = 1 * sizeof(double);
  auto det = (double*)malloc(nBytesDet);
  // init inv
  unsigned int nEntryInv = 3 * 3;
  unsigned int nBytesInv = nEntryInv * sizeof(double);
  auto inv = (double*)malloc(nBytesInv);
  // execute
  fem::det_33(a, det);
  if (abs(det[0]) < 1e-6) {
    free(a);
    free(det);
    free(inv);
    LOG(FATAL) << "singular matrix encountered, det: " << det[0];
  }
  fem::inv_33(a, det, inv);
  // validate
  if (layout) {
    LOG(INFO) << "matrix a layout";
    print_mat_dp(a, 3, 3);
    LOG(INFO) << "matrix det layout";
    print_mat_dp(det, 1, 1);
    LOG(INFO) << "matrix inv layout";
    print_mat_dp(inv, 3, 3);
  }
  // free
  free(a);
  free(det);
  free(inv);
  LOG(INFO) << "test_inv_33_cpu_succeed";
}

void test_mattile_diag_33_dp_cpu() {
  // init a
  unsigned int nEntryA = 3 * 3;
  unsigned int nBytesA = nEntryA * sizeof(double);
  auto a = (double*)malloc(nBytesA);
  init_rand<double>(a, nEntryA);
  LOG(INFO) << "matrix a layout";
  print_mat_dp(a, 3, 3);
  // init tile
  unsigned int nEntryTile = 9 * 9;
  unsigned int nBytesTile = nEntryTile * sizeof(double);
  auto tile = (double*)malloc(nBytesTile);
  // execute
  fem::mattile_diag_33(a, tile);
  LOG(INFO) << "matrix tile layout";
  print_mat_dp(tile, 9, 9);
  // free a
  free(a);
  // free tile
  free(tile);
  LOG(INFO) << "test_mattile_diag33_dp_cpu succeed";
}

void test_matmul2_3n6_66_63n_dp_cpu(bool layout=false, double tol=1e-6) {
  // init a
  unsigned int N = 8;
  unsigned int _3N = 3 * N;
  unsigned int nEntryA = 6 * _3N;
  unsigned int nBytesA = nEntryA * sizeof(double);
  auto a = (double*)malloc(nBytesA);
  init_rand<double>(a, nEntryA);
  // init b
  unsigned int nEntryB = 6 * 6;
  unsigned int nBytesB = nEntryB * sizeof(double);
  auto b = (double*)malloc(nBytesB);
  init_rand<double>(b, nEntryB);
  // init buffer
  unsigned int nEntryBuffer = 6;
  unsigned int nBytesBuffer = nEntryBuffer * sizeof(double);
  auto buffer = (double*)malloc(nBytesBuffer);
  // init c
  unsigned int nEntryC = _3N * _3N;
  unsigned int nBytesC = nEntryC * sizeof(double);
  auto c = (double*)malloc(nBytesC);
  // init d
  unsigned int nEntryD = _3N * 6;
  unsigned int nBytesD = nEntryD * sizeof(double);
  auto d = (double*)malloc(nBytesD);
  // init e
  unsigned int nEntryE = _3N * _3N;
  unsigned int nBytesE = nEntryE * sizeof(double);
  auto e = (double*)malloc(nBytesE);
  // execute
  fem::matmul2_3n6_66_63n(a, b, N, buffer, c);
  // validate
  cblas_dgemm(
    CblasRowMajor, CblasTrans, CblasNoTrans,
    _3N, 6, 6, 1.0, a, _3N, b, 6, 0.0, d, 6);
  cblas_dgemm(
    CblasRowMajor, CblasNoTrans, CblasNoTrans,
    _3N, _3N, 6, 1.0, d, 6, a, _3N, 0.0, e, _3N);
  if (layout) {
    LOG(INFO) << "matrix a layout";
    print_mat_dp(a, 6, _3N);
    LOG(INFO) << "matrix b layout";
    print_mat_dp(b, 6, 6);
    LOG(INFO) << "matrix c layout";
    print_mat_dp(c, _3N, _3N);
    LOG(INFO) << "matrix e layout";
    print_mat_dp(e, _3N, _3N);
  }
  if (validate_dp(c, e, _3N*_3N)) {
    LOG(INFO) << "test_matmul2_3n6_66_63n succeed";
  } else {
    LOG(WARNING) << "test_matmul2_3n6_66_63n failed";
  }
  // free
  free(a);
  free(b);
  free(buffer);
  free(c);
  free(d);
  free(e);
}

int main(int argc, char* argv[]) {
  // log init
  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr = 1;
  // tests
  // test_matmul_cblas();
  test_matmul_n333_dp_cpu();
  test_matmul_3nn3_dp_cpu();
  // test_inv_33_dp_cpu();
  // test_mattile_diag_33_dp_cpu();
  test_matmul2_3n6_66_63n_dp_cpu();
  return 0;
}
