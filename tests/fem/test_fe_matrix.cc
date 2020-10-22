#include <math.h>
#include <glog/logging.h>
#include <fem/fe_matrix.h>
#include "../common.h"

void test_matmul_n333_dp_cpu() {
  // init a
  unsigned int N = 8;
  unsigned int nEntryA = N * 3;
  unsigned int nBytesA = nEntryA * sizeof(double);
  auto a = (double*)malloc(nBytesA);
  rand_init<double>(a, nEntryA);
  LOG(INFO) << "matrix a layout";
  print_mat<double>(a, N, 3);
  // init b
  unsigned int nEntryB = 3 * 3;
  unsigned int nBytesB = nEntryB * sizeof(double);
  auto b = (double*)malloc(nBytesB);
  rand_init<double>(b, nEntryB);
  LOG(INFO) << "matrix b layout";
  print_mat<double>(b, 3, 3);
  // init c
  unsigned int nEntryC = N * 3;
  unsigned int nBytesC = nEntryC * sizeof(double);
  auto c = (double*)malloc(nBytesC);
  // execute
  fem::matmul_n333<double>(a, b, N, c);
  LOG(INFO) << "matrix c layout";
  print_mat<double>(c, N, 3);
  // free a
  free(a);
  // free b
  free(b);
  // free c
  free(c);
  LOG(INFO) << "test_matmul_n333_dp_cpu succeed";
}

void test_matmul_3nn3_dp_cpu() {
  // init a
  unsigned int N = 8;
  unsigned int nEntryA = 3 * N;
  unsigned int nBytesA = nEntryA * sizeof(double);
  auto a = (double*)malloc(nBytesA);
  rand_init<double>(a, nEntryA);
  LOG(INFO) << "matrix a layout";
  print_mat<double>(a, 3, N);
  // init b
  unsigned int nEntryB = N * 3;
  unsigned int nBytesB = nEntryB * sizeof(double);
  auto b = (double*)malloc(nBytesB);
  rand_init<double>(b, nEntryB);
  LOG(INFO) << "matrix b layout";
  print_mat<double>(b, N, 3);
  // init c
  unsigned int nEntryC = 3 * 3;
  unsigned int nBytesC = nEntryC * sizeof(double);
  auto c = (double*)malloc(nBytesC);
  // execute
  fem::matmul_3nn3<double>(a, b, N, c);
  LOG(INFO) << "matrix c layout";
  print_mat<double>(c, 3, 3);
  // free a
  free(a);
  // free b
  free(b);
  // free c
  free(c);
  LOG(INFO) << "test_matmul_3nn3_dp_cpu succeed";
}

void test_inv_33_dp_cpu() {
  // init a
  unsigned int nEntryA = 3 * 3;
  unsigned int nBytesA = nEntryA * sizeof(double);
  auto a = (double*)malloc(nBytesA);
  rand_init<double>(a, nEntryA);
  LOG(INFO) << "matrix a layout";
  print_mat<double>(a, 3, 3);
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
  LOG(INFO) << "matrix det layout";
  print_mat<double>(det, 1, 1);
  fem::inv_33(a, det, inv);
  LOG(INFO) << "matrix inv layout";
  print_mat<double>(inv, 3, 3);
  // free a
  free(a);
  // free det
  free(det);
  // free inv
  free(inv);
  LOG(INFO) << "test_inv_33_cpu_succeed";
}

int main(int argc, char* argv[]) {
  // log init
  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr = 1;
  // tests
  test_matmul_n333_dp_cpu();
  test_matmul_3nn3_dp_cpu();
  test_inv_33_dp_cpu();
  return 0;
}
