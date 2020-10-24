#include <glog/logging.h>
#include <fem/elem_stiff.h>
#include "../common/common.h"

template <typename T>
void test_lin_trans_mat_tl_cpu(bool layout=true) {
  // init h0
  unsigned int N = 8;
  unsigned int nEntryH0 = N * 3;
  unsigned int nBytesH0 = nEntryH0 * sizeof(T);
  auto h0 = (T*)malloc(nBytesH0);
  init_rand<T>(h0, nEntryH0);
  // init u0t
  unsigned int nEntryU0t = 3 * 3;
  unsigned int nBytesU0t = nEntryU0t * sizeof(T);
  auto u0t = (T*)malloc(nBytesU0t);
  init_rand<T>(u0t, nEntryU0t);
  // init B0t_L
  unsigned int nEntryB0t_L = 6 * 3 * N;
  unsigned int nBytesB0t_L = nEntryB0t_L * sizeof(T);
  auto B0t_L = (T*)malloc(nBytesB0t_L);
  // execute
  fem::lin_trans_mat_tl(h0, u0t, N, B0t_L);
  if (layout) {
    LOG(INFO) << "matrix h0 layout";
    print_mat<T>(h0, N, 3);
    LOG(INFO) << "matrix u0t layout";
    print_mat<T>(u0t, 3, 3);
    LOG(INFO) << "matrix B0t_L layout";
    print_mat<T>(B0t_L, 6, 3*N);
  }
  // free
  free(h0);
  free(u0t);
  free(B0t_L);
  LOG(INFO) << "test_lin_trans_mat_tl_cpu succeed";
}

template <typename T>
void test_nonlin_trans_mat_tl_cpu(bool layout=true) {
  // init h0
  unsigned int N = 8;
  unsigned int nEntryH0 = N * 3;
  unsigned int nBytesH0 = nEntryH0 * sizeof(T);
  auto h0 = (T*)malloc(nBytesH0);
  init_rand<T>(h0, nEntryH0);
  // init B0_NL
  unsigned int nEntryB0_NL = 9 * 3 * N;
  unsigned int nBytesB0_NL = nEntryB0_NL * sizeof(T);
  auto B0_NL = (T*)malloc(nBytesB0_NL);
  // execute
  fem::nonlin_trans_mat_tl(h0, N, B0_NL);
  if (layout) {
    LOG(INFO) << "matrix h0 layout";
    print_mat<T>(h0, N, 3);
    LOG(INFO) << "matrix B0t_L layout";
    print_mat<T>(B0_NL, 9, 3*N);
  }
  // free
  free(h0);
  free(B0_NL);
  LOG(INFO) << "test_nonlin_trans_mat_tl_cpu succeed";
}

int main(int argc, char* argv[]) {
  // log init
  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr = 1;
  // double precision tests
  test_lin_trans_mat_tl_cpu<double>(false);
  test_nonlin_trans_mat_tl_cpu<double>(false);
  LOG(INFO) << "double precision test passed";
  // float precision tests
  test_lin_trans_mat_tl_cpu<float>(false);
  test_nonlin_trans_mat_tl_cpu<float>(false);
  LOG(INFO) << "single precision test passed";
}
