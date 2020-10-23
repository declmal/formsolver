#include <glog/logging.h>
#include <fem/elem_stiff.h>
#include "../common.h"

void test_lin_trans_mat_tl_dp_cpu() {
  // init h0
  unsigned int N = 8;
  unsigned int nEntryH0 = N * 3;
  unsigned int nBytesH0 = nEntryH0 * sizeof(double);
  auto h0 = (double*)malloc(nBytesH0);
  rand_init<double>(h0, nEntryH0);
  LOG(INFO) << "matrix h0 layout";
  print_mat<double>(h0, N, 3);
  // init u0t
  unsigned int nEntryU0t = 3 * 3;
  unsigned int nBytesU0t = nEntryU0t * sizeof(double);
  auto u0t = (double*)malloc(nBytesU0t);
  rand_init<double>(u0t, nEntryU0t);
  LOG(INFO) << "matrix u0t layout";
  print_mat<double>(u0t, 3, 3);
  // init B0t_L
  unsigned int nEntryB0t_L = 6 * 3 * N;
  unsigned int nBytesB0t_L = nEntryB0t_L * sizeof(double);
  auto B0t_L = (double*)malloc(nBytesB0t_L);
  // execute
  fem::lin_trans_mat_tl(h0, u0t, N, B0t_L);
  LOG(INFO) << "matrix B0t_L layout";
  print_mat<double>(B0t_L, 6, 3*N);
  // free h0
  free(h0);
  // free u0t
  free(u0t);
  // free B0t_L
  free(B0t_L);
  LOG(INFO) << "test_lin_trans_mat_tl_dp_cpu succeed";
}

void test_nonlin_trans_mat_tl_dp_cpu() {
  // init h0
  unsigned int N = 8;
  unsigned int nEntryH0 = N * 3;
  unsigned int nBytesH0 = nEntryH0 * sizeof(double);
  auto h0 = (double*)malloc(nBytesH0);
  rand_init<double>(h0, nEntryH0);
  LOG(INFO) << "matrix h0 layout";
  print_mat<double>(h0, N, 3);
  // init B0_NL
  unsigned int nEntryB0_NL = 9 * 3 * N;
  unsigned int nBytesB0_NL = nEntryB0_NL * sizeof(double);
  auto B0_NL = (double*)malloc(nBytesB0_NL);
  // execute
  fem::nonlin_trans_mat_tl(h0, N, B0_NL);
  LOG(INFO) << "matrix B0t_L layout";
  print_mat<double>(B0_NL, 9, 3*N);
  // free h0
  free(h0);
  // free B0_NL
  free(B0_NL);
  LOG(INFO) << "test_nonlin_trans_mat_tl_dp_cpu succeed";
}

int main(int argc, char* argv[]) {
  // log init
  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr = 1;
  // tests
  test_lin_trans_mat_tl_dp_cpu();
  test_nonlin_trans_mat_tl_dp_cpu();
}
