#include <glog/logging.h>
#include <fem/fe_matrix.h>
#include "../common.h"

void test_matmul_n333_dp_cpu() {
  // init data a
  unsigned int N = 8;
  unsigned int nEntryA = N * 3;
  unsigned int nBytesA = nEntryA * sizeof(double);
  auto a = (double*)malloc(nBytesA);
  init_rand_data<double>(a, nEntryA);
  printf("matrix a layout\n");
  print_mat<double>(a, N, 3);
  // init data b
  unsigned int nEntryB = 3 * 3;
  unsigned int nBytesB = nEntryB * sizeof(double);
  auto b = (double*)malloc(nBytesB);
  init_rand_data<double>(b, nEntryB);
  printf("matrix b layout\n");
  print_mat<double>(b, 3, 3);
  // init data c
  unsigned int nEntryC = N * 3;
  unsigned int nBytesC = nEntryC * sizeof(double);
  auto c = (double*)malloc(nBytesC);
  // execute
  fem::matmul_n333<double>(a, b, N, c);
  printf("matrix c layout\n");
  print_mat<double>(c, N, 3);
  // free data a
  free(a);
  // free data b
  free(b);
  // free data c
  free(c);
}

int main(int argc, char* argv[]) {
  test_matmul_n333_dp_cpu();
  return 0;
}
