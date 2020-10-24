#include <math.h>
#include <glog/logging.h>
#include <common/gauss_legendre.h>
#include "common.h"

template <typename T, unsigned int N>
void test_gauss_interp(bool layout=true, double tol=1e-20) {
  // init roots
  auto roots = (T*)malloc(N*sizeof(T));
  // init weigths
  auto weights = (T*)malloc(N*sizeof(T));
  // execute
  gauss_interp<T, N>(roots, weights, tol);
  if (layout) {
    LOG(INFO) << "matrix roots layout";
    print_mat<T>(roots, N, 1);
    LOG(INFO) << "matrix weights layout";
    print_mat<T>(weights, N, 1);
  }
  // free
  free(roots);
  free(weights);
}

int main(int argc, char* argv[]) {
  // log init
  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr = 1;
  // double precision tests
  test_gauss_interp<double, 2>();
  test_gauss_interp<double, 3>();
  test_gauss_interp<double, 4>();
  test_gauss_interp<double, 5>();
  test_gauss_interp<double, 6>();
  test_gauss_interp<double, 7>();
  test_gauss_interp<double, 8>();
  test_gauss_interp<double, 9>();
  test_gauss_interp<double, 10>();
  test_gauss_interp<double, 11>();
  test_gauss_interp<double, 12>();
  test_gauss_interp<double, 13>();
  test_gauss_interp<double, 14>();
  test_gauss_interp<double, 15>();
  test_gauss_interp<double, 16>();
  // single precision tests
  test_gauss_interp<float, 2>();
  test_gauss_interp<float, 3>();
  test_gauss_interp<float, 4>();
  test_gauss_interp<float, 5>();
  test_gauss_interp<float, 6>();
  test_gauss_interp<float, 7>();
  test_gauss_interp<float, 8>();
  test_gauss_interp<float, 9>();
  test_gauss_interp<float, 10>();
  test_gauss_interp<float, 11>();
  test_gauss_interp<float, 12>();
  test_gauss_interp<float, 13>();
  test_gauss_interp<float, 14>();
  test_gauss_interp<float, 15>();
  test_gauss_interp<float, 16>();
  return 0;
}
