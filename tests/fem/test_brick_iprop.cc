#include <glog/logging.h>
#include <fem/iprop/brick_iprop.h>
#include <common/common.h>

template <typename T, unsigned int I, unsigned int N>
struct BrickInterpSum {
  static inline T interp_sum(T r0, T r1, T r2) {
    return
      BrickInterpSum<T,I-1,N>::interp_sum(r0, r1, r2) +
      fem::BrickInterpPoly<T,I,N>::compute(r0, r1, r2);
  }
};
template <typename T, unsigned int N>
struct BrickInterpSum<T, 0, N> {
  static inline T interp_sum(T r0, T r1, T r2) {
    return fem::BrickInterpPoly<T,0,N>::compute(r0, r1, r2);
  }
};

template <typename T, unsigned int N>
T brick_interp_sum(T r0, T r1, T r2) {
  return BrickInterpSum<T,N-1,N>::interp_sum(r0, r1, r2);
}

template <typename T, unsigned int N>
void test_brick_interp_sum(bool layout=false, double tol=1e-6) {
  T r0 = rand_gen<T>(-1, 1);
  T r1 = rand_gen<T>(-1, 1);
  T r2 = rand_gen<T>(-1, 1);
  auto res = brick_interp_sum<T,N>(r0, r1, r2);
  bool flag = (double)(res-1) <= tol;
  if (layout) {
    LOG(INFO) << "natural coordinates, r0: " 
      << r0 << ", r1: " << r1 << ", r2: " << r2;
    LOG(INFO) << "aggregation of interploation: " << res;
  }
  if (flag) {
    LOG(INFO) << "test_brick_interp_sum succeed, dtype: " 
      << typeid(T).name() << ", N: " << N;
  } else {
    LOG(FATAL) << "test_brick_interp_sum fail, dtype: "
      << typeid(T).name() << ", N: " << N;
  }
}
template <typename T, unsigned int N>
struct BrickInterpDerivValidate {
  static inline void interp_deriv(T* const h, T r0, T r1, T r2);
};
template <typename T>
struct BrickInterpDerivValidate<T,8> {
  static inline void interp_deriv(T* const h, T r0, T r1, T r2) {
    auto a0 = (T)1 + r0;
    auto m0 = (T)1 - r0;
    auto a1 = (T)1 + r1;
    auto m1 = (T)1 - r1;
    auto a2 = (T)1 + r2;
    auto m2 = (T)1 - r2;
    auto v = (T)0.125 * a1 * a2;
    h[0] = v;
    h[3] = -v;
    v = (T)0.125 * m1 * a2;
    h[6] = -v;
    h[9] = v;
    v = (T)0.125 * a1 * m2;
    h[12] = v;
    h[15] = -v;
    v = (T)0.125 * m1 * m2;
    h[18] = -v;
    h[21] = v;
    v = (T)0.125 * a0 * a2;
    h[1] = v;
    h[10] = -v;
    v = (T)0.125 * m0 * a2;
    h[4] = v;
    h[7] = -v;
    v = (T)0.125 * a0 * m2;
    h[13] = v;
    h[22] = -v;
    v = (T)0.125 * m0 * m2;
    h[16] = v;
    h[19] = -v;
    v = (T)0.125 * a0 * a1;
    h[2] = v;
    h[14] = -v;
    v = (T)0.125 * m0 * a1;
    h[5] = v;
    h[17] = -v;
    v = (T)0.125 * m0 * m1;
    h[8] = v;
    h[20] = -v;
    v = (T)0.125 * a0 * m1;
    h[11] = v;
    h[23] = -v;
  }
};

template <typename T, unsigned int N>
void test_brick_interp_deriv(bool layout=false, double tol=1e-6) {
  T r0 = rand_gen<T>(-1, 1);
  T r1 = rand_gen<T>(-1, 1);
  T r2 = rand_gen<T>(-1, 1);
  // init h0
  unsigned int nEntryH1 = N * 3;
  auto h0 = (T*)malloc(nEntryH1*sizeof(T));
  // init h1
  unsigned int nEntryH2 = N * 3;
  auto h1 = (T*)malloc(nEntryH2*sizeof(T));
  // execute
  fem::brick_interp_deriv<T,N>(h0, r0, r1, r2);
  // validate
  BrickInterpDerivValidate<T,N>::interp_deriv(h1, r0, r1, r2);
  if (layout) {
    LOG(INFO) << "natural coordinates, r0: " 
      << r0 << ", r1: " << r1 << ", r2: " << r2;
    LOG(INFO) << "matrix h0 layout";
    print_mat<T>(h0, N, 3);
    LOG(INFO) << "matrix h1 layout";
    print_mat<T>(h1, N, 3);
  }
  bool flag = validate<T>(h0, h1, N*3, tol);
  // free
  free(h0);
  free(h1);
  if (flag) {
    LOG(INFO) << "test_brick_interp_deriv succeed, dtype: " 
      << typeid(T).name() << ", N: " << N;
  } else {
    LOG(FATAL) << "test_brick_interp_deriv fail, dtype: "
      << typeid(T).name() << ", N: " << N;
  }
}

template <
  typename T, unsigned int N0, unsigned int N1, unsigned int N2, unsigned int N> 
void test_brick_interp_prop(bool layout=true) {
  fem::BrickIProp<T,N0,N1,N2,N> prop;
  if (layout) {
    LOG(INFO) << "tensor buf.h layout";
    auto NI = N0 * N1 * N2;
    T* h = prop.get_hbuf();
    for (unsigned i = 0; i < NI; ++i) {
      LOG(INFO) << "tensor buf.h[" << i << "] layout";
      print_mat<T>(h, N, 3);
      h += NI;
    }
  }
  LOG(INFO) << "test_brick_interp_prop succeed, dtype: "
    << typeid(T).name();
}

int main(int argc, char* argv[]) {
  // log init
  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr = 1;
  // double precison tests
  test_brick_interp_sum<double,8>();
  test_brick_interp_deriv<double,8>();
  test_brick_interp_sum<double, 20>();
  test_brick_interp_prop<double,2,2,2,8>();
  test_brick_interp_prop<double,1,1,1,8>();
  test_brick_interp_prop<double,3,3,3,20>();
  test_brick_interp_prop<double,2,2,2,20>();
  LOG(INFO) << "double precision test passed";
  // single precision tests
  test_brick_interp_sum<float,8>();
  test_brick_interp_deriv<float,8>();
  test_brick_interp_sum<float,20>();
  test_brick_interp_prop<float,2,2,2,8>();
  test_brick_interp_prop<float,1,1,1,8>();
  test_brick_interp_prop<float,3,3,3,20>();
  test_brick_interp_prop<float,2,2,2,20>();
  LOG(INFO) << "single precision test passed";
  return 0;
}
