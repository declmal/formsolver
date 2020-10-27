#include <glog/logging.h>
#include <fem/element/brick.h>
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
  typename T,
  template <typename> class BrickIPropType
>
void test_brick_interp_prop(bool layout=true) {
  BrickIPropType<T>::initialize();
  auto hbuf = BrickIPropType<T>::get_hbuf();
  auto weights = BrickIPropType<T>::get_weights();
  auto num_ipoints = BrickIPropType<T>::get_num_ipoints();
  auto num_nodes = BrickIPropType<T>::get_num_nodes();
  auto ndim = BrickIPropType<T>::get_ndim();
  auto stride_h = num_nodes * ndim;
  auto h = hbuf;
  if (layout) {
    LOG(INFO) << "tensor hbuf layout and weights";
    for (unsigned i = 0; i < num_ipoints; ++i) {
      LOG(INFO) << "tensor hbuf[" << i << "] layout";
      print_mat<T>(h, num_nodes, ndim);
      LOG(INFO) << "weights[" << i << "]: " << weights[i];
      h += stride_h;
    }
  }
  LOG(INFO) << "test_brick_interp_prop succeed, dtype: "
    << typeid(T).name() << " BrickIPropType: " 
    << typeid(BrickIPropType<T>).name();
}

template <
  typename T,
  template <typename> class BrickIPropType,
  template <typename> class BrickTLFormType
>
void test_brick_tl_form(bool layout=true) {
  auto Dim = BrickIPropType<T>::get_ndim();
  auto N = BrickIPropType<T>::get_num_nodes();
  // init X0
  auto nEntryX0 = Dim * N;
  auto X0 = (T*)malloc(nEntryX0*sizeof(T));
  init_rand(X0, nEntryX0);
  // init hbuf
  auto hbuf = BrickIPropType<T>::get_hbuf();
  // init Ut
  auto nEntryUt = Dim * N;
  auto Ut = (T*)malloc(nEntryUt*sizeof(T));
  init_rand<T>(Ut, nEntryUt);
  // init C0
  auto nRowC0 = 3*Dim - 3;
  auto nEntryC0 = nRowC0 * nRowC0;
  auto C0 = (T*)malloc(nEntryC0*sizeof(T));
  init_rand<T>(C0, nEntryC0);
  // init S0t
  auto nEntryS0t = Dim * Dim;
  auto S0t = (T*)malloc(nEntryS0t*sizeof(T));
  init_rand<T>(S0t, nEntryS0t);
  // init J0
  auto nEntryJ0 = Dim * Dim;
  auto J0 = (T*)malloc(nEntryJ0*sizeof(T));
  // init invJ0
  auto nEntryInvJ0 = Dim * Dim;
  auto invJ0 = (T*)malloc(nEntryInvJ0*sizeof(T));
  // init h0
  auto nEntryH0 = N * Dim;
  auto h0 = (T*)malloc(nEntryH0*sizeof(T));
  // init u0t
  auto nEntryU0t = Dim * Dim;
  auto u0t = (T*)malloc(nEntryU0t*sizeof(T));
  // init B0tL
  auto nRowB0tL = 3*Dim - 3;
  auto nColB0tL = Dim * N;
  auto nEntryB0tL = nRowB0tL * nColB0tL;
  auto B0tL = (T*)malloc(nEntryB0tL*sizeof(T));
  // init buf
  auto nEntryBuf = Dim * Dim;
  auto buf = (T*)malloc(nEntryBuf*sizeof(T));
  // init tmpK
  auto nRowTmpK = Dim * N;
  auto nEntryTmpK = nRowTmpK * nRowTmpK;
  auto tmpK = (T*)malloc(nEntryTmpK*sizeof(T));
  // init Ke
  auto nRowKe = Dim * N;
  auto nEntryKe = nRowKe * nRowKe;
  auto Ke = (T*)malloc(nEntryKe*sizeof(T));
  // init 
  BrickTLFormType<T>::form_elem_stiff(
    X0, hbuf, Ut, C0, S0t, J0, invJ0, h0, u0t, B0tL, buf, tmpK, Ke);
  // execute
  // free
  free(X0);
  free(J0);
  free(Ut);
  free(C0);
  free(S0t);
  free(invJ0);
  free(h0);
  free(u0t);
  free(B0tL);
  free(buf);
  free(tmpK);
  free(Ke);
  LOG(INFO) << "test_brick_tl_form passed, T: " << typeid(T).name()
    << ", BrickIPropType: " << typeid(BrickIPropType<T>).name()
    << ", BrickTLFormType: " << typeid(BrickTLFormType<T>).name();
}

int main(int argc, char* argv[]) {
  // log init
  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr = 1;
  bool layout = false;
  // double precison tests
  test_brick_interp_sum<double,8>(layout);
  test_brick_interp_deriv<double,8>(layout);
  test_brick_interp_sum<double, 20>(layout);
  test_brick_interp_prop<double,fem::C3D8IProp>(layout);
  // test_brick_interp_prop<double,fem::C3D8RIProp>(layout);
  test_brick_interp_prop<double,fem::C3D20IProp>(layout);
  test_brick_interp_prop<double,fem::C3D20RIProp>(layout);
  test_brick_tl_form<double,fem::C3D8IProp,fem::C3D8TLForm>();
  // test_brick_tl_form<double,fem::C3D8RIProp,fem::C3D8RTLForm>();
  test_brick_tl_form<double,fem::C3D20IProp,fem::C3D20TLForm>();
  test_brick_tl_form<double,fem::C3D20RIProp,fem::C3D20RTLForm>();
  LOG(INFO) << "double precision test passed";
  // single precision tests
  test_brick_interp_sum<float,8>(layout);
  test_brick_interp_deriv<float,8>(layout);
  test_brick_interp_sum<float,20>(layout);
  test_brick_interp_prop<float,fem::C3D8IProp>(layout);
  // test_brick_interp_prop<float,fem::C3D8RIProp>(layout);
  test_brick_interp_prop<float,fem::C3D20IProp>(layout);
  test_brick_interp_prop<float,fem::C3D20RIProp>(layout);
  test_brick_tl_form<float,fem::C3D8IProp,fem::C3D8TLForm>();
  // test_brick_tl_form<float,fem::C3D8RIProp,fem::C3D8RTLForm>();
  test_brick_tl_form<float,fem::C3D20IProp,fem::C3D20TLForm>();
  test_brick_tl_form<float,fem::C3D20RIProp,fem::C3D20RTLForm>();
  LOG(INFO) << "single precision test passed";
  return 0;
}
