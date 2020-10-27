#include <math.h>
#include <cblas.h>
#include <glog/logging.h>
#include <common/matrix.h>
#include <common/common.h>
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
  Matmul<T>::matmul(
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

template <typename T, unsigned int Dim, unsigned int N>
void test_matmul_nddd(bool layout=false, double tol=1e-6) {
  // init a
  unsigned int nEntryA = N * Dim;
  unsigned int nBytesA = nEntryA * sizeof(T);
  auto a = (T*)malloc(nBytesA);
  init_rand<T>(a, nEntryA);
  // init b
  unsigned int nEntryB = Dim * Dim;
  unsigned int nBytesB = nEntryB * sizeof(T);
  auto b = (T*)malloc(nBytesB);
  init_rand<T>(b, nEntryB);
  // init c
  unsigned int nEntryC = N * Dim;
  unsigned int nBytesC = nEntryC * sizeof(T);
  auto c = (T*)malloc(nBytesC);
  // init d
  unsigned int nEntryD = N * Dim;
  unsigned int nBytesD = nEntryD * sizeof(T);
  auto d = (T*)malloc(nBytesD);
  // execute
  FEMatrix<T,Dim>::matmul_nddd(a, b, N, c);
  // validate
  Matmul<T>::matmul(
    CblasRowMajor, CblasNoTrans, CblasNoTrans,
    N, Dim, Dim, (T)1.0, a, Dim, b, Dim, (T)0.0, d, Dim);
  if (layout) {
    LOG(INFO) << "matrix a layout";
    print_mat<T>(a, N, Dim);
    LOG(INFO) << "matrix b layout";
    print_mat<T>(b, Dim, Dim);
    LOG(INFO) << "matrix c layout";
    print_mat<T>(c, N, Dim);
    LOG(INFO) << "matrix d layout";
    print_mat<T>(d, N, Dim);
  }
  bool flag = validate<T>(c, d, nEntryC, tol);
  // free
  free(a);
  free(b);
  free(c);
  free(d);
  if (flag) {
    LOG(INFO) << "test_matmul_nddd succeed, T: " << typeid(T).name() <<
      ", Dim: " << Dim << ", N: " << N;
  } else {
    LOG(FATAL) << "test_matmul_nddd failed, T: " << typeid(T).name() <<
      ": Dim: " << Dim << ", N: " << N;
  }
}

template <typename T, unsigned int Dim, unsigned int N>
void test_matmul_dnnd(bool layout=false, double tol=1e-6) {
  // init a
  unsigned int nEntryA = Dim * N;
  unsigned int nBytesA = nEntryA * sizeof(T);
  auto a = (T*)malloc(nBytesA);
  init_rand<T>(a, nEntryA);
  // init b
  unsigned int nEntryB = N * Dim;
  unsigned int nBytesB = nEntryB * sizeof(T);
  auto b = (T*)malloc(nBytesB);
  init_rand<T>(b, nEntryB);
  // init c
  unsigned int nEntryC = Dim * Dim;
  unsigned int nBytesC = nEntryC * sizeof(T);
  auto c = (T*)malloc(nBytesC);
  // init d
  unsigned int nEntryD = Dim * Dim;
  unsigned int nBytesD = nEntryD * sizeof(T);
  auto d = (T*)malloc(nBytesD);
  // execute
  FEMatrix<T,Dim>::matmul_dnnd(a, b, N, c);
  // validate
  Matmul<T>::matmul(
    CblasRowMajor, CblasNoTrans, CblasNoTrans,
    Dim, Dim, N, (T)1.0, a, N, b, Dim, (T)0.0, d, Dim);
  if (layout) {
    LOG(INFO) << "matrix a layout";
    print_mat<T>(a, Dim, N);
    LOG(INFO) << "matrix b layout";
    print_mat<T>(b, N, Dim);
    LOG(INFO) << "matrix c layout";
    print_mat<T>(c, Dim, Dim);
    LOG(INFO) << "matrix d layout";
    print_mat<T>(d, Dim, Dim);
  }
  bool flag = validate<T>(c, d, nEntryC, tol);
  // free
  free(a);
  free(b);
  free(c);
  free(d);
  if (flag) {
    LOG(INFO) << "test_matmul_dnnd succeed, T: " << typeid(T).name() <<
      ", Dim: " << Dim << ", N: " << N;
  } else {
    LOG(FATAL) << "test_matmul_dnnd failed, T: " << typeid(T).name() <<
      ", Dim: " << Dim << ", N: " << N;
  }
}

template <typename T, unsigned int Dim>
void test_inv_dd(bool layout=false, double tol=1e-6) {
  // init a
  unsigned int nEntryA = Dim * Dim;
  unsigned int nBytesA = nEntryA * sizeof(T);
  auto a = (T*)malloc(nBytesA);
  init_rand<T>(a, nEntryA);
  // init inv
  unsigned int nEntryInv = Dim * Dim;
  unsigned int nBytesInv = nEntryInv * sizeof(T);
  auto inv = (T*)malloc(nBytesInv);
  // init mul
  unsigned int nEntryMul = Dim * Dim;
  unsigned int nBytesMul = nEntryMul * sizeof(T);
  auto mul = (T*)malloc(nBytesMul);
  // init unit
  unsigned int nEntryUnit = Dim * Dim;
  unsigned int nBytesUnit = nEntryUnit * sizeof(T);
  auto unit = (T*)malloc(nBytesUnit);
  init_diag_unit<T>(unit, Dim);
  // execute
  auto det = FEMatrix<T,Dim>::det_dd(a);
  if (abs(det) < 1e-6) {
    free(a);
    free(inv);
    free(mul);
    free(unit);
    LOG(FATAL) << "singular matrix encountered, det: " << det;
  }
  FEMatrix<T,Dim>::inv_dd(a, det, inv);
  // validate
  Matmul<T>::matmul(
    CblasRowMajor, CblasNoTrans, CblasNoTrans,
    Dim, Dim, Dim, (T)1.0, a, Dim, inv, Dim, (T)0.0, mul, Dim);
  if (layout) {
    LOG(INFO) << "matrix a layout";
    print_mat<T>(a, Dim, Dim);
    LOG(INFO) << "det: " << det;
    LOG(INFO) << "matrix inv layout";
    print_mat<T>(inv, Dim, Dim);
    LOG(INFO) << "matrix mul layout";
    print_mat<T>(mul, Dim, Dim);
  }
  bool flag = validate<T>(mul, unit, Dim*Dim, tol);
  // free
  free(a);
  free(inv);
  free(mul);
  free(unit);
  if (flag) {
    LOG(INFO) << "test_inv_dd succeed, T: " 
      << typeid(T).name() << ", Dim: " << Dim;
  } else {
    LOG(FATAL) << "test_inv_dd failed, T: " 
      << typeid(T).name() << ", Dim: " << Dim;
   }
}

template <typename T, unsigned int Dim>
struct DiagVal {
  static inline void diag_val(const T* const a, T* const out);
};
template <typename T>
struct DiagVal<T,3> {
  static inline void diag_val(const T* const a, T* const out) {
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
    // calcultate
    Matmul<T>::matmul(
      CblasRowMajor, CblasTrans, CblasNoTrans,
      9, 3, 3, (T)1.0, I0, 9, a, 3, (T)0.0, tmp, 3);
    Matmul<T>::matmul(
      CblasRowMajor, CblasNoTrans, CblasNoTrans,
      9, 9, 3, (T)1.0, tmp, 3, I0, 9, (T)0.0, mul0, 9);
    Matmul<T>::matmul(
      CblasRowMajor, CblasTrans, CblasNoTrans,
      9, 3, 3, (T)1.0, I1, 9, a, 3, (T)0.0, tmp, 3);
    Matmul<T>::matmul(
      CblasRowMajor, CblasNoTrans, CblasNoTrans,
      9, 9, 3, (T)1.0, tmp, 3, I1, 9, (T)0.0, mul1, 9);
    Matmul<T>::matmul(
      CblasRowMajor, CblasTrans, CblasNoTrans,
      9, 3, 3, (T)1.0, I2, 9, a, 3, (T)0.0, tmp, 3);
    Matmul<T>::matmul(
      CblasRowMajor, CblasNoTrans, CblasNoTrans,
      9, 9, 3, (T)1.0, tmp, 3, I2, 9, (T)0.0, mul2, 9);
    matadd<T>(mul0, mul1, 81, add0);
    matadd<T>(add0, mul2, 81, out);
    // free
    free(I0);
    free(I1);
    free(I2);
    free(tmp);
    free(mul0);
    free(mul1);
    free(mul2);
    free(add0);
  }
};

template <typename T, unsigned int Dim>
void test_mattile_diag_dd(bool layout=false, double tol=1e-6) {
  // init a
  unsigned int nEntryA = Dim * Dim;
  unsigned int nBytesA = nEntryA * sizeof(T);
  auto a = (T*)malloc(nBytesA);
  init_rand<T>(a, nEntryA);
  // init tile
  unsigned int Dim2 = Dim * Dim;
  unsigned int nEntryTile = Dim2 * Dim2;
  unsigned int nBytesTile = nEntryTile * sizeof(T);
  auto tile = (T*)malloc(nBytesTile);
  // init out
  auto out = (T*)malloc(nEntryTile*sizeof(T));
  // execute
  FEMatrix<T,Dim>::mattile_diag_dd(a, tile);
  // validate
  DiagVal<T,Dim>::diag_val(a, out);
  if (layout) {
    LOG(INFO) << "matrix a layout";
    print_mat<T>(a, Dim, Dim);
    LOG(INFO) << "matrix tile layout";
    print_mat<T>(tile, Dim2, Dim2);
    LOG(INFO) << "matrix out layout";
    print_mat<T>(out, Dim2, Dim2);
  }
  bool flag = validate<T>(tile, out, nEntryTile, tol);
  // free
  free(a);
  free(tile);
  free(out);
  if (flag) {
    LOG(INFO) << "test_mattile_diag33_cpu succeed, T: " 
      << typeid(T).name() << ", Dim: " << Dim;
  } else {
    LOG(FATAL) << "test_mattile_diag33_cpu failed, T: " 
      << typeid(T).name() << ", Dim: " << Dim;
  }
}

template <typename T, unsigned int Dim, unsigned int N>
void test_matmul2_dne_ee_edn(bool layout=false, double tol=1e-6) {
  // init a
  unsigned int _N = Dim * N;
  unsigned int _N2 = 3*Dim - 3;
  unsigned int nEntryA = _N2 * _N;
  unsigned int nBytesA = nEntryA * sizeof(T);
  auto a = (T*)malloc(nBytesA);
  init_rand<T>(a, nEntryA);
  // init b
  unsigned int nEntryB = _N2 * _N2;
  unsigned int nBytesB = nEntryB * sizeof(T);
  auto b = (T*)malloc(nBytesB);
  init_rand<T>(b, nEntryB);
  // init buffer
  unsigned int nEntryBuffer = _N2;
  unsigned int nBytesBuffer = nEntryBuffer * sizeof(T);
  auto buffer = (T*)malloc(nBytesBuffer);
  // init c
  unsigned int nEntryC = _N * _N;
  unsigned int nBytesC = nEntryC * sizeof(T);
  auto c = (T*)malloc(nBytesC);
  // init d
  unsigned int nEntryD = _N * _N2;
  unsigned int nBytesD = nEntryD * sizeof(T);
  auto d = (T*)malloc(nBytesD);
  // init e
  unsigned int nEntryE = _N * _N;
  unsigned int nBytesE = nEntryE * sizeof(T);
  auto e = (T*)malloc(nBytesE);
  // execute
  FEMatrix<T,Dim>::matmul2_dne_ee_edn(a, b, N, buffer, c);
  // validate
  Matmul<T>::matmul(
    CblasRowMajor, CblasTrans, CblasNoTrans,
    _N, _N2, _N2, (T)1.0, a, _N, b, _N2, (T)0.0, d, _N2);
  Matmul<T>::matmul(
    CblasRowMajor, CblasNoTrans, CblasNoTrans,
    _N, _N, _N2, (T)1.0, d, _N2, a, _N, (T)0.0, e, _N);
  if (layout) {
    LOG(INFO) << "matrix a layout";
    print_mat<T>(a, _N2, _N);
    LOG(INFO) << "matrix b layout";
    print_mat<T>(b, _N2, _N2);
    LOG(INFO) << "matrix c layout";
    print_mat<T>(c, _N, _N);
    LOG(INFO) << "matrix e layout";
    print_mat<T>(e, _N, _N);
  }
  bool flag = validate<T>(c, e, _N*_N, tol);
  // free
  free(a);
  free(b);
  free(buffer);
  free(c);
  free(d);
  free(e);
  if (flag) {
    LOG(INFO) << "test_matmul2_dne_ee_edn succeed, T: " << typeid(T).name()
      << ", Dim: " << Dim << ", N: " << N;
  } else {
    LOG(FATAL) << "test_matmul2_dne_ee_edn failed, T: " << typeid(T).name()
      << ", Dim: " << Dim << ", N: " << N;
  }
}

template <typename T, unsigned int Dim, unsigned int N>
void test_matmul2_dnf_ff_fdn(bool layout=false, double tol=1e-6) {
  // init a
  unsigned int _N = Dim * N;
  unsigned int _N2 = Dim * Dim;
  unsigned int nEntryA = _N2 * _N;
  unsigned int nBytesA = nEntryA * sizeof(T);
  auto a = (T*)malloc(nBytesA);
  init_rand<T>(a, nEntryA);
  // init b
  unsigned int nEntryB = _N2 * _N2;
  unsigned int nBytesB = nEntryB * sizeof(T);
  auto b = (T*)malloc(nBytesB);
  init_rand<T>(b, nEntryB);
  // init buffer
  unsigned int nEntryBuffer = _N2 * _N2;
  unsigned int nBytesBuffer = nEntryBuffer * sizeof(T);
  auto buffer = (T*)malloc(nBytesBuffer);
  // init c
  unsigned int nEntryC = _N * _N;
  unsigned int nBytesC = nEntryC * sizeof(T);
  auto c = (T*)malloc(nBytesC);
  // init d
  unsigned int nEntryD = _N * _N2;
  unsigned int nBytesD = nEntryD * sizeof(T);
  auto d = (T*)malloc(nBytesD);
  // init e
  unsigned int nEntryE = _N * _N;
  unsigned int nBytesE = nEntryE * sizeof(T);
  auto e = (T*)malloc(nBytesE);
  // execute
  FEMatrix<T,Dim>::matmul2_dnf_ff_fdn(a, b, N, buffer, c);
  // validate
  Matmul<T>::matmul(
    CblasRowMajor, CblasTrans, CblasNoTrans,
    _N, _N2, _N2, (T)1.0, a, _N, b, _N2, (T)0.0, d, _N2);
  Matmul<T>::matmul(
    CblasRowMajor, CblasNoTrans, CblasNoTrans,
    _N, _N, _N2, (T)1.0, d, _N2, a, _N, (T)0.0, e, _N);
  if (layout) {
    LOG(INFO) << "matrix a layout";
    print_mat<T>(a, _N2, _N);
    LOG(INFO) << "matrix b layout";
    print_mat<T>(b, _N2, _N2);
    LOG(INFO) << "matrix c layout";
    print_mat<T>(c, _N, _N);
    LOG(INFO) << "matrix e layout";
    print_mat<T>(e, _N, _N);
  }
  bool flag = validate<T>(c, e, _N*_N, tol);
  // free
  free(a);
  free(b);
  free(buffer);
  free(c);
  free(d);
  free(e);
  if (flag) {
    LOG(INFO) << "test_matmul2_dnf_ff_fdn succeed, T: " << typeid(T).name() 
      << ", Dim: " << Dim << ", N: " << N;
  } else {
    LOG(FATAL) << "test_matmul2_dnf_ff_fdn failed, T: " << typeid(T).name()
      << ", Dim: " << Dim << ", N: " << N;
  }
}

int main(int argc, char* argv[]) {
  // log init
  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr = 1;
  // double precison tests
  test_matmul_cblas<double>(false);
  test_matmul_nddd<double,3,8>();
  test_matmul_dnnd<double,3,8>();
  test_inv_dd<double,3>();
  test_mattile_diag_dd<double,3>();
  test_matmul2_dne_ee_edn<double,3,8>();
  test_matmul2_dnf_ff_fdn<double,3,8>();
  LOG(INFO) << "double precision test passed";
  // single precision tests
  test_matmul_cblas<float>(false);
  test_matmul_nddd<float,3,8>();
  test_matmul_dnnd<float,3,8>();
  test_inv_dd<float,3>();
  test_mattile_diag_dd<float,3>();
  test_matmul2_dne_ee_edn<float,3,8>();
  test_matmul2_dnf_ff_fdn<float,3,8>();
  LOG(INFO) << "single precision test passed";
  return 0;
}
