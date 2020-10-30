#ifndef COMMON_MATRIX_H_
#define COMMON_MATRIX_H_

#include <cblas.h>

template <typename T>
void matmul_dnnd_3d(
  const T* const a, const T* const b, const unsigned int N, T* const c);

template <typename T, unsigned int Dim>
struct MatmulDNND {
  /*!
   * \brief Matrix Multiplication, cij = aik * bkj, 
   *  i = 0,1,...,Dim-1; j = 0,1,...,Dim-1; k = 0,1,...,N-1
   *
   * \param a input variable, matrix of shape (Dim, N) 
   * \param b input variable, matrix of shape (N, Dim) 
   * \param N input variable, dimension variable
   * \param c output variable, matrix of shape (Dim, Dim)
   */
  static inline void matmul_dnnd(
    const T* const a, const T* const b, const unsigned int N, T* const c);
};
template <typename T>
struct MatmulDNND<T,3> {
  static inline void matmul_dnnd(
    const T* const a, const T* const b, const unsigned int N, T* const c) {
    matmul_dnnd_3d(a, b, N, c);
  }
};

template <typename T>
void matmul_nddd_3d(
  const T* const a, const T* const b, const unsigned int N, T* const c);

template <typename T, unsigned int Dim>
struct MatmulNDDD {
  /*!
   * \brief Matrix Multiplication, ckj = aki * bij, 
   *  i = 0,1,...,Dim-1; j = 0,1,...,Dim-1; k = 0,1,...,N-1
   *
   * \param a input variable, matrix of shape (N, Dim) 
   * \param b input variable, matrix of shape (Dim, Dim) 
   * \param N input variable, dimension variable
   * \param c output variable, matrix of shape (N, Dim)
   */
  static inline void matmul_nddd(
    const T* const a, const T* const b, const unsigned int N, T* const c);
};
template <typename T>
struct MatmulNDDD<T,3> {
  static inline void matmul_nddd(
    const T* const a, const T* const b, const unsigned int N, T* const c) {
    matmul_nddd_3d(a, b, N, c);
  }
};

template <typename T>
void matmul_ndddt_3d(
  const T* const a, const T* const b, const unsigned int N, T* const c);

template <typename T, unsigned int Dim>
struct MatmulNDDDT {
  /*!
   * \brief Matrix Multiplication, ckj = aki * transpose(b)ij, 
   *  i = 0,1,...,Dim-1; j = 0,1,...,Dim-1; k = 0,1,...,N-1
   *
   * \param a input variable, matrix of shape (N, Dim) 
   * \param b input variable, matrix of shape (Dim, Dim) 
   * \param N input variable, dimension variable
   * \param c output variable, matrix of shape (N, Dim)
   */
  static inline void matmul_ndddt(
    const T* const a, const T* const b, const unsigned int N, T* const c);
};
template <typename T>
struct MatmulNDDDT<T,3> {
  static inline void matmul_ndddt(
    const T* const a, const T* const b, const unsigned int N, T* const c) {
    matmul_ndddt_3d(a, b, N, c);
  }
};

template <typename T>
T det_dd_3d(const T* const a);

template <typename T, unsigned int Dim>
struct DetDD {
  /*!
   * \brief Determinant of Dim x Dim Matrix
   *
   * \param a input variable, matrix of shape (Dim, Dim)
   * \return det, determinant of matrix
   */
  static inline T det_dd(const T* const a);
};
template <typename T>
struct DetDD<T,3> {
  static constexpr double tol = 1e-6;
  static inline T det_dd(const T* const a) {
    return det_dd_3d(a);
  }
};

template <typename T>
void inv_dd_3d(const T* const a, const T det, T* const inv);

template <typename T, unsigned int Dim>
struct InvDD {
  /*!
   * \brief Inversion of Dim x Dim Matrix
   *
   * \param a input variable, matrix of shape (Dim, Dim)
   * \param det intput variable, determinant of matrix
   * \param inv output variable, inversion matrix, of shape (Dim, Dim),
   */
  static inline void inv_dd(const T* const a, const T det, T* const inv);
};
template <typename T>
struct InvDD<T,3> {
  static inline void inv_dd(const T* const a, const T det, T* const inv) {
    inv_dd_3d(a, det, inv);
  }
};

template <typename T>
void mattile_diag_dd_3d(const T* const a, T* const tile);

template <typename T, unsigned int Dim>
struct MattileDiagDD {
  /*!
   * \brief Diagnal Tile of Dim x Dim Matrix to Dim^2 x Dim^2
   *
   * \param a input variable, matrix of shape (Dim, Dim)
   * \param tile output variable, inversion matrix, of shape (Dim^2, Dim^2),
   */
  static inline void mattile_diag_dd(const T* const a, T* const tile);
};
template <typename T>
struct MattileDiagDD<T,3> {
  static inline void mattile_diag_dd(const T* const a, T* const tile) {
    mattile_diag_dd_3d(a, tile);
  }
};

template <typename T>
void matmul2_dne_ee_edn_3d(
  const T* const a, const T* const b, const unsigned int N,
  T* const buffer, T* const c);

template <typename T, unsigned int Dim>
struct Matmul2DNEEEEDN {
  /*!
   * \brief Matrix Multiplication, cij = ali blk akj
   *  i = 0,1,...,3N-1; j = 0,1,...,3N-1; 
   *  k = 0,1,...,3*Dim-4; l = 0,1,...,3*Dim-4
   *
   * \param a input variable, matrix of shape (3*Dim-3, 3N)
   * \param b input variable, matrix of shape (3*Dim-3, 3*Dim-3)
   * \param N input variable, dimension variable
   * \param buffer output variable, matrix of shape (3*Dim-3,)
   * \param c output variable, matrix of shape (3N, 3N)
   */
  static inline void matmul2_dne_ee_edn(
    const T* const a, const T* const b, const unsigned int N,
    T* const buffer, T* const c);
};
template <typename T>
struct Matmul2DNEEEEDN<T,3> {
  static inline void matmul2_dne_ee_edn(
    const T* const a, const T* const b, const unsigned int N,
    T* const buffer, T* const c) {
    matmul2_dne_ee_edn_3d(a, b, N, buffer, c);
  }
};

template <typename T>
void matmul2_dnf_ff_fdn_3d(
  const T* const a, const T* const b, const unsigned int N,
  T* const buffer, T* const c);

template <typename T, unsigned int Dim>
struct Matmul2DNFFFFDN {
  /*!
   * \brief Matrix Multiplication, cij = ali blk akj
   *  i = 0,1,...,3N-1; j = 0,1,...,3N-1; 
   *  k = 0,1,...,Dim^2-1; l = 0,1,...,Dim^2-1
   *
   * \param a input variable, matrix of shape (Dim^2, 3N)
   * \param b input variable, matrix of shape (Dim^2, Dim^2)
   * \param N input variable, dimension variable
   * \param buffer output variable, matrix of shape (Dim^2,)
   * \param c output variable, matrix of shape (3N, 3N)
   */
  static inline void matmul2_dnf_ff_fdn(
    const T* const a, const T* const b, const unsigned int N,
    T* const buffer, T* const c);
};
template <typename T>
struct Matmul2DNFFFFDN<T,3> {
  static inline void matmul2_dnf_ff_fdn(
    const T* const a, const T* const b, const unsigned int N,
    T* const buffer, T* const c) {
    matmul2_dnf_ff_fdn_3d(a, b, N, buffer, c);
  }
};

template <typename T>
void matadd(
  const T* const a, const T* const b, const unsigned int size, T* const c);

template <typename T>
void matinc_mul(const T* const a, const unsigned int size, T* const b, T k);

template <typename T>
void matinc(const T* const a, const unsigned int size, T* const b);

template <typename T>
void matdec(const T* const a, const unsigned int size, T* const b);

#define FORM_REGISTER_MATRIX_OP(T) \
  template \
  void matmul_dnnd_3d( \
    const T* const a, const T* const b, const unsigned int N, T* const c); \
  template \
  void matmul_nddd_3d( \
    const T* const a, const T* const b, const unsigned int N, T* const c); \
  template \
  void matmul_ndddt_3d( \
    const T* const a, const T* const b, const unsigned int N, T* const c); \
  template \
  T det_dd_3d(const T* const a); \
  template \
  void inv_dd_3d(const T* const a, const T det, T* const inv); \
  template \
  void mattile_diag_dd_3d(const T* const a, T* const tile); \
  template \
  void matmul2_dne_ee_edn_3d( \
    const T* const a, const T* const b, const unsigned int N, \
    T* const buffer, T* const c); \
  template \
  void matmul2_dnf_ff_fdn_3d( \
    const T* const a, const T* const b, const unsigned int N, \
    T* const buffer, T* const c); \
  template \
  void matadd( \
    const T* const a, const T* const b, const unsigned int size, T* const c); \
  template \
  void matinc_mul(const T* const a, const unsigned int size, T* const b, T k); \
  template \
  void matinc(const T* const a, const unsigned int size, T* const b); \
  template \
  void matdec(const T* const a, const unsigned int size, T* const b);

template <typename T>
struct Matmul {
  static inline void matmul(
    CBLAS_ORDER order, CBLAS_TRANSPOSE trans_a, CBLAS_TRANSPOSE trans_b,
    const unsigned int M, const unsigned int N, const unsigned int K,
    const double alpha, const T* a, const unsigned int lda, const T* b,
    const unsigned int ldb, const double beta, T* c, 
    const unsigned int ldc);
};
template <>
struct Matmul<double> {
  static inline void matmul(
    CBLAS_ORDER order, CBLAS_TRANSPOSE trans_a, CBLAS_TRANSPOSE trans_b,
    const unsigned int M, const unsigned int N, const unsigned int K,
    const double alpha, const double* a, const unsigned int lda, const double* b,
    const unsigned int ldb, const double beta, double* c, 
    const unsigned int ldc) {
    cblas_dgemm(
      order, trans_a, trans_b, M, N, K, alpha, a, lda, b, ldb, beta, c, ldc);
  }
};
template <>
struct Matmul<float> {
  static inline void matmul(
    CBLAS_ORDER order, CBLAS_TRANSPOSE trans_a, CBLAS_TRANSPOSE trans_b,
    const unsigned int M, const unsigned int N, const unsigned int K,
    const float alpha, const float* a, const unsigned int lda, const float* b,
    const unsigned int ldb, const float beta, float* c, 
    const unsigned int ldc) {
    cblas_sgemm(
      order, trans_a, trans_b, M, N, K, alpha, a, lda, b, ldb, beta, c, ldc);
  }
};

#endif // COMMON_MATRIX_H_
