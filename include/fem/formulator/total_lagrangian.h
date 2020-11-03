#ifndef FEM_FORMULATOR_TOTAL_LAGRANGIAN_H_
#define FEM_FORMULATOR_TOTAL_LAGRANGIAN_H_

#include <common/matrix.h>
#include <common/common.h>
#include <iostream> // for debug purpose

namespace fem {
template <typename T>
void lin_trans_mat_tl_3d(
  const T* const h0, const T* const u0t, 
  const unsigned int N, T* const B0tL);

template <typename T, unsigned int Dim>
struct LinTransMatTL {
  /*!
   * \brief Linear Transformation Matrix 
   *  for Total Lagrangian Formulation
   *
   * \param h0 input variable, interpolation derivative 
   *  with respect to initial global coordinate, of shape (N, Dim)
   * \param u0t input variable, temporal displacement derivative
   *  with respect to initial global coordinate, of shape (Dim, Dim)
   * \param N input variable, the number of nodes in an element
   * \param B0tL output variable, linear transformation matrix, 
   *  of shape (3*Dim-3, Dim*N)
   */
  static inline void lin_trans_mat_tl(
    const T* const h0, const T* const u0t, 
    const unsigned int N, T* const B0tL);
};
template <typename T>
struct LinTransMatTL<T,3> {
  static inline void lin_trans_mat_tl(
    const T* const h0, const T* const u0t, 
    const unsigned int N, T* const B0tL) {
    lin_trans_mat_tl_3d(h0, u0t, N, B0tL);
  }
};

template <typename T>
void nonlin_trans_mat_tl_3d(
  const T* const h0, const unsigned int N, T* const B0NL);

template <typename T, unsigned int Dim>
struct NonlinTransMatTL {
  /*!
   * \brief Nonlinear Transformation Matrix 
   *  for Total Lagrangian Formulation
   *
   * \param h0 input variable, interpolation derivative 
   *  with respect to initial global coordinate, of shape (N, Dim)
   * \param N input variable, the number of nodes in an element
   * \param B0NL output variable, nonlinear transformation matrix, 
   *  of shape (Dim^2, Dim*N)
   */
  static inline void nonlin_trans_mat_tl(
    const T* const h0, const unsigned int N, T* const B0NL);
};
template <typename T>
struct NonlinTransMatTL<T,3> {
  static inline void nonlin_trans_mat_tl(
    const T* const h0, const unsigned int N, T* const B0NL) {
    nonlin_trans_mat_tl_3d(h0, N, B0NL);
  }
};

template <typename T>
void lin_trans_mat_3d(
  const T* const h0, const unsigned int N, T* const B);

template <typename T, unsigned int Dim>
struct LinTransMat {
  static inline void lin_trans_mat(
    const T* const h0, const unsigned int N, T* const B);
};
template <typename T>
struct LinTransMat<T,3> {
  static inline void lin_trans_mat(
    const T* const h0, const unsigned int N, T* const B) {
    lin_trans_mat_3d<T>(h0, N, B);
  }
};

template <typename T>
void trans_dil_mat_3d(
  const T* const h0, const unsigned int N, T* const B);

template <typename T, unsigned int Dim>
struct TransDivMat {
  static inline void trans_dil_mat(
    const T* const h0, const unsigned int N, T* const B);
};
template <typename T>
struct TransDivMat<T,3> {
  static inline void trans_dil_mat(
    const T* const h0, const unsigned int N, T* const B) {
    trans_dil_mat_3d<T>(h0, N, B);
  }
};

#define FORM_REGISTER_TL_OP(T) \
  template \
  void lin_trans_mat_tl_3d( \
    const T* const h0, const T* const u0t, \
    const unsigned int N, T* const B0tL); \
  template \
  void nonlin_trans_mat_tl_3d( \
    const T* const h0, const unsigned int N, T* const B0NL); \
  template \
  void lin_trans_mat_3d( \
    const T* const h0, const unsigned int N, T* const B); \
  template \
  void trans_dil_mat_3d( \
    const T* const h0, const unsigned int N, T* const B);

/*!
 * \brief Total Lagrangian Formulation (Single Element Version)
 */
template <typename T, unsigned int NI, unsigned int N, unsigned int Dim>
struct TLForm {
  /*!
   * \brief 
   *
   * \param X0 input variable, initial element nodal coordiate,
   *  of shape (Dim, N)
   * \param hbuf input variable, buffer of interpolation derivative
   *  with respect to natural coordinate, of shape (NI, N, Dim)
   * \param Ut input variable, temporal element nodal displacement,
   *  of shape (Dim, N)
   * \param C0 input variable, ,
   *  of shape (3*Dim-3, 3*Dim-3)
   * \param S0t input variable, ,
   *  of shape (Dim, Dim)
   * \param J0, auxiliary variable, jacobian matrix with respect to 
   *  inital configuration, of shape (Dim, Dim)
   * \param invJ0, auxiliary variable, inversion of J0, of shape (Dim, Dim)
   * \param h0 auxiliary variable, interpolation derivative
   *  with respect to initial global coordinate, of shape (N, Dim)
   * \param u0t auxiliary variable, displacement derivative with respect to
   *  initial global coordinate, of shape (Dim, Dim)
   * \param B0tL auxiliary variable, linear strain incremental stiffness
   *  matrix, of shape (3*Dim-3, Dim*N)
   * \param buf auxiliary variable, of shape (Dim^2,)
   * \param tmpK auxiliary variable, temporary stiffness matrix component,
   *  of shape (Dim*N, Dim*N)
   * \param B0NL auxiliary variable, no- linear strain incremental stiffness
   *  matrix, of shape (Dim^2, Dim*N)
   * \param tile auxiliary variable, of shape (Dim^2, Dim^2)
   * \param Ke, output variable, element stiffness matrix, 
   *  of shape (Dim*N, Dim*N)
   */
  static int form_elem_stiff(
    const T* const X0, const T* const hbuf, const T* const weights,
    const T* const Ut, 
    const T* const C0, const T* const S0t, T* const J0, T* const invJ0, 
    T* const h0, T* const u0t, T* const B0tL, T* const buf, T* const tmpK,
    T* const B0NL, T* const tile, T* const Ke) {
    auto rowKe = Dim * N;
    auto nEntryKe = rowKe * rowKe;
    init_zero<T>(Ke, nEntryKe);
    auto h = hbuf;
    auto stride_h = Dim * N;
    for (unsigned int i = 0; i < NI; ++i) {
      MatmulDNND<T,Dim>::matmul_dnnd(X0, h, N, J0);
      auto detJ0 = DetDD<T,Dim>::det_dd(J0);
      if ((double)abs(detJ0) < DetDD<T,Dim>::tol) {
        printf("singular jacobian matrix\n");
        return -1;
      }
      InvDD<T,Dim>::inv_dd(J0, detJ0, invJ0);
      MatmulNDDD<T,Dim>::matmul_nddd(h, invJ0, N, h0);
      h += stride_h;
      MatmulDNND<T,Dim>::matmul_dnnd(Ut, h0, N, u0t);
      LinTransMatTL<T,Dim>::lin_trans_mat_tl(h0, u0t, N, B0tL);
      Matmul2DNEEEEDN<T,Dim>::matmul2_dne_ee_edn(B0tL, C0, N, buf, tmpK);
      matinc_mul<T>(tmpK, nEntryKe, Ke, abs(detJ0));
      NonlinTransMatTL<T,Dim>::nonlin_trans_mat_tl(h0, N, B0NL);
      MattileDiagDD<T,Dim>::mattile_diag_dd(S0t, tile);
      Matmul2DNFFFFDN<T,Dim>::matmul2_dnf_ff_fdn(B0NL, tile, N, buf, tmpK);
      matinc_mul<T>(tmpK, nEntryKe, Ke, weights[i]*abs(detJ0));
    }
    printf("hihihihhih\n\n\n\n\n");
    return 0;
  }

  /* \brief Linear Element Stiffness Matrix Formulation (Only for Debug Purpose)
   *
   */
  static int form_linear_elem_stiff(
    const T* const X0, const T* const hbuf, const T* const weights,
    const T* const C, T* const J0, T* const invJ0, 
    T* const h0, T* const B, T* const buf, T* const tmpK,
    T* const Ke) {
    auto rowKe = Dim * N;
    auto nEntryKe = rowKe * rowKe;
    init_zero<T>(Ke, nEntryKe);
    auto h = hbuf;
    auto stride_h = Dim * N;
    for (unsigned int i = 0; i < NI; ++i) {
      MatmulDNND<T,Dim>::matmul_dnnd(X0, h, N, J0);
      auto detJ0 = DetDD<T,Dim>::det_dd(J0);
      auto absDetJ0 = abs(detJ0);
      if ((double)absDetJ0 < DetDD<T,Dim>::tol) {
        printf("singular jacobian matrix\n");
        return -1;
      }
      InvDD<T,Dim>::inv_dd(J0, detJ0, invJ0);
      MatmulNDDD<T,Dim>::matmul_nddd(h, invJ0, N, h0);
      h += stride_h;
      LinTransMat<T,Dim>::lin_trans_mat(h0, N, B);
      Matmul2DNEEEEDN<T,Dim>::matmul2_dne_ee_edn(B, C, N, buf, tmpK);
      matinc_mul<T>(tmpK, nEntryKe, Ke, weights[i]*absDetJ0);
    }
    return 0;
  }

  static int form_linear_elem_stiff_Bbar(
    const T* const X0, const T* const hbuf, const T* const weights,
    const T* const C, T* const J0, T* const invJ0, 
    T* const BdilBar, T* const tmpB,
    T* const h0, T* const B, T* const buf, T* const tmpK,
    T* const Ke) {
    // TODO buf h0
    auto h = hbuf;
    auto stride_h = Dim * N;
    T V0 = (T)0;
    for (unsigned int i = 0; i < NI; ++i) {
      MatmulDNND<T,Dim>::matmul_dnnd(X0, h, N, J0);
      auto detJ0 = DetDD<T,Dim>::det_dd(J0);
      auto absDetJ0 = abs(detJ0);
      if ((double)absDetJ0 < DetDD<T,Dim>::tol) {
        printf("singular jacobian matrix 1\n");
        return -1;
      }
      V0 += weights[i] * absDetJ0;
      h += stride_h;
    }
    h = hbuf;
    auto nRowBdilBar = 3*Dim - 3;
    auto nColBdilBar = Dim * N;
    auto nEntryBdilBar = nRowBdilBar * nColBdilBar;
    init_zero<T>(BdilBar, nEntryBdilBar);
    for (unsigned int i = 0; i < NI; ++i) {
      MatmulDNND<T,Dim>::matmul_dnnd(X0, h, N, J0);
      auto detJ0 = DetDD<T,Dim>::det_dd(J0);
      auto absDetJ0 = abs(detJ0);
      if ((double)absDetJ0 < DetDD<T,Dim>::tol) {
        printf("singular jacobian matrix 2\n");
        return -1;
      }
      InvDD<T,Dim>::inv_dd(J0, detJ0, invJ0);
      MatmulNDDD<T,Dim>::matmul_nddd(h, invJ0, N, h0);
      h += stride_h;
      TransDivMat<T,Dim>::trans_dil_mat(h0, N, tmpB);
      matinc_mul<T>(tmpB, nEntryBdilBar, BdilBar, absDetJ0*weights[i]/V0);
    }
    h = hbuf;
    auto rowKe = Dim * N;
    auto nEntryKe = rowKe * rowKe;
    init_zero<T>(Ke, nEntryKe);
    for (unsigned int i = 0; i < NI; ++i) {
      MatmulDNND<T,Dim>::matmul_dnnd(X0, h, N, J0);
      auto detJ0 = DetDD<T,Dim>::det_dd(J0);
      auto absDetJ0 = abs(detJ0);
      if ((double)absDetJ0 < DetDD<T,Dim>::tol) {
        printf("singular jacobian matrix\n");
        return -1;
      }
      InvDD<T,Dim>::inv_dd(J0, detJ0, invJ0);
      MatmulNDDD<T,Dim>::matmul_nddd(h, invJ0, N, h0);
      h += stride_h;
      LinTransMat<T,Dim>::lin_trans_mat(h0, N, B);
      TransDivMat<T,Dim>::trans_dil_mat(h0, N, tmpB);
      matdec<T>(tmpB, nEntryBdilBar, B);
      matinc<T>(BdilBar, nEntryBdilBar, B);
      Matmul2DNEEEEDN<T,Dim>::matmul2_dne_ee_edn(B, C, N, buf, tmpK);
      matinc_mul<T>(tmpK, nEntryKe, Ke, absDetJ0*weights[i]);
    }
    return 0;
  }
};

#define FORM_REGISTER_ELEM_FORM(ElemType, NI, N, Dim) \
  template <typename T> \
  using ElemType##TLForm = TLForm<T,NI,N,Dim>;
} // namespace fem

#endif // FEM_FORMULATOR_TOTAL_LAGRANGIAN_H_
