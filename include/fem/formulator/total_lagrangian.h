#ifndef FEM_FORMULATOR_TOTAL_LAGRANGIAN_H_
#define FEM_FORMULATOR_TOTAL_LAGRANGIAN_H_

#include <common/matrix.h>
#include <common/common.h>

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
   *  of shape (Dim*Dim, Dim*N)
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

#define FORM_REGISTER_TL_OP(T) \
  template \
  void lin_trans_mat_tl_3d( \
    const T* const h0, const T* const u0t, \
    const unsigned int N, T* const B0tL); \
  template \
  void nonlin_trans_mat_tl_3d( \
    const T* const h0, const unsigned int N, T* const B0NL);

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
   * \param J0, auxiliary variable, jacobian matrix with respect to 
   *  inital configuration, of shape (Dim, Dim)
   * \param invJ0, auxiliary variable, inversion of J0, of shape (Dim, Dim)
   * \param h0 auxiliary variable, interpolation derivative
   *  with respect to initial global coordinate, of shape (N, Dim)
   * \param u0t auxiliary variable, displacement derivative with respect to
   *  initial global coordinate, of shape (Dim, Dim)
   * \param B0tL auxiliary variable, linear strain incremental stiffness
   *  matrix, of shape (3*Dim-3, Dim*N)
   * \param Ke, output variable, element stiffness matrix, 
   *  of shape (Dim*N, Dim*N)
   */
  static int form_elem_stiff(
    const T* const X0, const T* const hbuf, const T* const Ut, T* const J0, 
    T* const invJ0, T* const h0, T* const u0t, T* const B0tL, T* const Ke) {
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
      MatmulNDDDT<T,Dim>::matmul_ndddt(h, invJ0, N, h0);
      h += stride_h;
      MatmulDNND<T,Dim>::matmul_dnnd(Ut, h0, N, u0t);
    }
    printf("hihihihhih\n\n\n\n\n");
    return 0;
  }
};

#define FORM_REGISTER_ELEM_FORM(ElemType, NI, N, Dim) \
  template <typename T> \
  using ElemType##TLForm = TLForm<T,NI,N,Dim>;
} // namespace fem

#endif // FEM_FORMULATOR_TOTAL_LAGRANGIAN_H_
