#ifndef FEM_FORMULATOR_TOTAL_LAGRANGIAN_H_
#define FEM_FORMULATOR_TOTAL_LAGRANGIAN_H_

#include <common/matrix.h>
#include <common/common.h>

namespace fem {
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
   * \param J0, auxiliary variable, jacobian matrix with respect to 
   *  inital configuration, of shape (Dim, Dim)
   * \param invJ0, auxiliary variable, inversion of J0, of shape (Dim, Dim)
   * \param Ke, output variable, element stiffness matrix, 
   *  of shape (Dim*N, Dim*N)
   */
  static int form_elem_stiff(
    const T* const X0, const T* const hbuf, T* const J0, T* const invJ0,
    T* const Ke) {
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
      h += stride_h;
    }
    return 0;
  }
};

#define FORM_REGISTER_ELEM_FORM(ElemType, NI, N, Dim) \
  template <typename T> \
  using ElemType##TLForm = TLForm<T,NI,N,Dim>;
} // namespace fem

#endif // FEM_FORMULATOR_TOTAL_LAGRANGIAN_H_
