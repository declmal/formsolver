#ifndef FEM_ELEM_STIFF_H_
#define FEM_ELEM_STIFF_H_

namespace fem {
/*!
 * \brief 3d linear transformation matrix (Total Lagrangian Version)
 *
 * \param h0 input variable, interpolation derivative 
 *  with respect to initial global coordinate, of shape (N, 3)
 * \param u0t input variable, temporal displacement derivative
 *  with respect to initial global coordinate, of shape (3, 3)
 * \param N input variable, the number of nodes in an element
 * \param B0t_L output variable, linear transformation matrix, of shape (6, 3N)
 */
template <typename T>
void lin_trans_mat_tl_3d(
  const T* const h0, const T* const u0t, const unsigned int N, T* const B0t_L);

/*!
 * \brief 3d nonlinear transformation matrix (Total Lagrangian Version)
 *
 * \param h0 input variable, interpolation derivative 
 *  with respect to initial global coordinate, of shape (N, 3)
 * \param N input variable, the number of nodes in an element
 * \param B0_NL output variable, nonlinear transformation matrix, of shape (9, 3N)
 */
template <typename T>
void nonlin_trans_mat_tl_3d(
  const T* const h0, const unsigned int N, T* const B0_NL);
} // namespace fem

#endif // FEM_ELEM_STIFF_H_
