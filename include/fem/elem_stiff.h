#ifndef FEM_ELEM_STIFF_H_
#define FEM_ELEM_STIFF_H_

namespace fem {
/*!
 * \brief linear transformation matrix (Total Lagrangian Version)
 *
 * \param h0 input variable, interpolation derivative 
 *  with respect to initial global coordinate, of shape (N, 3)
 * \param u0t input variable, temporal displacement derivative
 *  with respect to initial global coordinate, of shape (3, 3)
 * \param N input variable, the number of nodes in an element
 * \param B0t_L output variable, linear transformation matrix, of shape (6, 3N)
*/
template <typename T>
void lin_trans_mat_tl(
  const T* const h0, const T* const u0t, const unsigned int N, T* const B0t_L);
} // namespace fem

#endif // FEM_ELEM_STIFF_H_
