#ifndef FEM_TRANS_MAT_H_
#define FEM_TRANS_MAT_H_

namespace fem {
/*!
 * \brief Determinant of Jacobian Matrix
 *
 * \param J input variable, Jacobi matrix, of shape (3, 3)
 * \param detJ output variable, determinant Jacobi matrix, of shape (1,)
*/
template <typename T>
void det_jacobi(const T* const J, T* const detJ);

/*!
 * \brief Interpolation Derivative with Respect to Global Coordinate
 *
 * \param Jinv input variable, inverse of Jacobi matrix with respect to
 *  initial (TL) or temporal (UL) global coordinate, of shape (3, 3)
 * \param h input variable, interpolation derivative 
 *  with respect to local natural coordinate, of shape (N, 3)
 * \param N input variable, the number of nodes in an element
 * \param htau output variable, interpolation derivative with respect to
 *  initial (TL) or temporal (UL) global coordinate, of shape (N, 3)
*/
template <typename T>
void intr_deriv(
  const T* const Jinv, const T* const h, const unsigned int N, T* const htau);

/*!
 * \brief Temporal Displacement Derivative with Respect to 
 * Initial Global Coordinate (Total Lagrangian Only)
 *
 * \param h0 input variable, interpolation derivative 
 *  with respect to initial global coordinate, of shape (N, 3)
 * \param Ut input variable, temporal displacement, of shape (3, N)
 *  with respect to global coordinate, of shape (3, 3)
 * \param N input variable, the number of nodes in an element
 * \param u0t output variable, temporal displacement derivative
 *  with respect to initial global coordinate, of shape (3, 3)
*/
template <typename T>
void disp_deriv_tl(
  const T* const h0, const T* const Ut, const unsigned int N, T* const u0t);

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

template <typename T>
void nonlin_trans_mat_tl();
} // namespace fem

#endif // FEM_TRANS_MAT_H_
