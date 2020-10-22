#ifndef FEM_FE_MATRIX_H_
#define FEM_FE_MATRIX_H_

namespace fem {
/*!
 * \brief Matrix Multiplication, cij = aik * bkj, 
 *  i = 0,1,2, j = 0,1,2, k = 0,1,...,N-1
 *
 * \param a input variable, matrix of shape (3, N) 
 * \param b input variable, matrix of shape (N, 3) 
 * \param N input variable, dimension variable
 * \param c output variable, matrix of shape (3, 3)
*/
template <typename T>
void matmul_3nn3(const T* const a, const T* const b, const unsigned int N, T* const c);

/*!
 * \brief Matrix Multiplication, ckj = aki * bij, 
 *  i = 0,1,2, j = 0,1,2, k = 0,1,...,N-1
 *
 * \param a input variable, matrix of shape (N, 3) 
 * \param b input variable, matrix of shape (3, 3) 
 * \param N input variable, dimension variable
 * \param c output variable, matrix of shape (N, 3)
*/
template <typename T>
void matmul_n333(const T* const a, const T* const b, const unsigned int N, T* const c);

/*!
 * \brief Determinant of Jacobian Matrix
 *
 * \param a input variable, matrix of shape (3, 3)
 * \param det output variable, determinant of matrix, of shape (1,)
*/
template <typename T>
void det_33(const T* const a, T* const det);

/*!
 * \brief Inversion of Jacobian Matrix
 *
 * \param a input variable, matrix of shape (3, 3)
 * \param det intput variable, determinant of matrix, of shape (1,)
 * \param inv output variable, inversion matrix, of shape (3, 3),
*/
template <typename T>
void inv_33(const T* const a, const T* const det, T* const inv);
} // namespace fem

#endif // FEM_FE_MATRIX_H_
