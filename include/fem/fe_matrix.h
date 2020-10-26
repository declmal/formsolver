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
 * \brief Determinant of 3x3 Matrix
 *
 * \param a input variable, matrix of shape (3, 3)
 * \param det output variable, determinant of matrix, of shape (1,)
 */
template <typename T>
void det_33(const T* const a, T* const det);

/*!
 * \brief Inversion of 3x3 Matrix
 *
 * \param a input variable, matrix of shape (3, 3)
 * \param det intput variable, determinant of matrix, of shape (1,)
 * \param inv output variable, inversion matrix, of shape (3, 3),
 */
template <typename T>
void inv_33(const T* const a, const T* const det, T* const inv);

/*!
 * \brief Diagnal Tile of 3x3 Matrix to 9x9
 *
 * \param a input variable, matrix of shape (3, 3)
 * \param tile output variable, inversion matrix, of shape (9, 9),
 */
template <typename T>
void mattile_diag_33(const T* const a, T* const tile);

/*!
 * \brief Matrix Multiplication, cij = ali blk akj
 *  i = 0,1,...,3N-1; j = 0,1,...,3N-1; k = 0,1,...,5, l = 0,1,...,5
 *
 * \param a input variable, matrix of shape (6, 3N)
 * \param b input variable, matrix of shape (6, 6)
 * \param N input variable, dimension variable
 * \param buffer output variable, matrix of shape (6,)
 * \param c output variable, matrix of shape (3N, 3N)
 */
template <typename T>
void matmul2_3n6_66_63n(
  const T* const a, const T* const b, const unsigned int N,
  T* const buffer, T* const c);

/*!
 * \brief Matrix Multiplication, cij = ali blk akj
 *  i = 0,1,...,3N-1; j = 0,1,...,3N-1; k = 0,1,...,8, l = 0,1,...,8
 *
 * \param a input variable, matrix of shape (9, 3N)
 * \param b input variable, matrix of shape (9, 9)
 * \param N input variable, dimension variable
 * \param buffer output variable, matrix of shape (9,)
 * \param c output variable, matrix of shape (3N, 3N)
 */
template <typename T>
void matmul2_3n9_99_93n(
  const T* const a, const T* const b, const unsigned int N,
  T* const buffer, T* const c);
} // namespace fem

#endif // FEM_FE_MATRIX_H_
