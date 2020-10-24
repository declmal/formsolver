#ifndef FEM_FE_ELEMENT_H_
#define FEM_FE_ELEMENT_H_

namespace fem {
/*!
 * \brief 3-d Isoparametric 8-node Element Shape Function Derivatives
 *
 * \param r input variable, natural coordinates of point, of shape (3,)
 * \param h output variable, shape function derivatives, of shape (3, 8)
 *  shape function derivatives at point r = [r0,r1,r2]
 *  h[j][n] = dhn/drj, j = 0,1,2
 *  hn(r0,r1,r2) = Gn0(r0) * Gn1(r1) * Gn2(r2)
 *  Gnj(rj) = 0.5 * (1+Rnj*rj)
 *  Rnj = [[1,1,1],
 *         [-1,1,1],
 *         [-1,-1,1],
 *         [1,-1,1],
 *         [1,1,-1],
 *         [-1,1,-1],
 *         [-1,-1,-1],
 *         [1,-1,-1]]
 */
template <typename T>
void brickn8_shp_deri(const T* const r, T* const h);

/*!
 * \brief Gausss-Legendre Integration for Isoparametric 3d Brick Element
 *
 * \param r output variable, interpolation points, of shape (3, 8)
 * \param w output variable, interpolation weights, of shape (8,)
 */
template <typename T>
void brickn8_gauss_i2(T* const r, T* const w) {
  for (unsigned int i = 0; i < 2; ++i) {
    for (unsigned int j = 0; j < 2; ++j) {
      for (unsigned int k = 0; k < 2; ++k) {
      }
    }
  }
}
} // namespace fem

#endif // FEM_FE_ELEMENT_H_
