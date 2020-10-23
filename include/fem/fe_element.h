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
void brick3d_shp_deri(const T* const r, T* const h);
} // namespace fem

#endif // FEM_FE_ELEMENT_H_
