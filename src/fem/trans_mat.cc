namespace fem {
/*!
 * \brief linear transformation matrix (total lagrangian version)
 *
 * \param h0 input variable, interpolation derivative 
 *  with respect to global coordinate, of shape (N, 3)
 * \param u0t input variable, displacement derivative
 *  with respect to global coordinate, of shape (3, 3)
 * \param N input variable, number of nodes in an element
 * \param B0t_L output variable, linear transformation matrix
*/
template <typename T>
void lin_trans_mat_tl(
  const T* h0, const T* u0t, const unsigned int N, T* B0t_L) {
  int _3N = 3 * N;
  T* B0 = B0t_L;
  T* B1 = B0 + _3N;
  T* B2 = B1 + _3N;
  T* B3 = B2 + _3N;
  T* B4 = B3 + _3N;
  T* B5 = B4 + _3N;
  T* H = h0;
  for (int i = 0; i < N; ++i) {
    B0[0] = H[0] * ((T)1 + u0t[0]);
    B0[1] = H[0] * u0t[3];
    B0[2] = H[0] * u0t[6];
    B1[0] = H[1] * u0t[1];
    B1[1] = H[1] * ((T)1 + u0t[4]);
    B1[2] = H[1] * u0t[7];
    B2[0] = H[2] * u0t[2];
    B2[1] = H[2] * u0t[5];
    B2[2] = H[2] * ((T)1 + u0t[8]);
    B0 += 3;
    B1 += 3;
    B2 += 3;
    B3 += 3;
    B4 += 3;
    B5 += 3;
    H += 5;
  }
}

template <typename T>
void nonlin_trans_mat_tl() {
}
} // namespace fem
