#include <fem/fe_brick.h>

namespace fem {
template <typename T>
void c3d8_shp_deri(const T* const r, T* const h) {
  auto a0 = (T)1 + r[0];
  auto m0 = (T)1 - r[0];
  auto a1 = (T)1 + r[1];
  auto m1 = (T)1 - r[1];
  auto a2 = (T)1 + r[2];
  auto m2 = (T)1 - r[2];
  auto v = (T)0.125 * a1 * a2;
  h[0] = v;
  h[3] = -v;
  v = (T)0.125 * m1 * a2;
  h[6] = -v;
  h[9] = v;
  v = (T)0.125 * a1 * m2;
  h[12] = v;
  h[15] = -v;
  v = (T)0.125 * m1 * m2;
  h[18] = -v;
  h[21] = v;
  v = (T)0.125 * a0 * a2;
  h[1] = v;
  h[10] = -v;
  v = (T)0.125 * m0 * a2;
  h[4] = v;
  h[7] = -v;
  v = (T)0.125 * a0 * m2;
  h[13] = v;
  h[22] = -v;
  v = (T)0.125 * m0 * m2;
  h[16] = v;
  h[19] = -v;
  v = (T)0.125 * a0 * a1;
  h[2] = v;
  h[14] = -v;
  v = (T)0.125 * m0 * a1;
  h[5] = v;
  h[17] = -v;
  v = (T)0.125 * m0 * m1;
  h[8] = v;
  h[20] = -v;
  v = (T)0.125 * a0 * m1;
  h[11] = v;
  h[23] = -v;
}
template void c3d8_shp_deri(const double* const r, double* const h);
template void c3d8_shp_deri(const float* const r, float* const h);
} // namespace fem
