#include <fem/fe_element.h>
#include <iostream>

namespace fem {
template <typename T, unsigned int NI, unsigned int N> 
void ElementVol<T, NI, N>::form_elem_stiff(const T* const Ke) {
  auto h = iprop.hbuf;
  auto weights = iprop.weights;
  for (unsigned int i = 0; i < NI; ++i) {
    std::cout << weights[i] << std::endl;
    h += NI;
  }
}
template void ElementVol<double, 8, 8>::form_elem_stiff(
  const double* const Ke);
} // namespace fem
