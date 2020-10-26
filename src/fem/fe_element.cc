#include <fem/fe_element.h>
#include <iostream>

namespace fem {
template <typename T, unsigned int NI, unsigned int N> 
void ElementVol<T,NI,N>::form_elem_stiff(
  const T* const Ke, const ElementVolIProp<T,NI,N> iprop) {
  auto h = iprop.hbuf;
  auto weights = iprop.weights;
  for (unsigned int i = 0; i < NI; ++i) {
    std::cout << weights[i] << std::endl;
    h += NI;
  }
}
template
void ElementVol<double,8,8>::form_elem_stiff(
  const double* const Ke, const ElementVolIProp<double,8,8> iprop);
template
void ElementVol<float,8,8>::form_elem_stiff(
  const float* const Ke, const ElementVolIProp<float,8,8> iprop);
} // namespace fem
