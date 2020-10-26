#include <fem/fe_element.h>
#include <fem/fe_matrix.h>
#include <common/common.h>
#include <iostream> // debug

namespace fem {
template <typename T, unsigned int NI, unsigned int N> 
void Element<T,NI,N>::form_elem_stiff_cpu(
  const ElementIProp<T,NI,N> iprop) {
  auto h = iprop.hbuf;
  auto weights = iprop.weights;
  init_zero<T>(Ke, 9*N*N);
  for (unsigned int i = 0; i < NI; ++i) {
    h += NI;
  }
}

template <typename T, unsigned int NI, unsigned int N>
void Element<T,NI,N>::init_coordinate() {
  // TODO
}

template <typename T, unsigned int NI, unsigned int N>
void Element<T,NI,N>::init_coordinate(
  const T* const data, const unsigned int size) {
  for (unsigned int i = 0; i < size; ++i) {
    X0[i] = data[i];
  }
}

// C3D8
REGISTER_ELEMENT(double, 8, 8)
REGISTER_ELEMENT(float, 8, 8)
// C3D8R
REGISTER_ELEMENT(double, 1, 8)
REGISTER_ELEMENT(float, 1, 8)
// C3D20
REGISTER_ELEMENT(double, 27, 20)
REGISTER_ELEMENT(float, 27, 20)
// C3D20R
REGISTER_ELEMENT(double, 8, 20)
REGISTER_ELEMENT(float, 8, 20)

} // namespace fem
