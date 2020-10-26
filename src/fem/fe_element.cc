#include <fem/fe_element.h>
#include <fem/fe_matrix.h>
#include <fem/fe_brick.h>
#include <common/common.h>
#include <iostream> // debug

namespace fem {
template <
  typename T, unsigned int NI, unsigned int N,
  template <typename> class IPropType>
IPropType<T> Element<T,NI,N,IPropType>::iprop;

template <
  typename T, unsigned int NI, unsigned int N,
  template <typename> class IPropType>
void Element<T,NI,N,IPropType>::form_elem_stiff_cpu() {
  auto h = iprop.hbuf;
  auto weights = iprop.weights;
  init_zero<T>(Ke, 9*N*N);
  for (unsigned int i = 0; i < NI; ++i) {
    std::cout << weights[i] << std::endl;
    h += NI;
  }
}

template <
  typename T, unsigned int NI, unsigned int N,
  template <typename> class IPropType>
void Element<T,NI,N,IPropType>::init_coordinate() {
  // TODO
}

template <
  typename T, unsigned int NI, unsigned int N,
  template <typename> class IPropType>
void Element<T,NI,N,IPropType>::init_coordinate(
  const T* const data, const unsigned int size) {
  for (unsigned int i = 0; i < size; ++i) {
    X0[i] = data[i];
  }
}

// C3D8
FORM_REGISTER_ELEMENT(double, 8, 8, C3D8IProp)
FORM_REGISTER_ELEMENT(float, 8, 8, C3D8IProp)
// // C3D8R
// FORM_REGISTER_ELEMENT(double, 1, 8, C3D8RIProp)
// FORM_REGISTER_ELEMENT(float, 1, 8, C3D8RIProp)
// C3D20
FORM_REGISTER_ELEMENT(double, 27, 20, C3D20IProp)
FORM_REGISTER_ELEMENT(float, 27, 20, C3D20IProp)
// C3D20R
FORM_REGISTER_ELEMENT(double, 8, 20, C3D20RIProp)
FORM_REGISTER_ELEMENT(float, 8, 20, C3D20RIProp)

} // namespace fem
