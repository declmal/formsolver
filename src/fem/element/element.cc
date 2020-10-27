#include <fem/element/element.h>
#include <fem/iprop/brick_iprop.h>
#include <common/common.h>

namespace fem {
FORM_REGISTER_ELEMENT_TEMPLATE()
IPropType<T> Element<T,N,Dim,IPropType>::iprop;

FORM_REGISTER_ELEMENT_TEMPLATE()
void Element<T,N,Dim,IPropType>::init_coordinate(
  const T* const data, const unsigned int size) {
  for (unsigned int i = 0; i < size; ++i) {
    X0[i] = data[i];
  }
}

FORM_REGISTER_ELEMENT_TEMPLATE()
T* Element<T,N,Dim,IPropType>::get_X0() {
  return X0;
}

FORM_REGISTER_ELEMENT_TEMPLATE()
T* Element<T,N,Dim,IPropType>::get_hbuf() {
  return iprop.hbuf;
}

FORM_REGISTER_ELEMENT_TEMPLATE()
T* Element<T,N,Dim,IPropType>::get_weights() {
  return iprop.weights;
}

FORM_REGISTER_ELEMENT_TEMPLATE()
unsigned int Element<T,N,Dim,IPropType>::get_num_ipoints() {
  return iprop.num_ipoints;
}

FORM_REGISTER_ELEMENT_TEMPLATE()
unsigned int Element<T,N,Dim,IPropType>::get_num_nodes() {
  return N;
}

FORM_REGISTER_ELEMENT_TEMPLATE()
void Element<T,N,Dim,IPropType>::init_coordinate() {
  // TODO
}

// C3D8
FORM_REGISTER_ELEMENT(double, 8, 3, C3D8IProp)
FORM_REGISTER_ELEMENT(float, 8, 3, C3D8IProp)
// // C3D8R
// FORM_REGISTER_ELEMENT(double, 8, 3, C3D8RIProp)
// FORM_REGISTER_ELEMENT(float, 8, 3, C3D8RIProp)
// C3D20
FORM_REGISTER_ELEMENT(double, 20, 3, C3D20IProp)
FORM_REGISTER_ELEMENT(float, 20, 3, C3D20IProp)
// C3D20R
FORM_REGISTER_ELEMENT(double, 20, 3, C3D20RIProp)
FORM_REGISTER_ELEMENT(float, 20, 3, C3D20RIProp)

} // namespace fem
