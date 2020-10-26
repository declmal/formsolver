#include <fem/fe_element.h>
#include <fem/fe_matrix.h>
#include <fem/fe_brick.h>
#include <fem/formulator/total_lagrangian.h>
#include <common/common.h>
#include <iostream> // for debug

namespace fem {
FORM_REGISTER_ELEMENT_TEMPLATE()
IPropType<T> Element<T,NI,N,IPropType,FormType>::iprop;

FORM_REGISTER_ELEMENT_TEMPLATE()
FormType<T> Element<T,NI,N,IPropType,FormType>::form;

FORM_REGISTER_ELEMENT_TEMPLATE()
void Element<T,NI,N,IPropType,FormType>::form_elem_stiff() {
  auto h = iprop.hbuf;
  auto weights = iprop.weights;
  init_zero<T>(Ke, 9*N*N);
  for (unsigned int i = 0; i < NI; ++i) {
    std::cout << weights[i] << std::endl;
    h += NI;
  }
}

FORM_REGISTER_ELEMENT_TEMPLATE()
void Element<T,NI,N,IPropType,FormType>::init_coordinate() {
  // TODO
}

FORM_REGISTER_ELEMENT_TEMPLATE()
void Element<T,NI,N,IPropType,FormType>::init_coordinate(
  const T* const data, const unsigned int size) {
  for (unsigned int i = 0; i < size; ++i) {
    X0[i] = data[i];
  }
}

// C3D8 TL
FORM_REGISTER_ELEMENT(double, 8, 8, C3D8IProp, TL3D)
FORM_REGISTER_ELEMENT(float, 8, 8, C3D8IProp, TL3D)
// // C3D8R TL
// FORM_REGISTER_ELEMENT(double, 1, 8, C3D8RIProp, TL3D)
// FORM_REGISTER_ELEMENT(float, 1, 8, C3D8RIProp, TL3D)
// C3D20 TL
FORM_REGISTER_ELEMENT(double, 27, 20, C3D20IProp, TL3D)
FORM_REGISTER_ELEMENT(float, 27, 20, C3D20IProp, TL3D)
// C3D20R TL
FORM_REGISTER_ELEMENT(double, 8, 20, C3D20RIProp, TL3D)
FORM_REGISTER_ELEMENT(float, 8, 20, C3D20RIProp, TL3D)

} // namespace fem
