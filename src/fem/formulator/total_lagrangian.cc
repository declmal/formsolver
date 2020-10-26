#include <fem/formulator/total_lagrangian.h>
#include <fem/element/brick.h>
#include <iostream> // for debug

namespace fem {
FORM_REGISTER_FORMULATOR_TEMPLATE()
void TotalLagrangian<T,Dim,EType>::load_elem_data() {
  X0 = Formulator<T,Dim,EType>::elem->get_X0();
  J = Formulator<T,Dim,EType>::elem->get_J();
  Ke = Formulator<T,Dim,EType>::elem->get_Ke();
}


FORM_REGISTER_FORMULATOR_TEMPLATE()
void TotalLagrangian<T,Dim,EType>::form_elem_stiff() {
}

FORM_REGISTER_TL(float, 3, C3D8)

} // namespace fem
