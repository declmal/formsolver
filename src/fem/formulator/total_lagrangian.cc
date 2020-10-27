#include <fem/formulator/total_lagrangian.h>
#include <fem/element/brick.h>
#include <common/matrix.h>
#include <common/common.h>
#include <iostream> // for debug

namespace fem {
FORM_REGISTER_FORMULATOR_TEMPLATE()
TotalLagrangian<T,Dim,EType>::TotalLagrangian(
  EType<T>* elem_) : Formulator<T,Dim,EType>(elem_) {
  X0 = Formulator<T,Dim,EType>::elem->get_X0();
  hbuf = Formulator<T,Dim,EType>::elem->get_hbuf();
  weights = Formulator<T,Dim,EType>::elem->get_weights();
  num_ipoints = Formulator<T,Dim,EType>::elem->get_num_ipoints();
  num_nodes = Formulator<T,Dim,EType>::elem->get_num_nodes();
}

FORM_REGISTER_FORMULATOR_TEMPLATE()
int TotalLagrangian<T,Dim,EType>::form_elem_stiff() {
  auto h = hbuf;
  // init Ke
  auto rowKe = Dim * num_nodes;
  auto nEntryKe = rowKe * rowKe;
  auto Ke = (T*)malloc(nEntryKe*sizeof(T));
  init_zero<T>(Ke, nEntryKe);
  // init J0
  auto nEntryJ0 = Dim * Dim;
  auto J0 = (T*)malloc(nEntryJ0*sizeof(T));
  for (unsigned int i = 0; i < num_ipoints; ++i) {
    FEMatrix<T,Dim>::matmul_dnnd(X0, h, num_nodes, J0);
    std::cout << "hihi" << std::endl;
    auto det = FEMatrix<T,Dim>::det_dd(J0);
    if (abs(det) < 1e-6) {
      free(Ke);
      free(J0);
    }
    h += num_ipoints;
  }
  // free
  free(Ke);
  free(J0);
  return 0;
}

// C3D8 TL
FORM_REGISTER_TL(double, 3, C3D8)
FORM_REGISTER_TL(float, 3, C3D8)
// // C3D8R TL
// FORM_REGISTER_TL(double, 3, C3D8R)
// FORM_REGISTER_TL(float, 3, C3D8R)
// C3D20 TL
FORM_REGISTER_TL(double, 3, C3D20)
FORM_REGISTER_TL(float, 3, C3D20)
// C3D20R TL
FORM_REGISTER_TL(double, 3, C3D20R)
FORM_REGISTER_TL(float, 3, C3D20R)

} // namespace fem
