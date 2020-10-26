#ifndef FEM_TOTAL_LAGRANGIAN_H_
#define FEM_TOTAL_LAGRANGIAN_H_

#include "formulator.h"

namespace fem {
template <typename T, unsigned int Dim>
struct TotalLagrangian : Formulator<T,Dim> {
  void form_elem_stiff();
};

template <typename T>
using TL3D = TotalLagrangian<T,3>;
} // namespace fem

#endif // FEM_TOTAL_LAGRANGIAN_H_
