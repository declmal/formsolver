#ifndef FEM_FORMULATOR_TOTAL_LAGRANGIAN_H_
#define FEM_FORMULATOR_TOTAL_LAGRANGIAN_H_

#include "formulator.h"

namespace fem {
FORM_REGISTER_FORMULATOR_TEMPLATE()
class TotalLagrangian : public Formulator<T,Dim,EType> {
  public:
    TotalLagrangian(EType<T>* elem_) : Formulator<T,Dim,EType>(elem_) {}
    void load_elem_data();
    void form_elem_stiff();
  private:
    T* X0;
    T* J;
    T* Ke;
};

template <
  typename T,
  template <typename> class EType
>
using TL3D = TotalLagrangian<T,3,EType>;

#define FORM_REGISTER_TL(T, Dim, EType) \
  template \
  void TotalLagrangian<T,Dim,EType>::load_elem_data(); \
  template \
  void TotalLagrangian<T,Dim,EType>::form_elem_stiff();
} // namespace fem

#endif // FEM_FORMULATOR_TOTAL_LAGRANGIAN_H_
