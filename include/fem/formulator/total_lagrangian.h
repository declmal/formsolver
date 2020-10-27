#ifndef FEM_FORMULATOR_TOTAL_LAGRANGIAN_H_
#define FEM_FORMULATOR_TOTAL_LAGRANGIAN_H_

#include "formulator.h"

namespace fem {
FORM_REGISTER_FORMULATOR_TEMPLATE()
class TotalLagrangian : public Formulator<T,Dim,EType> {
  public:
    TotalLagrangian(EType<T>* elem_);
    void form_elem_stiff();
  private:
    T* X0;
    T* hbuf;
    T* weights;
    unsigned int num_ipoints;
    unsigned int num_nodes;
};

template <
  typename T,
  template <typename> class EType
>
using TL3D = TotalLagrangian<T,3,EType>;

#define FORM_REGISTER_TL(T, Dim, EType) \
  template \
  TotalLagrangian<T,Dim,EType>::TotalLagrangian(EType<T>* elem_); \
  template \
  void TotalLagrangian<T,Dim,EType>::form_elem_stiff();
} // namespace fem

#endif // FEM_FORMULATOR_TOTAL_LAGRANGIAN_H_
