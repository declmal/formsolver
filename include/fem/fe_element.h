#ifndef FEM_FE_ELEMENT_H_
#define FEM_FE_ELEMENT_H_

#include <common/gauss_legendre.h>

namespace fem {
template <typename T, unsigned int NI, unsigned int N>
struct ElementVolIProp {
  T hbuf[NI*N*3];
  T weights[NI];
  // unsigned int num_nodes;
  // constexpr ElementVolIProp() {
    // num_nodes = N;
  // }
};

template <typename T, unsigned int NI, unsigned int N>
class ElementVol {
  public:
    virtual void form_elem_stiff(
      const T* const Ke, const ElementVolIProp<T,NI,N> iprop);
};
} // namespace fem

#endif // FEM_FE_ELEMENT_H_
