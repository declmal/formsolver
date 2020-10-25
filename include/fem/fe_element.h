#ifndef FEM_FE_ELEMENT_H_
#define FEM_FE_ELEMENT_H_

#include <common/gauss_legendre.h>

namespace fem {
template <typename T, unsigned int NI, unsigned int N>
struct ElementVolInterpProp {
  T hbuf[NI*N*3];
  T weights[NI];
  // unsigned int num_nodes;
  // constexpr ElementVolInterpProp() {
    // num_nodes = N;
  // }
};

template <typename T, unsigned int NI, unsigned int N> class ElementVol {
  public:
    static ElementVolInterpProp<T, NI, N> iprop;
    virtual void form_elem_stiff(const T* const Ke);
};
} // namespace fem

#endif // FEM_FE_ELEMENT_H_
