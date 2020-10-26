#ifndef FEM_IPROP_ELEMENT_IPROP_H_
#define FEM_IPROP_ELEMENT_IPROP_H_

namespace fem {
template <typename T, unsigned int NI, unsigned int N>
struct ElementIProp {
  T hbuf[NI*N*3];
  T weights[NI];
  unsigned int num_ipoints;
  constexpr ElementIProp() {
    num_ipoints = NI;
  }
};
} // namespace fem

#endif // FEM_IPROP_ELEMENT_IPROP_H_
