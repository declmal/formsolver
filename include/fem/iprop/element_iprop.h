#ifndef FEM_IPROP_ELEMENT_IPROP_H_
#define FEM_IPROP_ELEMENT_IPROP_H_

namespace fem {
template <typename T, unsigned int NI, unsigned int N>
class ElementIProp {
  public:
    T get_num_ipoints() { return NI; }
    T* get_hbuf() { return hbuf; }
    T* get_weights() { return weights; }
  protected:
    T hbuf[NI*N*3];
    T weights[NI];
};
} // namespace fem

#endif // FEM_IPROP_ELEMENT_IPROP_H_
