#ifndef FEM_IPROP_ELEMENT_IPROP_H_
#define FEM_IPROP_ELEMENT_IPROP_H_

namespace fem {
template <typename T>
struct IPropInitializer {
  static inline void initialize(T* const hbuf, T* const weights);
};

#define FORM_IPROP_TEMPLATE() \
  template < \
    typename T, unsigned int NI, unsigned int Dim, unsigned int N, \
    template <typename> class IPropInitializerType \
  >

FORM_IPROP_TEMPLATE()
struct IProp {
  T hbuf[NI*N*Dim];
  T weights[NI];
  static unsigned int get_num_ipoints() { return NI; }
  static unsigned int get_num_nodes() { return N; }
  static unsigned int get_ndim() { return Dim; }
  T* get_hbuf() { return hbuf; }
  T* get_weights() { return weights; }
  IProp() {
    IPropInitializerType<T>::initialize(hbuf, weights);
  }
};
} // namespace fem

#endif // FEM_IPROP_ELEMENT_IPROP_H_
