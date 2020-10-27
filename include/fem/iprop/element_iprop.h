#ifndef FEM_IPROP_ELEMENT_IPROP_H_
#define FEM_IPROP_ELEMENT_IPROP_H_

namespace fem {
template <typename T>
struct IPropInitializer {
  static inline void iprop_initialize(T* const hbuf, T* const weights);
};

#define FORM_IPROP_TEMPLATE() \
  template < \
    typename T, unsigned int NI, unsigned int Dim, unsigned int N, \
    template <typename> class IPropInitializerType \
  >

FORM_IPROP_TEMPLATE()
struct IProp {
  static T hbuf[NI*N*Dim];
  static T weights[NI];
  static unsigned int get_num_ipoints() { return NI; }
  static unsigned int get_num_nodes() { return N; }
  static unsigned int get_ndim() { return Dim; }
  static T* get_hbuf() { return hbuf; }
  static T* get_weights() { return weights; }
  static void initialize() {
    IPropInitializerType<T>::iprop_initialize(hbuf, weights);
  }
};

FORM_IPROP_TEMPLATE()
T IProp<T,NI,Dim,N,IPropInitializerType>::hbuf[NI*N*Dim];

FORM_IPROP_TEMPLATE()
T IProp<T,NI,Dim,N,IPropInitializerType>::weights[NI];

} // namespace fem

#endif // FEM_IPROP_ELEMENT_IPROP_H_
