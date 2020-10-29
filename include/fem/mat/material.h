#ifndef FEM_MAT_MATERIAL_H_
#define FEM_MAT_MATERIAL_H_

namespace fem {
template <typename T>
struct MatInitializer {
  static inline void initialize(T* const C, T* const p);
};

#define FORM_MAT_TEMPLATE() \
  template < \
    typename T, unsigned int NR, unsigned int NP, \
    template <typename> class MatIinitializerType \
  >

FORM_MAT_TEMPLATE()
struct Mat {
  static T C[NR*NR];
  static unsigned int get_nrows() { return NR; }
  static unsigned int get_num_params() { return NP; }
  static T* get_C() { return C; }
  static inline void initialize(T* p) {
    MatIinitializerType<T>::initialize(C, p);
  }
};

FORM_MAT_TEMPLATE()
T Mat<T,NR,NP,MatIinitializerType>::C[NR*NR];
} // namespace

#endif // FEM_MAT_MATERIAL_H_
