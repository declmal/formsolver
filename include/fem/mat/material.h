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
  T C[NR*NR];
  static unsigned int get_nrows() { return NR; }
  static unsigned int get_num_params() { return NP; }
  T* get_C() { return C; }
  Mat(T* p) {
    MatIinitializerType<T>::initialize(C, p);
  }
};
} // namespace

#endif // FEM_MAT_MATERIAL_H_
