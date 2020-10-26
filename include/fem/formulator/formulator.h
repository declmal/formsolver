#ifndef FEM_FORMULATOR_H_
#define FEM_FORMULATOR_H_

namespace fem {
template <typename T, unsigned int Dim>
struct Formulator {
  virtual void form_elem_stiff() = 0;
};
} // namespace fem

#endif // FEM_FORMULATOR_H_
