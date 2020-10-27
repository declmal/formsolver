#ifndef FEM_FORMULATOR_FORMULATOR_H_
#define FEM_FORMULATOR_FORMULATOR_H_

namespace fem {
#define FORM_REGISTER_FORMULATOR_TEMPLATE() \
  template < \
    typename T, unsigned int Dim, \
    template <typename> class EType \
  >

FORM_REGISTER_FORMULATOR_TEMPLATE()
class Formulator {
  public:
    Formulator(EType<T>* elem_) : elem(elem_) {}
    virtual int form_elem_stiff() = 0;
  protected:
    EType<T>* elem;
};
} // namespace fem

#endif // FEM_FORMULATOR_FORMULATOR_H_
