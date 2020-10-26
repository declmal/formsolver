#ifndef FEM_ELEMENT_BRICK_H_
#define FEM_ELEMENT_BRICK_H_

#include <fem/iprop/brick_iprop.h>
#include <fem/element/element.h>

namespace fem {
template <
  typename T, template <typename> class FormType>
class C3D8 : public Element<T,8,C3D8IProp,FormType> {
};

// template <
  // typename T, template <typename> class FormType>
// class C3D8R : public Element<T,8,C3D8RIProp,FormType> {
// };

template <
  typename T, template <typename> class FormType>
class C3D20 : public Element<T,20,C3D20IProp,FormType> {
};

template <
  typename T, template <typename> class FormType>
class C3D20R : public Element<T,20,C3D20RIProp,FormType> {
};
} // namespace fem

#endif // FEM_ELEMENT_BRICK_H_
