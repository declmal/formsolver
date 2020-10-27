#ifndef FEM_ELEMENT_BRICK_H_
#define FEM_ELEMENT_BRICK_H_

#include <fem/iprop/brick_iprop.h>
#include "element.h"

namespace fem {
template <typename T>
class C3D8 : public Element<T,8,3,C3D8IProp> {
};

// template <typename T>
// class C3D8R : public Element<T,8,3,C3D8RIProp> {
// };

template <typename T>
class C3D20 : public Element<T,20,3,C3D20IProp> {
};

template <typename T>
class C3D20R : public Element<T,20,3,C3D20RIProp> {
};
} // namespace fem

#endif // FEM_ELEMENT_BRICK_H_
