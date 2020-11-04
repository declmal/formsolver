#ifndef FEM_ELEMENT_BRICK_H_
#define FEM_ELEMENT_BRICK_H_

#include <fem/formulator/total_lagrangian.h>
#include <fem/iprop/brick_iprop.h>

namespace fem {
#define FORM_REGISTER_BRICK_ELEMENT_TL(BrickType, N0, N1, N2, N, Dim) \
  FORM_REGISTER_ELEM_FORM(BrickType, N0*N1*N2, N, Dim) \
  FORM_REGISTER_BRICK_IPROP(BrickType, N0, N1, N2, N)

FORM_REGISTER_BRICK_ELEMENT_TL(C3D8, 2, 2, 2, 8, 3)
// FORM_REGISTER_BRICK_ELEMENT_TL(C3D8R, 1, 1, 1, 8, 3)
FORM_REGISTER_BRICK_ELEMENT_TL(C3D20, 3, 3, 3, 20, 3)
FORM_REGISTER_BRICK_ELEMENT_TL(C3D20R, 2, 2, 2, 20, 3)
FORM_REGISTER_BRICK_ELEMENT_TL(C3D8I, 2, 2, 2, 11, 3)
} // namespace fem

#endif // FEM_ELEMENT_BRICK_H_
