#ifndef FEM_ELEMENT_ELEMENT_H_
#define FEM_ELEMENT_ELEMENT_H_

#include <fem/iprop/brick_iprop.h>

namespace fem {
#define FORM_REGISTER_ELEMENT_TEMPLATE() \
  template < \
    typename T, unsigned int N, unsigned int Dim, \
    template <typename> class IPropType>

FORM_REGISTER_ELEMENT_TEMPLATE()
class Element {
  public:
    Element(unsigned int id_) : id(id_) {}
    unsigned int const get_id();
    T* const get_hbuf();
    T* const get_weights();
    unsigned int const get_num_ipoints();
  protected:
    /*!
     * \brief Element Id
     */
    unsigned int id;
    /*!
     * \brief Element Interpolation Property
     */
    static IPropType<T> iprop;
    /*!
     * \brief Global Node Id List, of shape (N,)
     */
    unsigned int IdN[N];
};

#define FORM_REGISTER_ELEMENT(T, N, Dim, IPropType) \
  template \
  unsigned int const Element<T,N,Dim,IPropType>::get_id(); \
  template \
  T* const Element<T,N,Dim,IPropType>::get_hbuf(); \
  template \
  T* const Element<T,N,Dim,IPropType>::get_weights(); \
  template \
  unsigned int const Element<T,N,Dim,IPropType>::get_num_ipoints();

template <typename T>
using C3D8 = Element<T,8,3,C3D8IProp>;
// template <typename T>
// using C3D8R = Element<T,8,3,C3D8RIProp>;
template <typename T>
using C3D20 = Element<T,20,3,C3D20IProp>;
template <typename T>
using C3D20R = Element<T,20,3,C3D20RIProp>;

} // namespace fem

#endif // FEM_ELEMENT_ELEMENT_H_
