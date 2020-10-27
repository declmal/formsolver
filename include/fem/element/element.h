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
    void init_coordinate(const T* const data);
    unsigned int const get_id();
    unsigned int const get_ndim();
    T* const get_X0();
    T* const get_hbuf();
    T* const get_weights();
    unsigned int const get_num_ipoints();
    unsigned int const get_num_nodes();
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
    /*!
     * \brief Initial Global Coordinate Matrix, of shape (3, N)
     */
    T X0[Dim*N];
};

#define FORM_REGISTER_ELEMENT(T, N, Dim, IPropType) \
  template \
  void Element<T,N,Dim,IPropType>::init_coordinate( \
    const T* const data); \
  template \
  unsigned int const Element<T,N,Dim,IPropType>::get_id(); \
  template \
  unsigned int const Element<T,N,Dim,IPropType>::get_ndim(); \
  template \
  T* const Element<T,N,Dim,IPropType>::get_X0(); \
  template \
  T* const Element<T,N,Dim,IPropType>::get_hbuf(); \
  template \
  T* const Element<T,N,Dim,IPropType>::get_weights(); \
  template \
  unsigned int const Element<T,N,Dim,IPropType>::get_num_ipoints(); \
  template \
  unsigned int const Element<T,N,Dim,IPropType>::get_num_nodes();

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
