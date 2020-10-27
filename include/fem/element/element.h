#ifndef FEM_ELEMENT_ELEMENT_H_
#define FEM_ELEMENT_ELEMENT_H_

namespace fem {
#define FORM_REGISTER_ELEMENT_TEMPLATE() \
  template < \
    typename T, unsigned int N, unsigned int Dim, \
    template <typename> class IPropType>

FORM_REGISTER_ELEMENT_TEMPLATE()
class Element {
  public:
    void init_coordinate(const T* const data);
    void init_coordinate();
    unsigned int get_ndim();
    T* get_X0();
    T* get_hbuf();
    T* get_weights();
    unsigned int get_num_ipoints();
    unsigned int get_num_nodes();
  protected:
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
  unsigned int Element<T,N,Dim,IPropType>::get_ndim(); \
  template \
  T* Element<T,N,Dim,IPropType>::get_X0(); \
  template \
  T* Element<T,N,Dim,IPropType>::get_hbuf(); \
  template \
  T* Element<T,N,Dim,IPropType>::get_weights(); \
  template \
  unsigned int Element<T,N,Dim,IPropType>::get_num_ipoints(); \
  template \
  unsigned int Element<T,N,Dim,IPropType>::get_num_nodes();
} // namespace fem

#endif // FEM_ELEMENT_ELEMENT_H_