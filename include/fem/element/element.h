#ifndef FEM_ELEMENT_ELEMENT_H_
#define FEM_ELEMENT_ELEMENT_H_

namespace fem {
#define FORM_REGISTER_ELEMENT_TEMPLATE() \
  template < \
    typename T, unsigned int N, \
    template <typename> class IPropType>

FORM_REGISTER_ELEMENT_TEMPLATE()
class Element {
  public:
    void init_coordinate(const T* const data, const unsigned int size);
    void init_coordinate();
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
    T X0[3*N];
};

#define FORM_REGISTER_ELEMENT(T, N, IPropType) \
  template \
  void Element<T,N,IPropType>::init_coordinate( \
    const T* const data, const unsigned int size); \
  template \
  T* Element<T,N,IPropType>::get_X0(); \
  template \
  T* Element<T,N,IPropType>::get_hbuf(); \
  template \
  T* Element<T,N,IPropType>::get_weights(); \
  template \
  unsigned int Element<T,N,IPropType>::get_num_ipoints(); \
  template \
  unsigned int Element<T,N,IPropType>::get_num_nodes();
} // namespace fem

#endif // FEM_ELEMENT_ELEMENT_H_
