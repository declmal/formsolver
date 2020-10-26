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
    virtual void form_elem_stiff();
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
    /*!
     * \brief Jacobi Matrix, of shape (3, 3)
     */
    T J[9];
    /*!
     * \brief Element Stiffness Matrix, of shape (3N, 3N)
     */
    T Ke[9*N*N];
};

#define FORM_REGISTER_ELEMENT(T, N, IPropType) \
  template \
  void Element<T,N,IPropType>::form_elem_stiff(); \
  template \
  void Element<T,N,IPropType>::init_coordinate( \
    const T* const data, const unsigned int size);
} // namespace fem

#endif // FEM_ELEMENT_ELEMENT_H_
