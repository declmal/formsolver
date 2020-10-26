#ifndef FEM_ELEMENT_ELEMENT_H_
#define FEM_ELEMENT_ELEMENT_H_

namespace fem {
#define FORM_REGISTER_ELEMENT_TEMPLATE() \
  template < \
    typename T, unsigned int N, \
    template <typename> class IPropType, \
    template <typename> class FormType>

FORM_REGISTER_ELEMENT_TEMPLATE()
class Element {
  public:
    void init_coordinate(const T* const data, const unsigned int size);
    void init_coordinate();
    virtual void form_elem_stiff();
  private:
    static FormType<T> form;
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

#define FORM_REGISTER_ELEMENT(T, N, IPropType, FormType) \
  template \
  void Element<T,N,IPropType,FormType>::form_elem_stiff(); \
  template \
  void Element<T,N,IPropType,FormType>::init_coordinate( \
    const T* const data, const unsigned int size);
} // namespace fem

#endif // FEM_ELEMENT_ELEMENT_H_
