#ifndef FEM_FE_ELEMENT_H_
#define FEM_FE_ELEMENT_H_

#include <common/gauss_legendre.h>

namespace fem {
template <typename T, unsigned int NI, unsigned int N>
struct ElementIProp {
  T hbuf[NI*N*3];
  T weights[NI];
  // unsigned int num_nodes;
  // constexpr ElementIProp() {
    // num_nodes = N;
  // }
};

template <typename T, unsigned int NI, unsigned int N>
class Element {
  public:
    void init_coordinate(const T* const data, const unsigned int size);
    void init_coordinate();
    virtual void form_elem_stiff_cpu(
      const ElementIProp<T,NI,N> iprop);
  private:
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

#define REGISTER_ELEMENT(T, NI, N) \
  template \
  void Element<T,NI,N>::form_elem_stiff_cpu( \
    const ElementIProp<T,NI,N> iprop); \
  template \
  void Element<T,NI,N>::init_coordinate( \
    const T* const data, const unsigned int size);
} // namespace fem

#endif // FEM_FE_ELEMENT_H_
