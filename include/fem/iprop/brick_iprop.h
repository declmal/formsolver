#ifndef FEM_IPROP_BRICK_IPROP_H_
#define FEM_IPROP_BRICK_IPROP_H_

#include <common/gauss_legendre.h>
#include "element_iprop.h"

namespace fem {
template <typename T, int R>
struct BrickInterpFact {
  static inline T compute(T r);
  static inline T compute_deriv_1st(T r);
};
template <typename T>
struct BrickInterpFact<T,1> {
  static inline T compute(T r) {
    return (T)(0.5*(1+r));
  }
  static inline T compute_deriv_1st(T r) {
    return (T)0.5;
  }
};
template <typename T>
struct BrickInterpFact<T,-1> {
  static inline T compute(T r) {
    return (T)(0.5*(1-r));
  }
  static inline T compute_deriv_1st(T r) {
    return (T)-0.5;
  }
};
template <typename T>
struct BrickInterpFact<T,0> {
  static inline T compute(T r) {
    return (T)(1-r*r);
  }
  static inline T compute_deriv_1st(T r) {
    return (T)(-2*r);
  }
};

template <unsigned int I>
struct BrickInterpCoord {
  const static int R0;
  const static int R1;
  const static int R2;
};

#define FORM_REGISTER_BRICKINTERPCOORD(i, x, y, z) \
  template <> \
  struct BrickInterpCoord<i> { \
    const static int R0 = x; \
    const static int R1 = y; \
    const static int R2 = z; \
  };

FORM_REGISTER_BRICKINTERPCOORD(0, 1, 1, 1)
FORM_REGISTER_BRICKINTERPCOORD(1, -1, 1, 1)
FORM_REGISTER_BRICKINTERPCOORD(2, -1, -1, 1)
FORM_REGISTER_BRICKINTERPCOORD(3, 1, -1, 1)
FORM_REGISTER_BRICKINTERPCOORD(4, 1, 1, -1)
FORM_REGISTER_BRICKINTERPCOORD(5, -1, 1, -1)
FORM_REGISTER_BRICKINTERPCOORD(6, -1, -1, -1)
FORM_REGISTER_BRICKINTERPCOORD(7, 1, -1, -1)
FORM_REGISTER_BRICKINTERPCOORD(8, 0, 1, 1)
FORM_REGISTER_BRICKINTERPCOORD(9, -1, 0, 1)
FORM_REGISTER_BRICKINTERPCOORD(10, 0, -1, 1)
FORM_REGISTER_BRICKINTERPCOORD(11, 1, 0, 1)
FORM_REGISTER_BRICKINTERPCOORD(12, 0, 1, -1)
FORM_REGISTER_BRICKINTERPCOORD(13, -1, 0, -1)
FORM_REGISTER_BRICKINTERPCOORD(14, 0, -1, -1)
FORM_REGISTER_BRICKINTERPCOORD(15, 1, 0, -1)
FORM_REGISTER_BRICKINTERPCOORD(16, 0, 1, 0)
FORM_REGISTER_BRICKINTERPCOORD(17, -1, 0, 0)
FORM_REGISTER_BRICKINTERPCOORD(18, 0, -1, 0)
FORM_REGISTER_BRICKINTERPCOORD(19, 1, 0, 0)

template <typename T, unsigned int I>
struct BrickInterpMono {
  static inline T compute(T r0, T r1, T r2) {
    auto G0 = BrickInterpFact<T,BrickInterpCoord<I>::R0>::compute(r0);
    auto G1 = BrickInterpFact<T,BrickInterpCoord<I>::R1>::compute(r1);
    auto G2 = BrickInterpFact<T,BrickInterpCoord<I>::R2>::compute(r2);
    return G0 * G1 * G2;
  }
  static inline T compute_d1_r0(T r0, T r1, T r2) {
    auto G1 = BrickInterpFact<T,BrickInterpCoord<I>::R1>::compute(r1);
    auto G2 = BrickInterpFact<T,BrickInterpCoord<I>::R2>::compute(r2);
    auto dG0 = BrickInterpFact<
      T, BrickInterpCoord<I>::R0>::compute_deriv_1st(r0);
    return dG0 * G1 * G2;
  }
  static inline T compute_d1_r1(T r0, T r1, T r2) {
    auto G0 = BrickInterpFact<T,BrickInterpCoord<I>::R0>::compute(r0);
    auto G2 = BrickInterpFact<T,BrickInterpCoord<I>::R2>::compute(r2);
    auto dG1 = BrickInterpFact<
      T,BrickInterpCoord<I>::R1>::compute_deriv_1st(r1);
    return G0 * dG1 * G2;
  }
  static inline T compute_d1_r2(T r0, T r1, T r2) {
    auto G0 = BrickInterpFact<T,BrickInterpCoord<I>::R0>::compute(r0);
    auto G1 = BrickInterpFact<T,BrickInterpCoord<I>::R1>::compute(r1);
    auto dG2 = BrickInterpFact<
      T, BrickInterpCoord<I>::R2>::compute_deriv_1st(r2);
    return G0 * G1 * dG2;
  }
};

template <typename T, unsigned int I, unsigned int N>
struct BrickInterpPoly {
  static inline T compute(T r0, T r1, T r2);
  static inline T compute_d1_r0(T r0, T r1, T r2);
  static inline T compute_d1_r1(T r0, T r1, T r2);
  static inline T compute_d1_r2(T r0, T r1, T r2);
};
template <typename T, unsigned int I>
struct BrickInterpPoly<T,I,8> {
  static inline T compute(T r0, T r1, T r2) {
    return BrickInterpMono<T,I>::compute(r0, r1, r2);
  }
  static inline T compute_d1_r0(T r0, T r1, T r2) {
    return BrickInterpMono<T,I>::compute_d1_r0(r0, r1, r2);
  }
  static inline T compute_d1_r1(T r0, T r1, T r2) {
    return BrickInterpMono<T,I>::compute_d1_r1(r0, r1, r2);
  }
  static inline T compute_d1_r2(T r0, T r1, T r2) {
    return BrickInterpMono<T,I>::compute_d1_r2(r0, r1, r2);
  }
};

#define FORM_REGISTER_BRICKINTERPPOLY_N20_LOWER(i, j, k, l) \
  template <typename T> \
  struct BrickInterpPoly<T,i,20> { \
    static inline T compute(T r0, T r1, T r2) { \
      auto g##i = BrickInterpMono<T,i>::compute(r0, r1, r2); \
      auto g##j = BrickInterpMono<T,j>::compute(r0, r1, r2); \
      auto g##k = BrickInterpMono<T,k>::compute(r0, r1, r2); \
      auto g##l = BrickInterpMono<T,l>::compute(r0, r1, r2); \
      return (T)(g##i-(g##j+g##k+g##l)/2); \
    } \
    static inline T compute_d1_r0(T r0, T r1, T r2) { \
      auto dg##i = BrickInterpMono<T,i>::compute_d1_r0(r0, r1, r2); \
      auto dg##j = BrickInterpMono<T,j>::compute_d1_r0(r0, r1, r2); \
      auto dg##k = BrickInterpMono<T,k>::compute_d1_r0(r0, r1, r2); \
      auto dg##l = BrickInterpMono<T,l>::compute_d1_r0(r0, r1, r2); \
      return (T)(dg##i-(dg##j+dg##k+dg##l)/2); \
    } \
    static inline T compute_d1_r1(T r0, T r1, T r2) { \
      auto dg##i = BrickInterpMono<T,i>::compute_d1_r1(r0, r1, r2); \
      auto dg##j = BrickInterpMono<T,j>::compute_d1_r1(r0, r1, r2); \
      auto dg##k = BrickInterpMono<T,k>::compute_d1_r1(r0, r1, r2); \
      auto dg##l = BrickInterpMono<T,l>::compute_d1_r1(r0, r1, r2); \
      return (T)(dg##i-(dg##j+dg##k+dg##l)/2); \
    } \
    static inline T compute_d1_r2(T r0, T r1, T r2) { \
      auto dg##i = BrickInterpMono<T,i>::compute_d1_r2(r0, r1, r2); \
      auto dg##j = BrickInterpMono<T,j>::compute_d1_r2(r0, r1, r2); \
      auto dg##k = BrickInterpMono<T,k>::compute_d1_r2(r0, r1, r2); \
      auto dg##l = BrickInterpMono<T,l>::compute_d1_r2(r0, r1, r2); \
      return (T)(dg##i-(dg##j+dg##k+dg##l)/2); \
    } \
  };

FORM_REGISTER_BRICKINTERPPOLY_N20_LOWER(0, 8, 11, 16)
FORM_REGISTER_BRICKINTERPPOLY_N20_LOWER(1, 8, 9, 17)
FORM_REGISTER_BRICKINTERPPOLY_N20_LOWER(2, 9, 10, 18)
FORM_REGISTER_BRICKINTERPPOLY_N20_LOWER(3, 10, 11, 19)
FORM_REGISTER_BRICKINTERPPOLY_N20_LOWER(4, 12, 15, 16)
FORM_REGISTER_BRICKINTERPPOLY_N20_LOWER(5, 12, 13, 17)
FORM_REGISTER_BRICKINTERPPOLY_N20_LOWER(6, 13, 14, 18)
FORM_REGISTER_BRICKINTERPPOLY_N20_LOWER(7, 14, 15, 19)

#define FORM_REGISTER_BRICKINTERPPOLY_N20_HIGHER(i) \
template <typename T> \
  struct BrickInterpPoly<T,i,20> { \
  static inline T compute(T r0, T r1, T r2) { \
    return BrickInterpMono<T,i>::compute(r0, r1, r2); \
  } \
  static inline T compute_d1_r0(T r0, T r1, T r2) { \
    return BrickInterpMono<T,i>::compute_d1_r0(r0, r1, r2); \
  } \
  static inline T compute_d1_r1(T r0, T r1, T r2) { \
    return BrickInterpMono<T,i>::compute_d1_r1(r0, r1, r2); \
  } \
  static inline T compute_d1_r2(T r0, T r1, T r2) { \
    return BrickInterpMono<T,i>::compute_d1_r2(r0, r1, r2); \
  } \
};

FORM_REGISTER_BRICKINTERPPOLY_N20_HIGHER(8)
FORM_REGISTER_BRICKINTERPPOLY_N20_HIGHER(9)
FORM_REGISTER_BRICKINTERPPOLY_N20_HIGHER(10)
FORM_REGISTER_BRICKINTERPPOLY_N20_HIGHER(11)
FORM_REGISTER_BRICKINTERPPOLY_N20_HIGHER(12)
FORM_REGISTER_BRICKINTERPPOLY_N20_HIGHER(13)
FORM_REGISTER_BRICKINTERPPOLY_N20_HIGHER(14)
FORM_REGISTER_BRICKINTERPPOLY_N20_HIGHER(15)
FORM_REGISTER_BRICKINTERPPOLY_N20_HIGHER(16)
FORM_REGISTER_BRICKINTERPPOLY_N20_HIGHER(17)
FORM_REGISTER_BRICKINTERPPOLY_N20_HIGHER(18)
FORM_REGISTER_BRICKINTERPPOLY_N20_HIGHER(19)

template <typename T, unsigned int I, unsigned int N> struct BrickInterpDeriv {
  static inline void interp_deriv(T* const h, T r0, T r1, T r2) {
    BrickInterpDeriv<T,I-1,N>::interp_deriv(h, r0, r1, r2);
    auto ind = 3 * I;
    h[ind] = BrickInterpPoly<T,I,N>::compute_d1_r0(r0, r1, r2);
    h[ind+1] = BrickInterpPoly<T,I,N>::compute_d1_r1(r0, r1, r2);
    h[ind+2] = BrickInterpPoly<T,I,N>::compute_d1_r2(r0, r1, r2);
  }
};
template <typename T, unsigned int N> struct BrickInterpDeriv<T,0,N> {
  static inline void interp_deriv(T* const h, T r0, T r1, T r2) {
    h[0] = BrickInterpPoly<T,0,N>::compute_d1_r0(r0, r1, r2);
    h[1] = BrickInterpPoly<T,0,N>::compute_d1_r1(r0, r1, r2);
    h[2] = BrickInterpPoly<T,0,N>::compute_d1_r2(r0, r1, r2);
  } 
};

template <typename T, unsigned int N> void brick_interp_deriv(
  T* const h, T r0, T r1, T r2) {
  BrickInterpDeriv<T,N-1,N>::interp_deriv(h, r0, r1, r2);
}

template <
  typename T, unsigned int N0, unsigned int N1, unsigned int N2, unsigned int N> 
struct BrickIProp : ElementIProp<T,N0*N1*N2,N> {
  constexpr BrickIProp() : ElementIProp<T,N0*N1*N2,N>() {
    // init hbuf
    auto NI = N0 * N1 * N2;
    GaussRoots<T,3,N0,N1,N2> gr;
    T* r = gr.roots;
    T* h = this->hbuf;
    for (unsigned int i = 0; i < NI; ++i) {
      brick_interp_deriv<T, N>(h, r[0], r[1], r[2]);
      r += 3;
      h += NI;
    }
    // init weights
    GaussWeights<T,3,N0,N1,N2> gw;
    for (unsigned int i = 0; i < NI; ++i) {
      this->weights[i] = gw.weights[i];
    }
  }
};

template <typename T>
using C3D8IProp = BrickIProp<T,2,2,2,8>;
template <typename T>
using C3D8RIProp = BrickIProp<T,1,1,1,8>;
template <typename T>
using C3D20IProp = BrickIProp<T,3,3,3,20>;
template <typename T>
using C3D20RIProp = BrickIProp<T,2,2,2,20>;
} // namespace fem

#endif // FEM_IPROP_BRICK_IPROP_H_
