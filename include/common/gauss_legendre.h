#ifndef GAUSS_LEGENDRE_H_
#define GAUSS_LEGENDRE_H_

#include <math.h>

template <typename T, unsigned int N>
struct Legendre {
  static inline T impl(T x) {
    return (
      (T)(2*N-1) * x * Legendre<T,N-1>::impl(x) -
      (T)(N-1) * Legendre<T,N-2>::impl(x)
    ) / (T)N;
  }
}; 
template <typename T>
struct Legendre<T, 1> {
  static inline T impl(T x) {
    return x;
  }
};
template <typename T>
struct Legendre<T, 0> {
  static inline T impl(T x) {
    return (T)1;
  }
};

template <typename T, unsigned int N>
T dlegendre(T x) {
  return (T)N/(pow(x,2)-(T)1) * (
    x*Legendre<T,N>::impl(x) - Legendre<T,N-1>::impl(x)
  );
}

// template <typename T, unsigned int N>
// void legendre_roots_tmp(T* const roots, double tol=1e-20) {
  // unsigned int b = ((N-1)>>1) + 1;
  // unsigned int M = N>>1;
  // roots[b-1] = (T)0;
  // for (unsigned int i = 1; i <= M; ++i) {
    // T x = (T)cos(
      // M_PI * (T)(i-0.25) / (T)(N+0.5)
    // );
    // auto error = 10 * tol;
    // unsigned int iters = 0;
    // while ((error > tol) && (iters < 1000)) {
      // auto dx = -Legendre<T,N>::impl(x) / dlegendre<T,N>(x);
      // x += dx;
      // iters += 1;
      // error = (double)abs(dx);
    // }
    // roots[i-1] = -x;
    // roots[M-i+b] = x;
  // }
// }

// template <typename T, unsigned int N>
// void gauss_interp(T* const roots, T* const weights, double tol=1e-20) {
  // legendre_roots_tmp<T, N>(roots, tol);
  // for (unsigned int i = 0; i < N; ++i) {
    // weights[i] = (T)2 / (
      // (T)(1-pow(roots[i],2)) * pow(dlegendre<T,N>(roots[i]),2)
    // );
  // }
// }

template <typename T, unsigned int N>
void legendre_roots(T roots[], const double tol=1e-20) {
  unsigned int b = ((N-1)>>1) + 1;
  unsigned int M = N>>1;
  roots[b-1] = (T)0;
  for (unsigned int i = 1; i <= M; ++i) {
    T x = (T)cos(
      M_PI * (T)(i-0.25) / (T)(N+0.5)
    );
    auto error = 10 * tol;
    unsigned int iters = 0;
    while ((error > tol) && (iters < 1000)) {
      auto dx = -Legendre<T,N>::impl(x) / dlegendre<T,N>(x);
      x += dx;
      iters += 1;
      error = (double)abs(dx);
    }
    roots[i-1] = -x;
    roots[M-i+b] = x;
  }
}

template <typename T, unsigned int N>
struct Gauss1D {
  T roots[N];
  T weights[N];
  constexpr static double tol = 1e-20;
  constexpr Gauss1D() {
    legendre_roots<T, N>(roots, tol);
    for (unsigned int i = 0; i < N; ++i) {
      weights[i] = (T)2 / (
        (T)(1-pow(roots[i],2)) * pow(dlegendre<T,N>(roots[i]),2)
      );
    }
  }
};

template <typename T, unsigned int N0, unsigned int N1, unsigned int N2>
struct Gauss3D {
  /*!
   * \brief number of gauss points
   */
  constexpr static auto N = N0 * N1 * N2;
  /*!
   * \brief r interpolation points, of shape (N, 3)
   */
  T r[N*3];
  /*!
   * \brief w interpolation weights, of shape (N,)
   */
  T w[N];
  constexpr Gauss3D() {
    Gauss1D<T, N0> g0;
    Gauss1D<T, N1> g1;
    Gauss1D<T, N2> g2;
    auto stride0 = N1 * N2;
    auto stride1 = N2;
    for (unsigned int i = 0; i < N0; ++i) {
      for (unsigned int j = 0; j < N1; ++j) {
        for (unsigned int k = 0; k < N2; ++k) {
          auto ind = i*stride0 + j*stride1 + k;
          auto ind3 = ind * 3;
          r[ind3] = g0.roots[i];
          r[ind3+1] = g1.roots[j];
          r[ind3+2] = g2.roots[k];
          w[ind] = g0.weights[i] * g1.weights[j] * g2.weights[k];
        }
      }
    }
  }
};

#endif // GAUSS_LEGENFRE_H_
