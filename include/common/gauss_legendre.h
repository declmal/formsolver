#ifndef GAUSS_LEGENDRE_H_
#define GAUSS_LEGENDRE_H_

#include <math.h>

template <typename T, unsigned int N> struct Legendre {
  static inline T legendre(T x) {
    return (
      (T)(2*N-1) * x * Legendre<T,N-1>::legendre(x) -
      (T)(N-1) * Legendre<T,N-2>::legendre(x)
    ) / (T)N;
  }
}; 
template <typename T> struct Legendre<T, 1> {
  static inline T legendre(T x) {
    return x;
  }
};
template <typename T> struct Legendre<T, 0> {
  static inline T legendre(T x) {
    return (T)1;
  }
};

template <typename T, unsigned int N> T dlegendre(T x) {
  return (T)N/(pow(x,2)-(T)1) * (
    x*Legendre<T,N>::legendre(x) - Legendre<T,N-1>::legendre(x)
  );
}

template <typename T, unsigned int N> void legendre_roots(
  T roots[], const double tol=1e-20) {
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
      auto dx = -Legendre<T,N>::legendre(x) / dlegendre<T,N>(x);
      x += dx;
      iters += 1;
      error = (double)abs(dx);
    }
    roots[i-1] = -x;
    roots[M-i+b] = x;
  }
}

template <typename T, unsigned int N> struct GaussRoots1D {
  T roots[N];
  constexpr static double tol = 1e-20;
  constexpr GaussRoots1D() {
    legendre_roots<T, N>(roots, tol);
  }
};

template <typename T, unsigned int N> struct GaussWeights1D {
  T weights[N];
  constexpr GaussWeights1D() {
    GaussRoots1D<T, N> gr;
    for (unsigned int i = 0; i < N; ++i) {
      weights[i] = (T)2 / (
        (T)(1-pow(gr.roots[i],2)) * pow(dlegendre<T,N>(gr.roots[i]),2)
      );
    }
  }
};

template <typename T, unsigned int N0, unsigned int N1, unsigned int N2>
struct GaussRoots3D {
  constexpr static auto NI = N0 * N1 * N2;
  T roots[NI*3];
  constexpr GaussRoots3D() {
    GaussRoots1D<T, N0> gr0;
    GaussRoots1D<T, N1> gr1;
    GaussRoots1D<T, N2> gr2;
    auto stride0 = N1 * N2;
    auto stride1 = N2;
    for (unsigned int i = 0; i < N0; ++i) {
      for (unsigned int j = 0; j < N1; ++j) {
        for (unsigned int k = 0; k < N2; ++k) {
          auto ind = i*stride0 + j*stride1 + k;
          auto ind3 = ind * 3;
          roots[ind3] = gr0.roots[i];
          roots[ind3+1] = gr1.roots[j];
          roots[ind3+2] = gr2.roots[k];
        }
      }
    }
  }
};

template <typename T, unsigned int N0, unsigned int N1, unsigned int N2>
struct GaussWeights3D {
  constexpr static auto NI = N0 * N1 * N2;
  T weights[NI];
  constexpr GaussWeights3D() {
    GaussWeights1D<T, N0> gw0;
    GaussWeights1D<T, N1> gw1;
    GaussWeights1D<T, N2> gw2;
    auto stride0 = N1 * N2;
    auto stride1 = N2;
    for (unsigned int i = 0; i < N0; ++i) {
      for (unsigned int j = 0; j < N1; ++j) {
        for (unsigned int k = 0; k < N2; ++k) {
          auto ind = i*stride0 + j*stride1 + k;
          auto ind3 = ind * 3;
          weights[ind] = gw0.weights[i] * gw1.weights[j] * gw2.weights[k];
        }
      }
    }
  }
};

#endif // GAUSS_LEGENFRE_H_
