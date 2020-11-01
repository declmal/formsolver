#ifndef GAUSS_LEGENDRE_H_
#define GAUSS_LEGENDRE_H_

#include <math.h>

template <typename T, unsigned int N>
struct Legendre {
  static inline T legendre(T x) {
    return (
      (T)(2*N-1) * x * Legendre<T,N-1>::legendre(x) -
      (T)(N-1) * Legendre<T,N-2>::legendre(x)
    ) / (T)N;
  }
}; 
template <typename T>
struct Legendre<T, 1> {
  static inline T legendre(T x) {
    return x;
  }
};
template <typename T>
struct Legendre<T, 0> {
  static inline T legendre(T x) {
    return (T)1;
  }
};

template <typename T, unsigned int N>
T dlegendre(T x) {
  return (T)N/(pow(x,2)-(T)1) * (
    x*Legendre<T,N>::legendre(x) - Legendre<T,N-1>::legendre(x)
  );
}

template <typename T, unsigned int N>
void legendre_roots(
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

template <
  typename T, unsigned Dim, unsigned int N0, unsigned int N1, unsigned int N2>
struct GaussRoots {
};
template <typename T, unsigned int N0>
struct GaussRoots<T,1,N0,1,1> {
  constexpr static double tol = 1e-20;
  T roots[N0];
  constexpr GaussRoots() {
    legendre_roots<T,N0>(roots, tol);
  }
};
template <typename T, unsigned int N0, unsigned int N1, unsigned int N2>
struct GaussRoots<T,3,N0,N1,N2> {
  constexpr static auto NI = N0 * N1 * N2;
  T roots[NI*3];
  constexpr GaussRoots() {
    GaussRoots<T,1,N0,1,1> gr0;
    GaussRoots<T,1,N1,1,1> gr1;
    GaussRoots<T,1,N2,1,1> gr2;
    auto stride0 = N0 * N1;
    auto stride1 = N1;
    for (unsigned int i = 0; i < N2; ++i) {
      for (unsigned int j = 0; j < N1; ++j) {
        for (unsigned int k = 0; k < N0; ++k) {
          auto ind = i*stride0 + j*stride1 + k;
          auto ind3 = ind * 3;
          roots[ind3] = gr0.roots[k];
          roots[ind3+1] = gr1.roots[j];
          roots[ind3+2] = gr2.roots[i];
        }
      }
    }
  }
};

template <
  typename T, unsigned int Dim, unsigned int N0, unsigned N1, unsigned int N2>
struct GaussWeights {
};
template <typename T, unsigned int N0>
struct GaussWeights<T,1,N0,1,1>  {
  T weights[N0];
  constexpr GaussWeights() {
    GaussRoots<T,1,N0,1,1> gr;
    for (unsigned int i = 0; i < N0; ++i) {
      weights[i] = (T)2 / (
        (T)(1-pow(gr.roots[i],2)) * pow(dlegendre<T,N0>(gr.roots[i]),2)
      );
    }
  }
};
template <typename T, unsigned int N0, unsigned int N1, unsigned int N2>
struct GaussWeights<T,3,N0,N1,N2> {
  constexpr static auto NI = N0 * N1 * N2;
  T weights[NI];
  constexpr GaussWeights() {
    GaussWeights<T,1,N0,1,1> gw0;
    GaussWeights<T,1,N1,1,1> gw1;
    GaussWeights<T,1,N2,1,1> gw2;
    auto stride0 = N0 * N1;
    auto stride1 = N0;
    for (unsigned int i = 0; i < N2; ++i) {
      for (unsigned int j = 0; j < N1; ++j) {
        for (unsigned int k = 0; k < N0; ++k) {
          auto ind = i*stride0 + j*stride1 + k;
          weights[ind] = gw0.weights[k] * gw1.weights[j] * gw2.weights[i];
        }
      }
    }
  }
};

#endif // GAUSS_LEGENFRE_H_
