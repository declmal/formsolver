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

template <typename T, unsigned int N>
void legendre_roots(T* const roots, double tol=1e-20) {
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
void gauss_interp(T* const roots, T* const weights, double tol=1e-20) {
  legendre_roots<T, N>(roots, tol);
  for (unsigned int i = 0; i < N; ++i) {
    weights[i] = (T)2 / (
      (T)(1-pow(roots[i],2)) * pow(dlegendre<T,N>(roots[i]),2)
    );
  }
}

#endif // GAUSS_LEGENFRE_H_
