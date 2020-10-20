template <typename T, unsigned int x>
struct JacobiSum {
  static inline void impl(const T* Ja, const T* hn, T* BL0) {
    JacobiSum<T, x-1>::impl(Ja, hn, BL0) + Ja[x-1]*hn[x-1];
  }
};

template <typename T>
struct JacobiSum<T, 0> {
  static inline void impl(const T* Ja, const T* hn, T* BL0) {}
};
