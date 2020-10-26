#ifndef COMMON_COMMON_H_
#define COMMON_COMMON_H_

template <typename T>
void init_zero(T* const data, const unsigned int size) {
  for (unsigned int i = 0; i < size; ++i) {
    data[i] = (T)0;
  }
}

#endif // COMMON_COMMON_H_
