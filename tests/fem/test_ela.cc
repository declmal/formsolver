#include <glog/logging.h>
#include <fem/mat/elastic.h>
#include <common/common.h>

template <typename T>
void test_ela_3d(bool layout=true) {
  T p[2] = {(T)1e8, (T)0.3};
  unsigned int nEntryP = (unsigned int) (sizeof(p) / sizeof(T));
  if (nEntryP != fem::Ela3D<T>::get_num_params()) {
    LOG(FATAL) << "invalid params list";
  }
  fem::Ela3D<T> m(p);
  auto NR = fem::Ela3D<T>::get_nrows();
  auto C = m.get_C();
  if (layout) {
    LOG(INFO) << "matrix C layout";
    print_mat<T>(C, NR, NR);
  }
  LOG(INFO) << "test_ela_3d passed, T: " << typeid(T).name();
}

int main(int argc, char* argv[]) {
  // log init
  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr = 1;
  bool layout = false;
  // double precison tests
  test_ela_3d<double>();
  LOG(INFO) << "double precision test passed";
  // single precision tests
  test_ela_3d<float>();
  LOG(INFO) << "single precision test passed";
  return 0;
}
