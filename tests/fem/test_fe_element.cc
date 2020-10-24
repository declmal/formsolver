#include <glog/logging.h>
#include <fem/fe_element.h>
#include "../common/common.h"

template <typename T>
void test_brickn8_shp_deri_cpu(bool layout=true) {
  // init r
  unsigned int nEntryR = 3;
  auto r = (T*)malloc(nEntryR*sizeof(T));
  init_rand<T>(r, nEntryR);
  // init h
  unsigned int nEntryH = 8 * 3;
  auto h = (T*)malloc(nEntryH*sizeof(T));
  // execute
  fem::brickn8_shp_deri<T>(r, h);
  if (layout) {
    LOG(INFO) << "matrix r layout";
    print_mat<T>(r, 3, 1);
    LOG(INFO) << "matrix h layout";
    print_mat<T>(h, 8, 3);
  }
  // free
  free(r);
  free(h);
  LOG(INFO) << "test_brickn8_shp_deri_cpu succeed";
}

int main(int argc, char* argv[]) {
  // log init
  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr = 1;
  // double precison tests
  test_brickn8_shp_deri_cpu<double>(false);
  LOG(INFO) << "double precision test passed";
  // single precision tests
  test_brickn8_shp_deri_cpu<float>(false);
  LOG(INFO) << "single precision test passed";
  return 0;
}
