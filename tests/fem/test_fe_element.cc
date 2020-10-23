#include <glog/logging.h>
#include <fem/fe_element.h>
#include "../common.h"

void test_brick3d_shp_deri_dp_cpu(bool layout=true) {
  // init r
  unsigned int nEntryR = 3;
  auto r = (double*)malloc(nEntryR*sizeof(double));
  init_rand<double>(r, nEntryR);
  // init h
  unsigned int nEntryH = 8 * 3;
  auto h = (double*)malloc(nEntryH*sizeof(double));
  // execute
  fem::brick3d_shp_deri<double>(r, h);
  if (layout) {
    LOG(INFO) << "matrix r layout";
    print_mat_dp(r, 3, 1);
    LOG(INFO) << "matrix h layout";
    print_mat_dp(h, 8, 3);
  }
  // free
  free(r);
  free(h);
  LOG(INFO) << "test_brick3d_shp_deri_dp_cpu succeed";
}

int main(int argc, char* argv[]) {
  // log init
  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr = 1;
  // tests
  test_brick3d_shp_deri_dp_cpu(true);
  return 0;
}
