#include <glog/logging.h>
#include "common_brick.h"

extern "C" {
  void e_c3d_();
}

template <unsigned int Dim, unsigned int N>
void test() {
  double co[Dim*N];
}

int main(int argc, char* argv[]) {
  // log init
  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr = 1;
  // tests
  test<3,20>();
  return 0;
}
