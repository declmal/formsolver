#include <glog/logging.h>

extern "C" {
  void e_c3d_();
}

void test() {
}

int main(int argc, char* argv[]) {
  // log init
  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr = 1;
  // tests
  test();
  return 0;
}
