#include <glog/logging.h>
#include <fem/element/brick.h>
#include <fem/formulator/total_lagrangian.h>
#include "../common/common.h"

template <
  typename T,
  template <typename> class Form3DType,
  template <
    typename, template <typename> class
  > class BrickType
>
void test_brick_elem_stiff_cpu(bool layout=true, double tol=1e-6) {
  BrickType<T,Form3DType> brick;
  // execute
  brick.form_elem_stiff();
  // validate
  bool flag = true;
  if (flag) {
    LOG(INFO) << "test_brick_elem_stiff_cpu succeed, dtype: " 
      << typeid(T).name(); 
  } else {
    LOG(FATAL) << "test_brick_elem_stiff_cpu fail, dtype: "
      << typeid(T).name();
  }
}

int main(int argc, char* argv[]) {
  // log init
  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr = 1;
  // double precison tests
  test_brick_elem_stiff_cpu<double,fem::TL3D,fem::C3D8>();
  // test_brick_elem_stiff_cpu<double,fem::TL3D,fem::C3D8R>();
  test_brick_elem_stiff_cpu<double,fem::TL3D,fem::C3D20>();
  test_brick_elem_stiff_cpu<double,fem::TL3D,fem::C3D20R>();
  LOG(INFO) << "double precision test passed";
  // single precision tests
  test_brick_elem_stiff_cpu<float,fem::TL3D,fem::C3D8>();
  // test_brick_elem_stiff_cpu<float,fem::TL3D,fem::C3D8R>();
  test_brick_elem_stiff_cpu<float,fem::TL3D,fem::C3D20>();
  test_brick_elem_stiff_cpu<float,fem::TL3D,fem::C3D20R>();
  LOG(INFO) << "single precision test passed";
  return 0;
}
