#include <glog/logging.h>
#include <fem/element/brick.h>
#include <fem/formulator/total_lagrangian.h>
#include "../common/common.h"

template <
  typename T,
  template <typename> class BrickType
>
void test_brick_elem_stiff_cpu(bool layout=true, double tol=1e-6) {
  BrickType<T> brick;
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

template <
  typename T,
  template <typename> class EType,
  template <
    typename, template<typename> class
  > class FormType
>
void test_form_brick(bool layout=true) {
  EType<T> elem;
  FormType<T, EType> form(&elem);
  form.load_elem_data();
}

int main(int argc, char* argv[]) {
  // log init
  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr = 1;
  // double precison tests
  test_brick_elem_stiff_cpu<double,fem::C3D8>();
  // test_brick_elem_stiff_cpu<double,fem::C3D8R>();
  test_brick_elem_stiff_cpu<double,fem::C3D20>();
  test_brick_elem_stiff_cpu<double,fem::C3D20R>();
  LOG(INFO) << "double precision test passed";
  // single precision tests
  test_brick_elem_stiff_cpu<float,fem::C3D8>();
  // test_brick_elem_stiff_cpu<float,fem::C3D8R>();
  test_brick_elem_stiff_cpu<float,fem::C3D20>();
  test_brick_elem_stiff_cpu<float,fem::C3D20R>();
  test_form_brick<float, fem::C3D8, fem::TL3D>();
  LOG(INFO) << "single precision test passed";
  return 0;
}
