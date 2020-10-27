#include <glog/logging.h>
#include <fem/element/brick.h>
#include <fem/formulator/total_lagrangian.h>
#include "../common/common.h"

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
  form.form_elem_stiff();
}

int main(int argc, char* argv[]) {
  // log init
  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr = 1;
  // double precison tests
  test_form_brick<double, fem::C3D8, fem::TL3D>();
  // test_form_brick<double, fem::C3D8R, fem::TL3D>();
  test_form_brick<double, fem::C3D20, fem::TL3D>();
  test_form_brick<double, fem::C3D20R, fem::TL3D>();
  LOG(INFO) << "double precision test passed";
  // single precision tests
  test_form_brick<float, fem::C3D8, fem::TL3D>();
  // test_form_brick<float, fem::C3D8R, fem::TL3D>();
  test_form_brick<float, fem::C3D20, fem::TL3D>();
  test_form_brick<float, fem::C3D20R, fem::TL3D>();
  LOG(INFO) << "single precision test passed";
  return 0;
}
