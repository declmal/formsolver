#include <glog/logging.h>
#include <fem/element/element.h>
#include <fem/formulator/total_lagrangian.h>
#include <common/common.h>

template <
  typename T,
  template <typename> class EType,
  template <
    typename, template<typename> class
  > class FormType
>
void test_form_brick(bool layout=true) {
  // init X0
  EType<T> elem(233);
  auto rowX0 = elem.get_ndim();
  auto colX0 = elem.get_num_nodes();
  auto nEntryX0 = rowX0 * colX0;
  auto X0 = (T*)malloc(nEntryX0*sizeof(T));
  init_rand<T>(X0, nEntryX0);
  // execute
  elem.init_coordinate(X0);
  FormType<T, EType> form(&elem);
  form.form_elem_stiff();
  if (layout) {
    LOG(INFO) << "matrix X0 layout";
    print_mat<T>(X0, rowX0, colX0);
  }
  // free
  free(X0);
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
