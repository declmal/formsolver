#include <glog/logging.h>
#include <common/gauss_legendre.h>
#include "common.h"
#include <stdio.h>

template <typename T, unsigned int N>
struct ValidateGauss {
  static bool impl(
    const T* const roots, const T* const weights, double tol=1e-6);
};

template <typename T>
struct ValidateGauss<T, 2> {
  static bool impl(
    const T* const roots, const T* const weights, double tol=1e-6) {
    const static T ref_r[2] = {
      -0.577350269189626,
      0.577350269189626,
    };
    const static T ref_w[2]= {
      1.000000000000000,
      1.000000000000000,
    };
    return
      validate<T>(roots, ref_r, 2, tol) && 
      validate<T>(weights, ref_w, 2, tol);
  }
};

template <typename T>
struct ValidateGauss<T, 3> {
  static bool impl(
    const T* const roots, const T* const weights, double tol=1e-6) {
    const static T ref_r[3] = {
      -0.774596669241483,
      0.000000000000000,
      0.774596669241483,
    };
    const static T ref_w[3]= {
      0.555555555555556,
      0.888888888888889,
      0.555555555555556,
    };
    return 
      validate<T>(roots, ref_r, 3, tol) && 
      validate<T>(weights, ref_w, 3, tol);
  }
};

template <typename T>
struct ValidateGauss<T, 4> {
  static bool impl(
    const T* const roots, const T* const weights, double tol=1e-6) {
    const static T ref_r[4] = {
      -0.861136311594053,
      -0.339981043584856,
      0.339981043584856,
      0.861136311594053,
    };
    const static T ref_w[4]= {
      0.347854845137454,
      0.652145154862546,
      0.652145154862546,
      0.347854845137454,
    };
    return 
      validate<T>(roots, ref_r, 4, tol) && 
      validate<T>(weights, ref_w, 4, tol);
  }
};

template <typename T>
struct ValidateGauss<T, 5> {
  static bool impl(
    const T* const roots, const T* const weights, double tol=1e-6) {
    const static T ref_r[5] = {
      -0.906179845938664,
      -0.538469310105683,
      0.000000000000000,
      0.538469310105683,
      0.906179845938664,
    };
    const static T ref_w[5]= {
      0.236926885056189,
      0.478628670499366,
      0.568888888888889,
      0.478628670499366,
      0.236926885056189,
    };
    return 
      validate<T>(roots, ref_r, 5, tol) && 
      validate<T>(weights, ref_w, 5, tol);
  }
};

template <typename T>
struct ValidateGauss<T, 6> {
  static bool impl(
    const T* const roots, const T* const weights, double tol=1e-6) {
    const static T ref_r[6] = {
      -0.932469514203152,
      -0.661209386466265,
      -0.238619186083197,
      0.238619186083197,
      0.661209386466265,
      0.932469514203152,
    };
    const static T ref_w[6]= {
      0.171324492379170,
      0.360761573048139,
      0.467913934572691,
      0.467913934572691,
      0.360761573048139,
      0.171324492379170,
    };
    return 
      validate<T>(roots, ref_r, 6, tol) && 
      validate<T>(weights, ref_w, 6, tol);
  }
};

template <typename T>
struct ValidateGauss<T, 7> {
  static bool impl(
    const T* const roots, const T* const weights, double tol=1e-6) {
    const static T ref_r[7] = {
      -0.949107912342759,
      -0.741531185599394,
      -0.405845151377397,
      0.000000000000000,
      0.405845151377397,
      0.741531185599394,
      0.949107912342759,
    };
    const static T ref_w[7]= {
      0.129484966168870,
      0.279705391489277,
      0.381830050505119,
      0.417959183673469,
      0.381830050505119,
      0.279705391489277,
      0.129484966168870,
    };
    return 
      validate<T>(roots, ref_r, 7, tol) && 
      validate<T>(weights, ref_w, 7, tol);
  }
};

template <typename T>
struct ValidateGauss<T, 8> {
  static bool impl(
    const T* const roots, const T* const weights, double tol=1e-6) {
    const static T ref_r[8] = {
      -0.960289856497536,
      -0.796666477413627,
      -0.525532409916329,
      -0.183434642495650,
      0.183434642495650,
      0.525532409916329,
      0.796666477413627,
      0.960289856497536,
    };
    const static T ref_w[8]= {
      0.101228536290376,
      0.222381034453374,
      0.313706645877887,
      0.362683783378362,
      0.362683783378362,
      0.313706645877887,
      0.222381034453374,
      0.101228536290376,
    };
    return 
      validate<T>(roots, ref_r, 8, tol) && 
      validate<T>(weights, ref_w, 8, tol);
  }
};

template <typename T>
struct ValidateGauss<T, 9> {
  static bool impl(
    const T* const roots, const T* const weights, double tol=1e-6) {
    const static T ref_r[9] = {
      -0.968160239507626,
      -0.836031107326636,
      -0.613371432700590,
      -0.324253423403809,
      0.000000000000000,
      0.324253423403809,
      0.613371432700590,
      0.836031107326636,
      0.968160239507626,
    };
    const static T ref_w[9]= {
      0.081274388361574,
      0.180648160694857,
      0.260610696402935,
      0.312347077040003,
      0.330239355001260,
      0.312347077040003,
      0.260610696402935,
      0.180648160694857,
      0.081274388361574,
    };
    return 
      validate<T>(roots, ref_r, 9, tol) && 
      validate<T>(weights, ref_w, 9, tol);
  }
};

template <typename T>
struct ValidateGauss<T, 10> {
  static bool impl(
    const T* const roots, const T* const weights, double tol=1e-6) {
    const static T ref_r[10] = {
      -0.973906528517172,
      -0.865063366688985,
      -0.679409568299024,
      -0.433395394129247,
      -0.148874338981631,
      0.148874338981631,
      0.433395394129247,
      0.679409568299024,
      0.865063366688985,
      0.973906528517172,
    };
    const static T ref_w[10]= {
      0.066671344308688,
      0.149451349150581,
      0.219086362515982,
      0.269266719309996,
      0.295524224714753,
      0.295524224714753,
      0.269266719309996,
      0.219086362515982,
      0.149451349150581,
      0.066671344308688,
    };
    return 
      validate<T>(roots, ref_r, 10, tol) && 
      validate<T>(weights, ref_w, 10, tol);
  }
};

template <typename T>
struct ValidateGauss<T, 11> {
  static bool impl(
    const T* const roots, const T* const weights, double tol=1e-6) {
    const static T ref_r[11] = {
      -0.978228658146057,
      -0.887062599768095,
      -0.730152005574049,
      -0.519096129110681,
      -0.269543155952345,
      0.000000000000000,
      0.269543155952345,
      0.519096129110681,
      0.730152005574049,
      0.887062599768095,
      0.978228658146057,
    };
    const static T ref_w[11]= {
      0.055668567116174,
      0.125580369464905,
      0.186290210927734,
      0.233193764591990,
      0.262804544510247,
      0.272925086777901,
      0.262804544510247,
      0.233193764591990,
      0.186290210927734,
      0.125580369464905,
      0.055668567116174,
    };
    return 
      validate<T>(roots, ref_r, 11, tol) && 
      validate<T>(weights, ref_w, 11, tol);
  }
};

template <typename T>
struct ValidateGauss<T, 12> {
  static bool impl(
    const T* const roots, const T* const weights, double tol=1e-6) {
    const static T ref_r[12] = {
      -0.981560634246719,
      -0.904117256370475,
      -0.769902674194305,
      -0.587317954286617,
      -0.367831498918180,
      -0.125333408511469,
      0.125333408511469,
      0.367831498918180,
      0.587317954286617,
      0.769902674194305,
      0.904117256370475,
      0.981560634246719,
    };
    const static T ref_w[12]= {
      0.047175336386512,
      0.106939325995318,
      0.160078328543346,
      0.203167426723066,
      0.233492536538355,
      0.249147045813403,
      0.249147045813403,
      0.233492536538355,
      0.203167426723066,
      0.160078328543346,
      0.106939325995318,
      0.047175336386512,
    };
    return 
      validate<T>(roots, ref_r, 12, tol) && 
      validate<T>(weights, ref_w, 12, tol);
  }
};

template <typename T>
struct ValidateGauss<T, 13> {
  static bool impl(
    const T* const roots, const T* const weights, double tol=1e-6) {
    const static T ref_r[13] = {
      -0.984183054718588,
      -0.917598399222978,
      -0.801578090733310,
      -0.642349339440340,
      -0.448492751036447,
      -0.230458315955135,
      0.000000000000000,
      0.230458315955135,
      0.448492751036447,
      0.642349339440340,
      0.801578090733310,
      0.917598399222978,
      0.984183054718588,
    };
    const static T ref_w[13]= {
      0.040484004765316,
      0.092121499837728,
      0.138873510219787,
      0.178145980761946,
      0.207816047536889,
      0.226283180262897,
      0.232551553230874,
      0.226283180262897,
      0.207816047536889,
      0.178145980761946,
      0.138873510219787,
      0.092121499837728,
      0.040484004765316,
    };
    return 
      validate<T>(roots, ref_r, 13, tol) && 
      validate<T>(weights, ref_w, 13, tol);
  }
};

template <typename T>
struct ValidateGauss<T, 14> {
  static bool impl(
    const T* const roots, const T* const weights, double tol=1e-6) {
    const static T ref_r[14] = {
      -0.986283808696812,
      -0.928434883663574,
      -0.827201315069765,
      -0.687292904811685,
      -0.515248636358154,
      -0.319112368927890,
      -0.108054948707344,
      0.108054948707344,
      0.319112368927890,
      0.515248636358154,
      0.687292904811685,
      0.827201315069765,
      0.928434883663574,
      0.986283808696812,
    };
    const static T ref_w[14]= {
      0.035119460331752,
      0.080158087159760,
      0.121518570687903,
      0.157203167158194,
      0.185538397477938,
      0.205198463721290,
      0.215263853463158,
      0.215263853463158,
      0.205198463721290,
      0.185538397477938,
      0.157203167158194,
      0.121518570687903,
      0.080158087159760,
      0.035119460331752,
    };
    return 
      validate<T>(roots, ref_r, 14, tol) && 
      validate<T>(weights, ref_w, 14, tol);
  }
};

template <typename T>
struct ValidateGauss<T, 15> {
  static bool impl(
    const T* const roots, const T* const weights, double tol=1e-6) {
    const static T ref_r[15] = {
      -0.987992518020485,
      -0.937273392400706,
      -0.848206583410427,
      -0.724417731360170,
      -0.570972172608539,
      -0.394151347077563,
      -0.201194093997435,
      0.000000000000000,
      0.201194093997435,
      0.394151347077563,
      0.570972172608539,
      0.724417731360170,
      0.848206583410427,
      0.937273392400706,
      0.987992518020485,
    };
    const static T ref_w[15]= {
      0.030753241996117,
      0.070366047488108,
      0.107159220467172,
      0.139570677926154,
      0.166269205816994,
      0.186161000015562,
      0.198431485327111,
      0.202578241925561,
      0.198431485327111,
      0.186161000015562,
      0.166269205816994,
      0.139570677926154,
      0.107159220467172,
      0.070366047488108,
      0.030753241996117,
    };
    return 
      validate<T>(roots, ref_r, 15, tol) && 
      validate<T>(weights, ref_w, 15, tol);
  }
};

template <typename T>
struct ValidateGauss<T, 16> {
  static bool impl(
    const T* const roots, const T* const weights, double tol=1e-6) {
    const static T ref_r[16] = {
      -0.989400934991650,
      -0.944575023073233,
      -0.865631202387832,
      -0.755404408355003,
      -0.617876244402644,
      -0.458016777657227,
      -0.281603550779259,
      -0.095012509837637,
      0.095012509837637,
      0.281603550779259,
      0.458016777657227,
      0.617876244402644,
      0.755404408355003,
      0.865631202387832,
      0.944575023073233,
      0.989400934991650,
    };
    const static T ref_w[16]= {
      0.027152459411754,
      0.062253523938648,
      0.095158511682493,
      0.124628971255534,
      0.149595988816577,
      0.169156519395003,
      0.182603415044924,
      0.189450610455069,
      0.189450610455069,
      0.182603415044924,
      0.169156519395003,
      0.149595988816577,
      0.124628971255534,
      0.095158511682493,
      0.062253523938648,
      0.027152459411754,
    };
    return 
      validate<T>(roots, ref_r, 16, tol) && 
      validate<T>(weights, ref_w, 16, tol);
  }
};

template <typename T, unsigned int N>
void test_gauss_1d(bool layout=true, double tol=1e-6) {
  GaussRoots1D<T, N> gr;
  GaussWeights1D<T, N> gw;
  if (layout) {
    LOG(INFO) << "matrix roots layout";
    print_mat<T>(gr.roots, N, 1);
    LOG(INFO) << "matrix weights layout";
    print_mat<T>(gw.weights, N, 1);
  }
  bool flag = ValidateGauss<T, N>::impl(gr.roots, gw.weights, tol);
  if (flag) {
    LOG(INFO) << "test_gauss_1d succeed, type: " 
      << typeid(T).name() << ", N: " << N;
  } else {
    LOG(FATAL) << "test_gauss_1d failed, type: " 
      << typeid(T).name() << ", N: " << N;
  }
}

template <typename T, unsigned int N0, unsigned int N1, unsigned int N2>
void test_gauss_3d(bool layout=true) {
  GaussRoots3D<T, N0, N1, N2> gr;
  GaussWeights3D<T, N0, N1, N2> gw;
  if (layout) {
    LOG(INFO) << "matrix roots layout";
    print_mat<T>(gr.roots, N0*N1*N2, 3);
    LOG(INFO) << "matrix weights layout";
    print_mat<T>(gw.weights, N0*N1*N2, 1);
  }
  LOG(INFO) << "test_gauss_3d succeed, type: " << typeid(T).name() 
    << ", N0: " << N0 << ", N1: " << N1 << ", N2: " << N2;
}

int main(int argc, char* argv[]) {
  // log init
  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr = 1;
  // double precision tests
  test_gauss_1d<double, 2>(false);
  test_gauss_1d<double, 3>(false);
  test_gauss_1d<double, 4>(false);
  test_gauss_1d<double, 5>(false);
  test_gauss_1d<double, 6>(false);
  test_gauss_1d<double, 7>(false);
  test_gauss_1d<double, 8>(false);
  test_gauss_1d<double, 9>(false);
  test_gauss_1d<double, 10>(false);
  test_gauss_1d<double, 11>(false);
  test_gauss_1d<double, 12>(false, 1e-3);
  test_gauss_1d<double, 13>(false);
  test_gauss_1d<double, 14>(false);
  test_gauss_1d<double, 15>(false);
  test_gauss_1d<double, 16>(false);
  test_gauss_3d<double, 4, 3, 2>();
  LOG(INFO) << "double precision test passed";
  // single precision tests
  test_gauss_1d<float, 2>(false);
  test_gauss_1d<float, 3>(false);
  test_gauss_1d<float, 4>(false);
  test_gauss_1d<float, 5>(false);
  test_gauss_1d<float, 6>(false);
  test_gauss_1d<float, 7>(false);
  test_gauss_1d<float, 8>(false);
  test_gauss_1d<float, 9>(false);
  test_gauss_1d<float, 10>(false);
  test_gauss_1d<float, 11>(false);
  test_gauss_1d<float, 12>(false, 1e-3);
  test_gauss_1d<float, 13>(false, 1e-5);
  test_gauss_1d<float, 14>(false, 1e-5);
  test_gauss_1d<float, 15>(false, 1e-5);
  test_gauss_1d<float, 16>(false);
  test_gauss_3d<float, 4, 3, 2>();
  LOG(INFO) << "double precision test passed";
  LOG(INFO) << "single precision test passed";
  return 0;
}
