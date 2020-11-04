#include <glog/logging.h>
#include <fem/mat/elastic.h>
#include <fem/element/brick.h>
#include "common_brick.h"
#include <iostream>

extern "C" {
  void e_c3d_(
    double* co,

    int* kon,
    char* lakonl,
    double* p1,
    double* p2,
    double* omx,
    double* bodyfx,
    int* nbody,
    double* s,
    double* sm,
    double* ff,

    int* nelem,
    int* nmethod,
    double* elcon,
    int* nelcon,
    double* rhcon,
    int* nrhcon,
    double* alcon,
    int* nalcon,
    double* alzero,
    int* ielmat,

    int* ielorien,
    int* norien,
    double* orab,
    int* ntmat_,
    double* t0,
    double* t1,
    int* ithermal,
    double* vold,
    int* iperturb,
    int* nelemload,

    char* sideload,
    double* xload,
    int* nload,
    int* idist,
    double* sti,
    double* stx,
    int* iexpl,
    double* plicon,
    int* nplicon,
    double* plkcon,

    int* nplkcon,
    double* xstiff,
    int* npmat_,
    double* dtime,
    char* matname,
    int* mi,
    int* ncmat_,
    int* mass,
    int* stiffness,
    int* buckling,

    int* rhsi,
    int* intscheme,
    double* ttime,
    double* time,
    int* istep,
    int* iinc,
    int* coriolis,
    double* xloadold,
    double* reltime,
    int* ipompc,

    int* nodempc,
    double* coefmpc,
    int* nmpc,
    int* ikmpc,
    int* ilmpc,
    double* veold,
    double* springarea,
    int* nstate_,
    double* xstateini,
    double* xstate,

    int* ne0,
    int* ipkon,
    double* thicke,
    int* integerglob,
    double* doubleglob,
    char* tieset,
    int* istartset,
    int* iendset,
    int* ialset,
    int* ntie,

    int* nasym,
    double* pslavsurf,
    double* pmastsurf,
    int* mortar,
    double* clearini,
    int* ielprop,
    double* prop,
    int* kscale,
    double* smscalel,
    int* mscalmethod
  );
}

template <
  template <typename> class MatType,
  template <typename> class BrickIPropType,
  template <typename> class BrickTLFormType
>
void test_brick_form(
  bool layout=true, double tol=1e-6, unsigned int form=0) {
  double E = 2e5;
  double nu = 0.3;

  // execute
  double p[2] = {E, nu};
  MatType<double> m(p);
  BrickIPropType<double> bi;
  auto Dim = BrickIPropType<double>::get_ndim();
  auto N = BrickIPropType<double>::get_num_nodes();
  // init X0
  auto nEntryX0 = Dim * N;
  auto X0 = (double*)malloc(nEntryX0*sizeof(double));
  // for debug purpose
  // print_mat<double>(X0, Dim, N);
  init_x<double>(X0, Dim, N);
  // init hbuf
  auto hbuf = bi.get_hbuf();
  // init weights
  auto weights = bi.get_weights();
  // init Ut
  auto nEntryUt = Dim * N;
  auto Ut = (double*)malloc(nEntryUt*sizeof(double));
  init_zero<double>(Ut, nEntryUt);
  // init C0
  auto C0 = m.get_C();
  // init S0t
  auto nEntryS0t = Dim * Dim;
  auto S0t = (double*)malloc(nEntryS0t*sizeof(double));
  init_rand<double>(S0t, nEntryS0t);
  // init J0
  auto nEntryJ0 = Dim * Dim;
  auto J0 = (double*)malloc(nEntryJ0*sizeof(double));
  // init invJ0
  auto nEntryInvJ0 = Dim * Dim;
  auto invJ0 = (double*)malloc(nEntryInvJ0*sizeof(double));
  // init BdilBar
  auto nRowBdilBar = 3*Dim - 3;
  auto nColBdilBar = Dim * N;
  auto nEntryBdilBar = nRowBdilBar * nColBdilBar;
  auto BdilBar = (double*)malloc(nEntryBdilBar*sizeof(double));
  // init tmpB
  auto nRowTmpB = 3*Dim - 3;
  auto nColTmpB = Dim * N;
  auto nEntryTmpB = nRowTmpB * nColTmpB;
  auto tmpB = (double*)malloc(nEntryTmpB*sizeof(double));
  // init h0
  auto nEntryH0 = N * Dim;
  auto h0 = (double*)malloc(nEntryH0*sizeof(double));
  // init u0t
  auto nEntryU0t = Dim * Dim;
  auto u0t = (double*)malloc(nEntryU0t*sizeof(double));
  // init B0tL
  auto nRowB0tL = 3*Dim - 3;
  auto nColB0tL = Dim * N;
  auto nEntryB0tL = nRowB0tL * nColB0tL;
  auto B0tL = (double*)malloc(nEntryB0tL*sizeof(double));
  // init buf
  auto nEntryBuf = Dim * Dim;
  auto buf = (double*)malloc(nEntryBuf*sizeof(double));
  // init tmpK
  auto nRowTmpK = Dim * N;
  auto nEntryTmpK = nRowTmpK * nRowTmpK;
  auto tmpK = (double*)malloc(nEntryTmpK*sizeof(double));
  // init B0NL
  auto nRowB0NL = Dim * Dim;
  auto nColB0NL = Dim * N;
  auto nEntryB0NL = nRowB0NL * nColB0NL;
  auto B0NL = (double*)malloc(nEntryB0NL*sizeof(double));
  // init tile
  auto nRowTile = Dim * Dim;
  auto nEntryTile = nRowTile * nRowTile;
  auto tile = (double*)malloc(nEntryTile*sizeof(double));
  // init Ke
  auto nRowKe = Dim * N;
  auto nEntryKe = nRowKe * nRowKe;
  auto Ke = (double*)malloc(nEntryKe*sizeof(double));
  // run
  int ret;
  if (form == 0) {
    ret = BrickTLFormType<double>::form_elem_stiff(
      X0, hbuf, weights, Ut, C0, S0t, J0, invJ0, h0, u0t, B0tL, buf, tmpK, B0NL, tile, Ke);
  } else if (form == 1) {
    ret = BrickTLFormType<double>::form_linear_elem_stiff(
      X0, hbuf, weights, C0, J0, invJ0, h0, B0tL, buf, tmpK, Ke);
  } else {
    ret = BrickTLFormType<double>::form_linear_elem_stiff_Bbar(
      X0, hbuf, weights, C0,
      J0, invJ0, 
      BdilBar, tmpB, h0, B0tL, buf, tmpK, Ke);
  }
  if (!check_sym<double>(Ke, nRowKe, tol)) {
    ret = -1;
  } else {
    LOG(INFO) << "symmetrical Ke acquired";
  } 
  if (ret == 0) {
    if (layout) {
      LOG(INFO) << "matrix Ke layout"; 
      print_mat<double>(Ke, nRowKe, nRowKe);
    }
  } 
  std::cout << Ke[59*60+59] << std::endl;
  
  // validate
  double _co[3*20];
  init_x<double>(_co, 3, 20);
  double co[60]; /* input */
  transpose<double>(_co, 3, 20, co);

  int kon[20] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20}; /* input */
  char lakonl[8] = "C3D20R"; /* input */
  double* p1 = NULL;
  double* p2 = NULL;
  double omx = 0;
  double* bodyfx = NULL;
  int nbody = 0;
  double s[3600]; /* output */
  double* sm = NULL;
  double* ff = NULL;

  int nelem = 1; /* input */
  int nmethod = 0;
  double* elcon = NULL; /* input, not runtime */
  int nelcon[2] = {2,1}; /* input */
  double rhcon[2] = {0,0}; /* input */
  int* nrhcon = NULL;
  double* alcon = NULL;
  int* nalcon = NULL;
  double* alzero = NULL;
  int ielmat[1] = {1}; /* input */

  int* ielorien = NULL;
  int norien = 0;
  double* orab = NULL;
  int ntmat_ = 1; /* input */
  double* t0 = NULL;
  double* t1 = NULL;
  int ithermal[2] = {0,0}; /* input */
  double* vold = NULL;
  int iperturb[2] = {0,0}; /* input */
  int* nelemload = NULL;

  char* sideload = NULL;
  double* xload = NULL;
  int nload = 0;
  int idist = 0;
  double* sti = NULL;
  double* stx = NULL;
  int iexpl = 0;
  double* plicon = NULL;
  int* nplicon = NULL;
  double* plkcon = NULL;

  int* nplkcon = NULL;
  double _xstiff[216] = {
    E, E, E, E, E, E, E, E,
    nu, nu, nu, nu, nu, nu, nu, nu,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0
  };
  double xstiff[216]; /* input */
  transpose<double>(_xstiff, 27, 8, xstiff);
  int npmat_ = 0;
  double dtime = 0;
  char matname[80] = "STEEL"; /* input */
  int mi[3] = {8,3,1}; /* input */
  int ncmat_ = 2; /* input */
  int mass = 0; /* input */
  int stiffness = 1; /* input */
  int buckling = 0; /* input */

  int rhsi = 0;
  int intscheme = 0; /* input */
  double ttime = 0;
  double time = 0;
  int istep = 0;
  int iinc = 0;
  int coriolis = 0; /* input */
  double* xloadold = NULL;
  double reltime = 0;
  int* ipompc = NULL;

  int* nodempc = NULL;
  double* coefmpc = NULL;
  int nmpc = 0;
  int* ikmpc = NULL;
  int* ilmpc = NULL;
  double* veold = NULL;
  double* springarea = NULL;
  int nstate_ = 0;
  double* xstateini = NULL;
  double* xstate = NULL;

  int ne0 = 0;
  int ipkon[1] = {0}; /* input */
  double* thicke = NULL;
  int* integerglob = NULL;
  double* doubleglob = NULL;
  char* tieset = NULL;
  int* istartset = NULL;
  int* iendset = NULL;
  int* ialset = NULL;
  int ntie = 0;

  int nasym = 0;
  double* pslavsurf = NULL;
  double* pmastsurf = NULL;
  int mortar = 0;
  double* clearini = NULL;
  int* ielprop = NULL;
  double* prop = NULL;
  int kscale = 0;
  double smscalel = 0;
  int mscalmethod = 0;

  // validate
  e_c3d_(
    co,

    kon,
    lakonl,
    p1,
    p2,
    &omx,
    bodyfx,
    &nbody,
    s,
    sm,
    ff,

    &nelem,
    &nmethod,
    elcon,
    nelcon,
    rhcon,
    nrhcon,
    alcon,
    nalcon,
    alzero,
    ielmat,

    ielorien,
    &norien,
    orab,
    &ntmat_,
    t0,
    t1,
    ithermal,
    vold,
    iperturb,
    nelemload,

    sideload,
    xload,
    &nload,
    &idist,
    sti,
    stx,
    &iexpl,
    plicon,
    nplicon,
    plkcon,

    nplkcon,
    xstiff,
    &npmat_,
    &dtime,
    matname,
    mi,
    &ncmat_,
    &mass,
    &stiffness,
    &buckling,

    &rhsi,
    &intscheme,
    &ttime,
    &time,
    &istep,
    &iinc,
    &coriolis,
    xloadold,
    &reltime,
    ipompc,

    nodempc,
    coefmpc,
    &nmpc,
    ikmpc,
    ilmpc,
    veold,
    springarea,
    &nstate_,
    xstateini,
    xstate,

    &ne0,
    ipkon,
    thicke,
    integerglob,
    doubleglob,
    tieset,
    istartset,
    iendset,
    ialset,
    &ntie,

    &nasym,
    pslavsurf,
    pmastsurf,
    &mortar,
    clearini,
    ielprop,
    prop,
    &kscale,
    &smscalel,
    &mscalmethod
  );
  full_sym<double>(s, 60);
  std::cout << s[59*60+59] << std::endl;
  bool flag  = validate<double>(Ke, s, 3600, tol);
  // free
  free(X0);
  free(Ut);
  free(S0t);
  free(J0);
  free(invJ0);
  free(BdilBar);
  free(tmpB);
  free(h0);
  free(u0t);
  free(B0tL);
  free(buf);
  free(tmpK);
  free(B0NL);
  free(tile);
  free(Ke);
  if (flag) {
    LOG(INFO) << "test_brick_tl_form succeed, T: " << typeid(double).name()
      << ", BrickIPropType: " << typeid(BrickIPropType<double>).name()
      << ", BrickTLFormType: " << typeid(BrickTLFormType<double>).name();
  } else {
    LOG(INFO) << "test_brick_tl_form fail, T: " << typeid(double).name()
      << ", BrickIPropType: " << typeid(BrickIPropType<double>).name()
      << ", BrickTLFormType: " << typeid(BrickTLFormType<double>).name();
  }
}

template <
  template <typename> class MatType,
  template <typename> class BrickIPropType,
  template <typename> class BrickTLFormType
>
void test_brick_form_gen(
  bool layout=true, double tol=1e-6, unsigned int form=0) {
  double E = 2e5;
  double nu = 0.3;

  // execute
  double p[2] = {E, nu};
  MatType<double> m(p);
  BrickIPropType<double> bi;
  auto Dim = BrickIPropType<double>::get_ndim();
  auto N = BrickIPropType<double>::get_num_nodes();
  std::cout << "Dim: " << Dim << ", N: " << N << std::endl;
  // init X0
  auto nEntryX0 = Dim * N;
  auto X0 = (double*)malloc(nEntryX0*sizeof(double));
  // for debug purpose
  // print_mat<double>(X0, Dim, N);
  init_x<double>(X0, Dim, N);
  // init hbuf
  auto hbuf = bi.get_hbuf();
  // init weights
  auto weights = bi.get_weights();
  // init Ut
  auto nEntryUt = Dim * N;
  auto Ut = (double*)malloc(nEntryUt*sizeof(double));
  init_zero<double>(Ut, nEntryUt);
  // init C0
  auto C0 = m.get_C();
  // init S0t
  auto nEntryS0t = Dim * Dim;
  auto S0t = (double*)malloc(nEntryS0t*sizeof(double));
  init_rand<double>(S0t, nEntryS0t);
  // init J0
  auto nEntryJ0 = Dim * Dim;
  auto J0 = (double*)malloc(nEntryJ0*sizeof(double));
  // init invJ0
  auto nEntryInvJ0 = Dim * Dim;
  auto invJ0 = (double*)malloc(nEntryInvJ0*sizeof(double));
  // init BdilBar
  auto nRowBdilBar = 3*Dim - 3;
  auto nColBdilBar = Dim * N;
  auto nEntryBdilBar = nRowBdilBar * nColBdilBar;
  auto BdilBar = (double*)malloc(nEntryBdilBar*sizeof(double));
  // init tmpB
  auto nRowTmpB = 3*Dim - 3;
  auto nColTmpB = Dim * N;
  auto nEntryTmpB = nRowTmpB * nColTmpB;
  auto tmpB = (double*)malloc(nEntryTmpB*sizeof(double));
  // init h0
  auto nEntryH0 = N * Dim;
  auto h0 = (double*)malloc(nEntryH0*sizeof(double));
  // init u0t
  auto nEntryU0t = Dim * Dim;
  auto u0t = (double*)malloc(nEntryU0t*sizeof(double));
  // init B0tL
  auto nRowB0tL = 3*Dim - 3;
  auto nColB0tL = Dim * N;
  auto nEntryB0tL = nRowB0tL * nColB0tL;
  auto B0tL = (double*)malloc(nEntryB0tL*sizeof(double));
  // init buf
  auto nEntryBuf = Dim * Dim;
  auto buf = (double*)malloc(nEntryBuf*sizeof(double));
  // init tmpK
  auto nRowTmpK = Dim * N;
  auto nEntryTmpK = nRowTmpK * nRowTmpK;
  auto tmpK = (double*)malloc(nEntryTmpK*sizeof(double));
  // init B0NL
  auto nRowB0NL = Dim * Dim;
  auto nColB0NL = Dim * N;
  auto nEntryB0NL = nRowB0NL * nColB0NL;
  auto B0NL = (double*)malloc(nEntryB0NL*sizeof(double));
  // init tile
  auto nRowTile = Dim * Dim;
  auto nEntryTile = nRowTile * nRowTile;
  auto tile = (double*)malloc(nEntryTile*sizeof(double));
  // init Ke
  auto nRowKe = Dim * N;
  auto nEntryKe = nRowKe * nRowKe;
  auto Ke = (double*)malloc(nEntryKe*sizeof(double));
  // run
  int ret;
  if (form == 0) {
    ret = BrickTLFormType<double>::form_elem_stiff(
      X0, hbuf, weights, Ut, C0, S0t, J0, invJ0, h0, u0t, B0tL, buf, tmpK, B0NL, tile, Ke);
  } else if (form == 1) {
    ret = BrickTLFormType<double>::form_linear_elem_stiff(
      X0, hbuf, weights, C0, J0, invJ0, h0, B0tL, buf, tmpK, Ke);
  } else {
    ret = BrickTLFormType<double>::form_linear_elem_stiff_Bbar(
      X0, hbuf, weights, C0,
      J0, invJ0, 
      BdilBar, tmpB, h0, B0tL, buf, tmpK, Ke);
  }
  if (!check_sym<double>(Ke, nRowKe, tol)) {
    ret = -1;
  } else {
    LOG(INFO) << "symmetrical Ke acquired";
  } 
  if (ret == 0) {
    if (layout) {
      LOG(INFO) << "matrix Ke layout"; 
      print_mat<double>(Ke, nRowKe, nRowKe);
    }
  } 
  std::cout << Ke[23*60+23] << std::endl;
  
  // validate
  double _co[3*8];
  init_x<double>(_co, 3, 8);
  double co[24]; /* input */
  transpose<double>(_co, 3, 8, co);

  int kon[20] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20}; /* input */
  char lakonl[8] = "C3D8I"; /* input */
  double* p1 = NULL;
  double* p2 = NULL;
  double omx = 0;
  double* bodyfx = NULL;
  int nbody = 0;
  double s[3600]; /* output */
  double* sm = NULL;
  double* ff = NULL;

  int nelem = 1; /* input */
  int nmethod = 0;
  double* elcon = NULL; /* input, not runtime */
  int nelcon[2] = {2,1}; /* input */
  double rhcon[2] = {0,0}; /* input */
  int* nrhcon = NULL;
  double* alcon = NULL;
  int* nalcon = NULL;
  double* alzero = NULL;
  int ielmat[1] = {1}; /* input */

  int* ielorien = NULL;
  int norien = 0;
  double* orab = NULL;
  int ntmat_ = 1; /* input */
  double* t0 = NULL;
  double* t1 = NULL;
  int ithermal[2] = {0,0}; /* input */
  double* vold = NULL;
  int iperturb[2] = {0,0}; /* input */
  int* nelemload = NULL;

  char* sideload = NULL;
  double* xload = NULL;
  int nload = 0;
  int idist = 0;
  double* sti = NULL;
  double* stx = NULL;
  int iexpl = 0;
  double* plicon = NULL;
  int* nplicon = NULL;
  double* plkcon = NULL;

  int* nplkcon = NULL;
  double _xstiff[216] = {
    E, E, E, E, E, E, E, E,
    nu, nu, nu, nu, nu, nu, nu, nu,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0
  };
  double xstiff[216]; /* input */
  transpose<double>(_xstiff, 27, 8, xstiff);
  int npmat_ = 0;
  double dtime = 0;
  char matname[80] = "STEEL"; /* input */
  int mi[3] = {8,3,1}; /* input */
  int ncmat_ = 2; /* input */
  int mass = 0; /* input */
  int stiffness = 1; /* input */
  int buckling = 0; /* input */

  int rhsi = 0;
  int intscheme = 0; /* input */
  double ttime = 0;
  double time = 0;
  int istep = 0;
  int iinc = 0;
  int coriolis = 0; /* input */
  double* xloadold = NULL;
  double reltime = 0;
  int* ipompc = NULL;

  int* nodempc = NULL;
  double* coefmpc = NULL;
  int nmpc = 0;
  int* ikmpc = NULL;
  int* ilmpc = NULL;
  double* veold = NULL;
  double* springarea = NULL;
  int nstate_ = 0;
  double* xstateini = NULL;
  double* xstate = NULL;

  int ne0 = 0;
  int ipkon[1] = {0}; /* input */
  double* thicke = NULL;
  int* integerglob = NULL;
  double* doubleglob = NULL;
  char* tieset = NULL;
  int* istartset = NULL;
  int* iendset = NULL;
  int* ialset = NULL;
  int ntie = 0;

  int nasym = 0;
  double* pslavsurf = NULL;
  double* pmastsurf = NULL;
  int mortar = 0;
  double* clearini = NULL;
  int* ielprop = NULL;
  double* prop = NULL;
  int kscale = 0;
  double smscalel = 0;
  int mscalmethod = 0;

  // validate
  e_c3d_(
    co,

    kon,
    lakonl,
    p1,
    p2,
    &omx,
    bodyfx,
    &nbody,
    s,
    sm,
    ff,

    &nelem,
    &nmethod,
    elcon,
    nelcon,
    rhcon,
    nrhcon,
    alcon,
    nalcon,
    alzero,
    ielmat,

    ielorien,
    &norien,
    orab,
    &ntmat_,
    t0,
    t1,
    ithermal,
    vold,
    iperturb,
    nelemload,

    sideload,
    xload,
    &nload,
    &idist,
    sti,
    stx,
    &iexpl,
    plicon,
    nplicon,
    plkcon,

    nplkcon,
    xstiff,
    &npmat_,
    &dtime,
    matname,
    mi,
    &ncmat_,
    &mass,
    &stiffness,
    &buckling,

    &rhsi,
    &intscheme,
    &ttime,
    &time,
    &istep,
    &iinc,
    &coriolis,
    xloadold,
    &reltime,
    ipompc,

    nodempc,
    coefmpc,
    &nmpc,
    ikmpc,
    ilmpc,
    veold,
    springarea,
    &nstate_,
    xstateini,
    xstate,

    &ne0,
    ipkon,
    thicke,
    integerglob,
    doubleglob,
    tieset,
    istartset,
    iendset,
    ialset,
    &ntie,

    &nasym,
    pslavsurf,
    pmastsurf,
    &mortar,
    clearini,
    ielprop,
    prop,
    &kscale,
    &smscalel,
    &mscalmethod
  );
  full_sym<double>(s, 60);
  std::cout << "1: " << s[23*60+23] << std::endl;
  std::cout << "2: " << s[2*60+13] << std::endl;
  std::cout << "3: " << s[9*60+12] << std::endl;
  std::cout << "4: " << s[2*60+3] << std::endl;
  std::cout << "5: " << s[22*60+18] << std::endl;
  bool flag  = validate<double>(Ke, s, 3600, tol);
  // free
  free(X0);
  free(Ut);
  free(S0t);
  free(J0);
  free(invJ0);
  free(BdilBar);
  free(tmpB);
  free(h0);
  free(u0t);
  free(B0tL);
  free(buf);
  free(tmpK);
  free(B0NL);
  free(tile);
  free(Ke);
  if (flag) {
    LOG(INFO) << "test_brick_tl_form succeed, T: " << typeid(double).name()
      << ", BrickIPropType: " << typeid(BrickIPropType<double>).name()
      << ", BrickTLFormType: " << typeid(BrickTLFormType<double>).name();
  } else {
    LOG(INFO) << "test_brick_tl_form fail, T: " << typeid(double).name()
      << ", BrickIPropType: " << typeid(BrickIPropType<double>).name()
      << ", BrickTLFormType: " << typeid(BrickTLFormType<double>).name();
  }
}

int main(int argc, char* argv[]) {
  // log init
  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr = 1;
  // tests
  // test_brick_form<
    // fem::Ela3D,fem::C3D20RIProp,fem::C3D20RTLForm>(false, 1e-6, 1);
  test_brick_form_gen<
    fem::Ela3D,fem::C3D8IProp,fem::C3D8TLForm>(false, 1e-6, 1);
  return 0;
}
