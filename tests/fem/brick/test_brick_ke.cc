#include <glog/logging.h>
#include <fem/mat/elastic.h>
#include "common_brick.h"

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

void test() {
  double E = 2e5;
  double nu = 0.3;
  // double p[2] = {E, nu};
  // fem::Ela3D<double> ela(p);
  
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
}

int main(int argc, char* argv[]) {
  // log init
  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr = 1;
  // tests
  test();
  return 0;
}
