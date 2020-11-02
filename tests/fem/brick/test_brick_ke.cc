#include <glog/logging.h>
#include "common_brick.h"

extern "C" {
  void e_c3d_(
    double* co,

    int* kon,
    char* lakonl,
    double* p1,
    double* p2,
    double omx,
    double* bodyfx,
    int nbody,
    double* s,
    double* sm,
    double* ff,

    int nelem,
    int nmethod,
    double* elcon,
    double* nelcon,
    double* rhcon,
    int* nrhcon,
    double* alcon,
    int* nalcon,
    double* alzero,
    int* ielmat,

    int* ielorien,
    int norien,
    double* orab,
    int ntmat_,
    double* t0,
    double* t1,
    int* ithermal,
    double* vold,
    int* iperturb,
    int* nelemload,

    char* sideload,
    double* xload,
    int nload,
    int idist,
    double* sti,
    double* stx,
    int iexpl,
    double* plicon,
    int* nplicon,
    double* plkcon,

    int* nplkcon,
    double* xstiff,
    int npmat_,
    double dtime,
    char* matname,
    int* mi,
    int ncmat_,
    int mass,
    int stiffness,
    int buckling,

    int rhsi,
    int intscheme,
    double ttime,
    double time,
    int istep,
    int iinc,
    int coriolis,
    double* xloadold,
    double reltime,
    int* ipompc,

    int* nodempc,
    double* coefmpc,
    int nmpc,
    int* ikmpc,
    int* ilmpc,
    double* veold,
    double* springarea,
    int nstate_,
    double* xstateini,
    double* xstate,

    int ne0,
    int* ipkon,
    double* thicke,
    int* integerglob,
    double* doubleglob,
    char* tieset,
    int* istartset,
    int* iendset,
    int* ialset,
    int ntie,

    int nasym,
    double* pslavsurf,
    double* pmastsurf,
    int mortar,
    double* clearini,
    int* ielprop,
    double* prop,
    int kscale,
    double smscalel,
    int mscalmethod
  );
}

template <unsigned int Dim, unsigned int N>
void test() {
    double co[Dim*N];
    init_x<double>(co, Dim, N);

    // int* kon;
    // char lakonl[8] = "C3D20R";
    double* p1 = NULL;
    double* p2 = NULL;
    double omx = 0;
    double* bodyfx = NULL;
    int nbody = 0;
    // double s[3600];
    double* sm = NULL;
    double* ff = NULL;

    // int nelem;
    int nmethod = 0;
    // double* elcon;
    // double* nelcon;
    double* rhcon = NULL;
    int* nrhcon = NULL;
    double* alcon = NULL;
    int* nalcon = NULL;
    double* alzero = NULL;
    // int* ielmat;

    int* ielorien = NULL;
    int norien = 0;
    double* orab = NULL;
    // int ntmat_;
    double* t0 = NULL;
    double* t1 = NULL;
    // int* ithermal;
    double* vold = NULL;
    // int* iperturb;
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
    // double* xstiff;
    int npmat_ = 0;
    double dtime = 0;
    char* matname = NULL;
    // int* mi;
    // int ncmat_ = 0;
    // int mass;
    // int stiffness;
    // int buckling;

    int rhsi = 0;
    // int intscheme;
    double ttime = 0;
    double time = 0;
    int istep = 0;
    int iinc = 0;
    // int coriolis;
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
    // int* ipkon;
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
}

int main(int argc, char* argv[]) {
  // log init
  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr = 1;
  // tests
  test<3,20>();
  return 0;
}
