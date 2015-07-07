#ifndef __Cuts_h_
#define __Cuts_h_
#include "Charged_Particle.h"
#include "Hist_dir.h"
#include "user_NEW.h"
#include "Cuts.h"


class Cuts{
public:
    /* Cuts(superCmpEvent* sevt, superBurst *sburst); */
    Cuts();
    ~Cuts();
    double DCH_e1e2;
    double DCH_mue1;
    double DCH_mue2;

    double Hod_e1e2;
    double Hod_mue1;
    double Hod_mue2;
    double Mu_P;
    double E1_P;
    double E2_P;
    double muee_P;
    double mee;
    double muee_Pt;
    double zvtx_e1e2;
    double zvtx_mue2;
    double zvtx_mue1;
    double yvtx_e1e2;
    double yvtx_mue2;
    double yvtx_mue1;
    double xvtx_e1e2;
    double xvtx_mue2;
    double xvtx_mue1;
};
#endif
