#ifndef __Cuts_pi0d_h_
#define __Cuts_pi0d_h_
#include "Charged_Particle.h"
#include "MC_Charged_Particle.h"
#include "Hist_dir_pi0d.h"
#include "user_NEW.h"



class Cuts{
public:
    /* Cuts(superCmpEvent* sevt, superBurst *sburst); */
    Cuts();
    ~Cuts();
    double DCH_e1e2;
    double DCH_pie1;
    double DCH_pie2;
    double DCH_Radius_pi;
    double DCH_Radius_el1;
    double DCH_Radius_el2;
    double DCH_Radius_pi1;
    double DCH_Radius_pi2;
    double DCH_Radius_pi3;

    double Lkr_x_el1 ;
    double Lkr_x_el2 ;
    double Lkr_x_pi  ;
    double Lkr_y_el1 ;
    double Lkr_y_el2 ;
    double Lkr_y_pi  ;
    double Lkr_cut_el1;
    double Lkr_cut_el2;
    double Lkr_cut_pi;

    double PIV_x_el1 ;
    double PIV_x_el2 ;
    double PIV_x_pi  ;
    double PIV_y_el1 ;
    double PIV_y_el2 ;
    double PIV_y_pi  ;



    double Hod_e1e2;
    double cluster_e1e2;
    double cluster_pie1;
    double cluster_pie2;
    double Hod_pie1;
    double Hod_pie2;
    double Pi_P;
    double E1_P;
    double E2_P;
    double mee;
    double piee_P;
    double piee_M;
    double piee_Pt;
    double zvtx_pie1_pie2;
    double zvtx_pie1_e1e2;
    double zvtx_pie2_e1e2;
    double yvtx_pie1_pie2;
    double yvtx_pie1_e1e2;
    double yvtx_pie2_e1e2;
    double xvtx_pie1_pie2;
    double xvtx_pie1_e1e2;
    double xvtx_pie2_e1e2;

};
#endif
