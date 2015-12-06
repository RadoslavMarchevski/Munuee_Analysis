#include "cmpio.h"
//#include "user.h"
#include "user_NEW.h"
#include "reader.h"

#include <constants.h>

#include <TAxis.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TDirectory.h>

#include <iostream>
#include <vector>
#include <stdio.h>
#include <math.h>
//User defined libraries
#include "Charged_Particle.h"
#include "Neutral_Particle.h"
/* #include "Hist_dir.h" */
#include "Hist_dir_pi0d.h"
/* #include "Cuts.h" */
#include "Cuts_pi0d.h"
using namespace std;

vector<TVector3*> clcorrpos;
TVector3* clpos;
//float __lda3(superCmpEvent*, superCmpEvent*, int);

void cleanup();

void FillK2pi(superCmpEvent* sevt, Hist_dir* dir, TLorentzVector pion,TLorentzVector el1,TLorentzVector el2,TLorentzVector gamma);
void Filling(superCmpEvent* sevt, superBurst* sbur, Hist_dir* dir,int ipi,int iel1,int iel2, int igamma);
Cuts make_cuts(Hist_dir* dir1, Charged_Particle& muon, Charged_Particle& electron1, Charged_Particle& electron2,double* Vertex_mu_e1,double* Vertex_mu_e2,double* Vertex_e1_e2, std::string decay_type);


double ldaMC(float p, float lda);
double ldaC(superCmpEvent *sevt,superBurst *sbur,int i);
float sigmoidl(double x);

int user_superCmpEvent(superBurst *sbur,superCmpEvent *sevt) {
    /* WARNING: do not alter things before this line */
    /*---------- Add user C code here ----------*/
    //cout << "Start" << endl;
    static int nuserevt=0;
    nuserevt++;
    //Particle Indexes
    int imu=-1;
    int iel1=-1;
    int iel2=-1;
    int ipi=-1;
    int pi1=-1;
    int pi2=-1;
    int pi3=-1;
    int igamma=-1;


    bool LKrCalCorr=true;

    //--------------------------------------------------------//
    //Definition of Variables
    int    ntrack = sevt->Ntrack; // number of tracks
    int    nclust = sevt->Ncluster; // number of clusters
    int    nvtx   = sevt->Nvtx;   // number of vtx
    double DCHbz = Geom->DCH.bz;         // z before magnet
    double DCHz = Geom->DCH.z;         // z after magnet
    double LKrz=Geom->Lkr.z;
    double COmPaCt_Z_Vertex = sevt->vtx[0].z;
    double chi2_Z_Vertex = sevt->vtx[0].chi2;
    static double massKaonC = 0.493677;
    static double massPi0C = 0.1349766;
    const double Electron_EoverP_up = 1.05;
    const double Electron_EoverP_down = 0.95;
    const double Muon_EoverP_up = 0.1;
    int Kcharge=0;
    int ngoodtrack=0;
    int Nelectrons=0;
    int NK3pi_pions=0;
    int Nmuons=0;
    int Track_icl;
    int Track_imu;
    int Event_Type=-1;
    int K3pi_Event_Type=-1;
    //Kinematic variables
    //float lda3;
    float lda3_e1;
    float lda3_e2;
    float lda3_pi1;
    float lda3_pi2;
    float lda3_pi3;
    double Track_Momentum;
    double Track_Charge;
    double Track_Quality;
    double Track_DeadCell_Distance;
    double Track_Energy;
    double Track_EoverP;
    double PiTrack_EoverP;
    double Track_Lkr_status;
    double Temp;
    TLorentzVector Gamma_P;
    TLorentzVector Pi0d;
    TLorentzVector K2Pi0d;
    /// apply non linearity corrections
    /// Correct LKr non-linearity (effective for energies < 11GeV)
    if (LKrCalCorr && IS_DATA) {
        user_lkrcalcor_SC (sbur,sevt,1);
    }

    if(COmPaCt_Z_Vertex < -1800. || COmPaCt_Z_Vertex > 8000.){return 0;}

    if(nvtx != 1){return 0;}
    if(chi2_Z_Vertex > 20){return 0;}
    if(ntrack != 3 ) {return 0;}
    if(nclust < 3 ) {return 0;}

    for (int j=0; j<ntrack;j++){
        Kcharge += sevt->track[j].q;
    }

    //Looping over the number of tracks for each event
    for (int i=0; i<ntrack; i++) {

        Track_Momentum          = sevt->track[i].p;
        Track_Charge            = sevt->track[i].q;
        Track_Quality           = sevt->track[i].quality;
        Track_DeadCell_Distance = sevt->track[i].dDeadCell;
        Track_icl               = sevt->track[i].iClus;
        Track_imu               = sevt->track[i].iMuon;
        Track_Energy  = sevt->cluster[Track_icl].energy;
        Track_EoverP  = Track_Energy / Track_Momentum;

        //Quality controll

        if(Track_icl > -1)
            Track_Lkr_status     = sevt->cluster[Track_icl].status;

        if (Track_Momentum < 3.) {continue;}
        if (Track_imu != - 1.){continue;}
        if (Track_DeadCell_Distance < 2.){continue;}
        if (Track_Quality < 0.8){continue;}
        ngoodtrack++;

        //Particle identification for K to mu nu e e
        if(Kcharge*Track_Charge == -1. ){

            if( iel1 < 0){
                iel1 = i;

            }

        } else if ( iel2 < 0 || ipi < 0 ){

            if(ipi < 0 && Track_EoverP < 0.85){
                ipi = i;
                //cout << "iel2 = " << iel2 << " ipi =" << ipi << endl;
                //cout << "Kcharge = " << Kcharge << "Track_charge = " << Track_Charge << " ipi =" << ipi << endl;
                PiTrack_EoverP= sevt->cluster[sevt->track[ipi].iClus].energy/sevt->track[ipi].p;
                //cout << " Pion E/p =  " << PiTrack_EoverP << "Electron E/p =  " << Track_EoverP << endl;
                //cout << "NEXT _________+_+_+)_+)_+)_+)_+)_++__++" << endl;
            } else if(ipi > 0 && Track_EoverP > PiTrack_EoverP && Track_EoverP > 0.95){
                iel2 = i;
                //lda3_e2 = ldaC(sevt,sbur,i);
                //cout << " iel2 =  " << iel2 << "ipi =  " <<  ipi << endl;
                //cout << " Pion E/p =  " << PiTrack_EoverP << "Electron E/p =  " << Track_EoverP << endl;
            }else {
                if(Track_EoverP > 0.95){
                    iel2 = i;
                    lda3_e2 = ldaC(sevt,sbur,i);
                    //     ipi  = i;
                }

            }
        }

    }//end for i
    //K2piD selection
    Charged_Particle pion(sevt,sbur,211,ipi);
    Charged_Particle electron1(sevt,sbur,-11,iel1);
    Charged_Particle electron2(sevt,sbur,-11,iel2);

//Requiring only one high energetic cluster without associated track
    double ngamma = 0;
    for (int j = 0; j< nclust; j++){

        if( sevt->cluster[j].iTrack < 0 && sevt->cluster[j].energy > 2. ){
            igamma = j;
            ngamma++;
        }
    }

    if (ngamma != 1) return 0;


    //Checking for goodness of tracks
    if(ngoodtrack!= 3){return 0;}
    //-------------------------KPi0D SIGNAL SELECTION ----------------------------------------
    //Signal particle identification
    if(ipi == -1 || iel1 == -1 || iel2 == -1 || igamma ){return 0;}
    if(!electron1.cluster_exists || !electron2.cluster_exists){return 0;}

    Neutral_Particle gamma(sevt,sbur,igamma);


    K2Pi0d = pion.Momentum + electron1.Momentum + electron2.Momentum + gamma.Momentum;
    Pi0d   = electron1.Momentum + electron2.Momentum + gamma.Momentum;

    if (fabs(K2Pi0d.P() - 60.) > 5.){return 0;}
    if (K2Pi0d.Pt() > 0.01){return 0;}
    if (fabs(K2Pi0d.M() - massKaonC) > 0.008){return 0;}
    if (fabs(Pi0d.M() - massPi0C) > 0.008){return 0;}
    //--Vertex Reconstruction --
    //Vertex reconstruction for each pair of tracks:
    //correct slopes are taken from the compact
    //three track vertex reconstruction routine
    double Vertex_pi_e1[3]= {0.};
    double cda_pi_e1 = 0;
    double Vertex_pi_e2[3]= {0.};
    double cda_pi_e2 = 0;
    double Vertex_e1_e2[3]= {0.};
    double cda_e1_e2 = 0;
    if(!closap_double_(pion.Position,electron1.Position,pion.Slopes,electron1.Slopes,&cda_pi_e1,Vertex_pi_e1) ||
       !closap_double_(pion.Position,electron2.Position,pion.Slopes,electron2.Slopes,&cda_pi_e2,Vertex_pi_e2) ||
       !closap_double_(electron1.Position,electron2.Position,electron1.Slopes,electron2.Slopes,&cda_e1_e2,Vertex_e1_e2)){
        return 0;
    }


    Cuts cutting = make_cuts(dir1, pion, electron1, electron2, Vertex_pi_e1, Vertex_pi_e2, Vertex_e1_e2, "pinuee");
    ////Cuts
    //-- CUT1 DCH Geometry and Time Cut --
    if(cutting.DCH_Radius_pi < 14  ||
       cutting.DCH_Radius_pi > 110 ||
       cutting.DCH_Radius_el1 < 14 ||
       cutting.DCH_Radius_el1 > 110||
       cutting.DCH_Radius_el2 < 14 ||
       cutting.DCH_Radius_el2 > 110
       ){return 0;}

    if(cutting.Lkr_cut_el1  < 15.     ||
       cutting.Lkr_cut_el2  < 15.     ||
       cutting.Lkr_cut_pi   < 15.
       ){return 0;}
    if(   fabs(cutting.Lkr_x_el1) > 113 ||
          fabs(cutting.Lkr_x_el2) > 113 ||
          fabs(cutting.Lkr_x_pi ) > 113 ||
          fabs(cutting.Lkr_y_el1) > 113 ||
          fabs(cutting.Lkr_y_el2) > 113 ||
          fabs(cutting.Lkr_y_pi ) > 113
          ) {return 0;}
    if( (fabs(cutting.Lkr_x_el1) + fabs(cutting.Lkr_y_el1) ) > 159.8 ||
        (fabs(cutting.Lkr_x_el2) + fabs(cutting.Lkr_y_el2) ) > 159.8 ||
        (fabs(cutting.Lkr_x_pi) +  fabs(cutting.Lkr_y_pi ) ) > 159.8
        ){return 0;}


    if(IS_DATA)
        if(fabs(cutting.DCH_e1e2) > 10. ||
           fabs(cutting.DCH_pie1) > 10. ||
           fabs(cutting.DCH_pie2) > 10. ||
           fabs(cutting.Hod_e1e2) > 2.  ||
           fabs(cutting.Hod_pie1) > 2.  ||
           fabs(cutting.Hod_pie2) > 2.
           ){return 0;}
    //-- CUT1 DCH Geometry and Time Cut --

    Filling(sevt, sbur, dir1, ipi, iel1, iel2, igamma);
    FillK2pi(sevt, dir1, pion.Momentum, electron1.Momentum, electron2.Momentum, gamma.Momentum);
    dir1->FillVertexHist(Vertex_pi_e1, cda_pi_e1 , Vertex_pi_e2, cda_pi_e2, Vertex_e1_e2, cda_e1_e2,"k2pid");
    dir1->FillHist(electron1,"electron1");
    dir1->FillHist(electron2,"electron1");
    dir1->FillHist(pion,electron1,"pie1");
    dir1->FillHist(pion,electron2,"pie2");
    dir1->FillHist(electron1,electron2,"e1e2");
    dir1->fh_Kaon_Charge    ->Fill(Kcharge);
    dir1->fh_el1_Charge     ->Fill(sevt->track[iel1].q);
    dir1->fh_el2_Charge     ->Fill(sevt->track[iel2].q);
    dir1->fh_pion_Charge    ->Fill(sevt->track[ipi].q);
    //FillK2pi(sevt,dir1,pion.Momentum,electron1.Momentum,electron2.Momentum,gamma.Momentum);
    //--End of Vertex Reconstruction --

    //if(lda3_e1 < 0.8 || lda3_e2 < 0.8){return 0;}


    return 0;
}

//Calculates and returns blue tube corrections for slopes of tracks on x-y plane
void cleanup()    {


    return;
}
//Function that takes as input the charged particles, together with the cut
//directory and reconstructed vertexes dependent on the selection (munuee or K3pi)
//and fills the variables of the class Cuts with all the variables that it would be cut on.
//For K3pi selection - muon(pion1), electron1(pion2), electron2(pion3)
Cuts make_cuts(Hist_dir* dir1,Charged_Particle& pion, Charged_Particle& electron1, Charged_Particle& electron2,double* Vertex_pi_e1,double* Vertex_pi_e2,double* Vertex_e1_e2,std::string decay_type){
    Cuts cutting;

    cutting.DCH_e1e2 = electron2.GetDCHtime() - electron1.GetDCHtime();
    cutting.DCH_pie1 = pion.GetDCHtime() - electron1.GetDCHtime();
    cutting.DCH_pie2 = pion.GetDCHtime() - electron2.GetDCHtime();
    cutting.Hod_e1e2 = electron2.GetHodTime() - electron1.GetHodTime();
    cutting.Hod_pie1 = pion.GetHodTime() - electron1.GetHodTime();
    cutting.Hod_pie2 = pion.GetHodTime() - electron2.GetHodTime();
    cutting.cluster_e1e2 = electron1.GetClusterTime() - electron2.GetClusterTime();
    //cutting.cluster_pie1 = pion.GetClusterTime() - electron1.GetClusterTime();
    //cutting.cluster_pie2 = pion.GetClusterTime() - electron2.GetClusterTime();
    if(decay_type.compare("pinuee") == 0 ){
        cutting.zvtx_pie1_pie2 = Vertex_pi_e1[2] - Vertex_pi_e2[2];
        cutting.zvtx_pie1_e1e2 = Vertex_pi_e1[2] - Vertex_e1_e2[2];
        cutting.zvtx_pie2_e1e2 = Vertex_pi_e2[2] - Vertex_e1_e2[2];
        cutting.yvtx_pie1_pie2 = Vertex_pi_e1[1] - Vertex_pi_e2[1];
        cutting.yvtx_pie1_e1e2 = Vertex_pi_e1[1] - Vertex_e1_e2[1];
        cutting.yvtx_pie2_e1e2 = Vertex_pi_e2[1] - Vertex_e1_e2[1];
        cutting.xvtx_pie1_pie2 = Vertex_pi_e1[0] - Vertex_pi_e2[0];
        cutting.xvtx_pie1_e1e2 = Vertex_pi_e1[0] - Vertex_e1_e2[0];
        cutting.xvtx_pie2_e1e2 = Vertex_pi_e2[0] - Vertex_e1_e2[0];
        cutting.DCH_Radius_pi  = pion.GetDCHradius();
        cutting.DCH_Radius_el1 = electron1.GetDCHradius();
        cutting.DCH_Radius_el2 = electron2.GetDCHradius();

        cutting.Lkr_x_el1 = electron1.extrapolated_track_Lkr[0];
        cutting.Lkr_x_el2 = electron2.extrapolated_track_Lkr[0];
        cutting.Lkr_x_pi  = pion.extrapolated_track_Lkr[0];
        cutting.Lkr_y_el1 = electron1.extrapolated_track_Lkr[1];
        cutting.Lkr_y_el2 = electron2.extrapolated_track_Lkr[1];
        cutting.Lkr_y_pi  = pion.extrapolated_track_Lkr[1];
        cutting.Lkr_cut_el1 = electron1.GetLkrRadius();
        cutting.Lkr_cut_el2 = electron2.GetLkrRadius();
        cutting.Lkr_cut_pi  = pion.GetLkrRadius();

    }

    return cutting;
}



void FillK2pi(superCmpEvent* sevt, Hist_dir* dir, TLorentzVector pion,TLorentzVector el1,TLorentzVector el2,TLorentzVector gamma){

    TLorentzVector Pi0d_Momentum;
    TLorentzVector Kaon_Momentum;
    //int Kch =;
    Pi0d_Momentum = el1 + el2 + gamma;
    Kaon_Momentum = pion + el1 + el2 + gamma;
    dir->fh_Ntracks->Fill(sevt->Ntrack);
    dir->fh_Nclusters->Fill(sevt->Ncluster);
    dir->fh_Nvtx->Fill(sevt->Nvtx);

    dir->fh_gamma_momentum ->Fill(gamma.P());
    dir->fh_pi_momentum    ->Fill(pion.P());
    dir->fh_el1_momentum   ->Fill(el1.P());
    dir->fh_el2_momentum   ->Fill(el2.P());
    dir->fh_k2pi0d_P       ->Fill(Kaon_Momentum.P());
    dir->fh_k2pi0d_Pt      ->Fill(Kaon_Momentum.Pt());
    dir->fh_k2pi0d_M       ->Fill(Kaon_Momentum.M());
    dir->fh_pi0d_P         ->Fill(Pi0d_Momentum.P());
    dir->fh_pi0d_Pt        ->Fill(Pi0d_Momentum.Pt());
    dir->fh_pi0d_M         ->Fill(Pi0d_Momentum.M());

    //dir->fh_muee_M->Fill(Pi0d_Momentum.M());
    //dir->fh_Muee_M_3pi_assumption->Fill(Kaon_Momentum.M());
    //dir->fh_muee_Pt->Fill(Kaon_Momentum.Pt());
    //dir->fh_muee_P->Fill(Kaon_Momentum.P());
}

void Filling(superCmpEvent* sevt,superBurst* sbur, Hist_dir* dir,int ipi,int iel1,int iel2,int igamma){
    double EovP_el1   = sevt->cluster[sevt->track[iel1].iClus].energy/sevt->track[iel1].p;
    int charge_el1 = sevt->track[iel1].q;
    double EovP_el2   = sevt->cluster[sevt->track[iel2].iClus].energy/sevt->track[iel2].p;
    int charge_el2 = sevt->track[iel2].q;
    double EovP_pi    = sevt->cluster[sevt->track[ipi].iClus].energy/sevt->track[ipi].p  ;
    float lda3_e1;
    float lda3_e2;



    dir->fh_EoP_el1->Fill (EovP_el1);
    dir->fh_EoP_el2->Fill (EovP_el2);
    dir->fh_EoP_pion->Fill(EovP_pi);

    if(EovP_el2 > 0.95 && EovP_el2 < 1.05){
        lda3_e2 = ldaC(sevt,sbur,iel2);
        dir->fh_lda3_e2->Fill(lda3_e2);
    }



    if (charge_el1 == -1){
        dir->fh_el1_minus_Charge->Fill (charge_el1);
        dir->fh_EoP_el1_minus->Fill (EovP_el1);
        dir->fh_EoP_vs_p_el1_minus->Fill (EovP_el1,sevt->track[iel1].p);
        if(EovP_el1 > 0.95 && EovP_el1 < 1.05){
            lda3_e1 = ldaC(sevt,sbur,iel1);
            dir->fh_lda3_e1_minus->Fill (lda3_e1);
            dir->fh_lda3_vs_p_e1_minus->Fill (lda3_e1,sevt->track[iel1].p);
        }

    }

    if (charge_el1 == 1){
        dir->fh_el1_plus_Charge->Fill (charge_el1);
        dir->fh_EoP_el1_plus->Fill (EovP_el1);
        dir->fh_EoP_vs_p_el1_plus->Fill (EovP_el1,sevt->track[iel1].p);
        if(EovP_el1 > 0.95 && EovP_el1 < 1.05){
            lda3_e1 = ldaC(sevt,sbur,iel1);
            dir->fh_lda3_e1_plus->Fill (lda3_e1);
            dir->fh_lda3_vs_p_e1_plus->Fill (lda3_e1,sevt->track[iel1].p);
        }

    }

    dir->fh_EoP_vs_p_el1->Fill (EovP_el1,sevt->track[iel1].p);
    if(EovP_el1 > 0.95 && EovP_el1 < 1.05){
        lda3_e1 = ldaC(sevt,sbur,iel1);
        dir->fh_lda3_e1->Fill(lda3_e1);
        dir->fh_lda3_vs_p_e1->Fill (lda3_e1,sevt->track[iel1].p);
    }


}

double ldaMC(float p, float lda){
    double result;
    if(lda==0.85)
        result=0.498379+0.137042*p-0.0151076*pow(p,2.)+0.000751304*pow(p,3.)-1.06634e-05*pow(p,4.)-4.62566e-07*pow(p,5.)+1.91446e-08*pow(p,6.)-2.00753e-10*pow(p,7.);
    else if(lda==0.9)
        result=0.410062+0.166293*p-0.0209528*pow(p,2.)+0.00144389*pow(p,3.)-5.79143e-05*pow(p,4.)+1.34417e-06*pow(p,5.)-1.67e-08*pow(p,6.)+8.57178e-11*pow(p,7.);
    else if(lda==0.95)
        result=0.362081+0.144883*p-0.0140851*pow(p,2.)+0.000622962*pow(p,3.)-6.70367e-06*pow(p,4.)-4.22388e-07*pow(p,5.)+1.51773e-08*pow(p,6.)-1.49003e-10*pow(p,7.);
    else
        result=0;
    return result;
}


double ldaC(superCmpEvent *sevt,superBurst *sbur,int i) {
    double lda3;
    int itr, iclu,n;
    double eovp, geo_lkrZ, Zdch, pit, rmsr, xlkr, ylkr, distx, disty, disn;
    double RIN1, out1;
    int MC = 0;
    if(sbur->brtype == 2){MC = 1;}
    n = sevt->Ntrack;
    float x1[n], y1[n];
    geo_lkrZ = Geom->Lkr.z;
    Zdch = Geom->DCH.z;

    lda3 = 0.;
    iclu = sevt->track[i].iClus;
    if(iclu>-1){
        pit = p_corr_ab(sevt->track[i].p, sevt->track[i].q);
        eovp = sevt->cluster[iclu].energy/pit;
        rmsr = sqrt(sevt->cluster[iclu].rmsx*sevt->cluster[iclu].rmsx+sevt->cluster[iclu].rmsy*sevt->cluster[iclu].rmsy);
        if(MC == 1) {
            //  z1[i] = Geom->Lkr.z+16.5+4.3*log(sevt->cluster[i].energy);
            x1[iclu] = (sevt->cluster[iclu].x-0.013);//*(1+(z1[i]-Geom->Lkr.z)/10998);
            y1[iclu] = sevt->cluster[iclu].y;//*(1+(z1[i]-Geom->Lkr.z)/10998);
        }
        else if(MC == 0) {
            //  z1[i] = Geom->Lkr.z+16.5+4.3*log(sevt->cluster[i].energy);
            x1[iclu] = (sevt->cluster[iclu].x+0.136+0.00087*sevt->cluster[iclu].y);//*(1+(z1[i]-Geom->Lkr.z)/10998);
            y1[iclu] = (sevt->cluster[iclu].y+0.300+0.00087*sevt->cluster[iclu].x);//*(1+(z1[i]-Geom->Lkr.z)/10998);
        }
        xlkr = sevt->track[i].x+sevt->track[i].dxdz*(geo_lkrZ-Zdch);
        ylkr = sevt->track[i].y+sevt->track[i].dydz*(geo_lkrZ-Zdch);
        distx = xlkr-x1[iclu];
        disty = ylkr-y1[iclu];
        disn = sqrt(pit*(distx*distx+disty*disty));

        RIN1 = 39.49729-82.20312*eovp+32.76240*rmsr+1.597024*disn;
        //cout<<"Pit: "<<pit<<endl<<"Iclu: "<<iclu<<endl<<"eovp: "<<eovp<<endl<<"rmsr: "<<rmsr<<endl<<"disn: "<<disn<<endl<<"-------------"<<endl;

        out1 = sigmoidl(RIN1);
        lda3 = 1.000964-0.9863838*out1;
    }
    //cout<<"lda3= "<<lda3<<endl;
    return lda3;
}

float sigmoidl(double x){
    float sigmoidl;

    if ( x > 37.){
        sigmoidl = -1;
    } else if ( x < -37){
        sigmoidl = 0;
    } else {
        sigmoidl = 1./(1. + exp(-x));
    }
    return sigmoidl;

}
