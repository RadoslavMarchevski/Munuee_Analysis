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
#include "Hist_dir.h"
#include "Cuts.h"
using namespace std;

vector<TVector3*> clcorrpos;
TVector3* clpos;

void cleanup();
Cuts make_cuts(Hist_dir* dir1, Charged_Particle& muon, Charged_Particle& electron1, Charged_Particle& electron2,double* Vertex_mu_e1,double* Vertex_mu_e2,double* Vertex_e1_e2);


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
    int pi1=-1;
    int pi2=-1;
    int pi3=-1;


    bool LKrCalCorr=true;

    //--------------------------------------------------------//
    //Definition of Variables
    int    ntrack = sevt->Ntrack; // number of tracks
    int    nclust = sevt->Ncluster; // number of tracks
    int    nvtx   = sevt->Nvtx;   // number of vtx
    double DCHbz = Geom->DCH.bz;         // z before magnet
    double DCHz = Geom->DCH.z;         // z after magnet
    double LKrz=Geom->Lkr.z;
    double COmPaCt_Z_Vertex = sevt->vtx[0].z;
    static double massKaonC = 0.493677;
    const double Electron_EoverP_up = 1.05;
    const double Electron_EoverP_down = 0.95;
    const double Muon_EoverP_up = 0.95;

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
    double Track_Momentum;
    double Track_Quality;
    double Track_DeadCell_Distance;
    double Track_Energy;
    double Track_EoverP;
    /// apply non linearity corrections
    /// Correct LKr non-linearity (effective for energies < 11GeV)
    if (LKrCalCorr && IS_DATA) {
        user_lkrcalcor_SC (sbur,sevt,1);
    }

    if (IS_DATA)
        if( sbur->BadB.Lkr == 1|| sbur->BadB.Dch == 1||
            sbur->BadB.Mbx == 1|| sbur->BadB.Muv == 1||
            sbur->BadB.HodC== 1|| sbur->BadB.Phys== 1){
            sbur->BadB.Skip=1;
        }
    Initial_dir->fh_Ntracks->Fill(ntrack);
    Initial_dir->fh_Nvtx->Fill(nvtx);
    Initial_dir->fh_COmPaCt_Z_Vertex->Fill(COmPaCt_Z_Vertex);
    if(nvtx != 1){return 0;}
    if(ntrack != 3 ) {return 0;}
    if(COmPaCt_Z_Vertex < -2000. || COmPaCt_Z_Vertex > 8000.){return 0;}

    //Looping over the number of tracks for each event
    for (int i=0; i<ntrack; i++) {
        Kcharge += sevt->track[i].q;
        Track_Momentum          = sevt->track[i].p;
        Track_Quality           = sevt->track[i].quality;
        Track_DeadCell_Distance = sevt->track[i].dDeadCell;
        Track_icl               = sevt->track[i].iClus;
        Track_imu               = sevt->track[i].iMuon;
        Track_Energy  = sevt->cluster[Track_icl].energy;
        Track_EoverP  = Track_Energy / Track_Momentum;

        ngoodtrack++;
        Initial_dir->fh_eop->Fill(Track_EoverP);
        Initial_dir->fh_Track_Momentum->Fill(Track_Momentum);
        //Particle identification for K to mu nu e e
        if(IS_DATA){
            if(Track_EoverP > Electron_EoverP_down && Track_EoverP < Electron_EoverP_up){
                if(iel1<0){
                    iel1 = i;
                } else if(iel2<0){
                    iel2 = i;
                } else {;}
                Nelectrons++;
            } else if (imu < 0 && Track_imu!= -1 && Track_EoverP < Muon_EoverP_up ){
                imu = i;
                Nmuons++;
            }


            //Particle identification for K to 3 pi charged
            if( Track_Momentum > 10. ){
                if(pi1<0){
                    pi1 = i;
                } else if(pi2<0){
                    pi2 = i;
                } else if(pi3<0){
                    pi3 = i;
                }

                NK3pi_pions++;
            }//ENDIF K3pi selection


        }// ENDIF IS_DATA
    }//end for i
    //Kmunuee selection
    Charged_Particle muon(sevt,sbur,-13,imu);
    Charged_Particle electron1(sevt,sbur,-11,iel1);
    Charged_Particle electron2(sevt,sbur,-11,iel2);

    //K3pi selection
    Charged_Particle pion1(sevt,sbur,211,pi1);
    Charged_Particle pion2(sevt,sbur,211,pi2);
    Charged_Particle pion3(sevt,sbur,211,pi3);

    if(Nelectrons+Nmuons==3)    {
        if(Nelectrons==2 && Nmuons==1 ){
            if((muon.GetCharge() + electron1.GetCharge() + electron2.GetCharge() )==3.)    {// mu+, e+, e+
                Event_Type = 1;
            } else if ((muon.GetCharge() + electron1.GetCharge() + electron2.GetCharge())== -3.){ // mu-, e-, e-
                Event_Type = 5;
            } else if ((muon.GetCharge() + electron1.GetCharge() + electron2.GetCharge())== 1. && muon.GetCharge() == 1. ){// mu+, e+, e-
                Event_Type = 0;
            } else if ((muon.GetCharge() + electron1.GetCharge() + electron2.GetCharge())== 1. && muon.GetCharge() == -1.){// mu-, e+, e+
                Event_Type = 4;
            } else if ((muon.GetCharge() + electron1.GetCharge() + electron2.GetCharge())== -1.&& muon.GetCharge() == -1.){// mu-, e+, e-
                Event_Type = 3;
            } else if ((muon.GetCharge() + electron1.GetCharge() + electron2.GetCharge())== -1.&& muon.GetCharge() == 1. ){// mu+, e-, e-
                Event_Type = 2;
            }
        }
    }

    if(NK3pi_pions==3){
        if((pion1.GetCharge() + pion2.GetCharge() + pion3.GetCharge() )==3.)    {// pi+, pi+, pi+
            K3pi_Event_Type = 1;
        } else if ((pion1.GetCharge() + pion2.GetCharge() + pion3.GetCharge())== -3.){ // pi-, pi-, pi-
            K3pi_Event_Type = 5;
        } else if ((pion1.GetCharge() + pion2.GetCharge() + pion3.GetCharge())== 1. ){// pi+, pi+, pi-
            K3pi_Event_Type = 0;
        } else if ((pion1.GetCharge() + pion2.GetCharge() + pion3.GetCharge())== -1.){// pi+, pi-, pi-
            K3pi_Event_Type = 2;
        }
    }
    //printf("abcog_params.cogX1p=%f, abcog_params.cogX1n=%f\n",abcog_params.cogX1p,abcog_params.cogX1n);
    //printf("abcog_params.cogX4p=%f,abcog_params.cogX4n=%f\n",abcog_params.cogX4p,abcog_params.cogX4n);
    //printf("abcog_params.cogY1p=%f, abcog_params.cogY1n=%f\n",abcog_params.cogY1p,abcog_params.cogY1n);
    //printf("abcog_params.cogY4p=%f,abcog_params.cogY4n=%f\n",abcog_params.cogY4p,abcog_params.cogY4n);


    //Normalization channel K3pi selection
    if(pi1 != -1 && pi2 != -1 && pi3 != -1){

        K3pi_selection->fh_Ntracks->Fill(ntrack);
        K3pi_selection->fh_Nvtx->Fill(nvtx);
        K3pi_selection->fh_Kaon_Charge->Fill(Kcharge);
        K3pi_selection->fh_Event_Type->Fill(K3pi_Event_Type);
        K3pi_selection->fh_COmPaCt_Z_Vertex->Fill(COmPaCt_Z_Vertex);
        K3pi_selection->FillHist(pion1,"pion1");
        K3pi_selection->FillHist(pion2,"pion2");
        K3pi_selection->FillHist(pion3,"pion3");

        //K3pi_selection->fh_Pion_Momentum->Fill(pion1.GetMomentum());
        //K3pi_selection->fh_Pion_Momentum->Fill(pion2.GetMomentum());
        //K3pi_selection->fh_Pion_Momentum->Fill(pion3.GetMomentum());
        //K3pi_selection->fh_eop->Fill(pion1.GetEnergyLeftInEcal()/ pion1.GetMomentum() );
        //K3pi_selection->fh_eop->Fill(pion2.GetEnergyLeftInEcal()/ pion2.GetMomentum() );
        //K3pi_selection->fh_eop->Fill(pion3.GetEnergyLeftInEcal()/ pion3.GetMomentum() );


        if(Kcharge*pion1.GetCharge()==-1)
            K3pi_selection->fh_odd_eop->Fill(pion1.GetEnergyLeftInEcal()/ pion1.GetMomentum());
        else if(Kcharge*pion2.GetCharge()==-1)
            K3pi_selection->fh_odd_eop->Fill(pion2.GetEnergyLeftInEcal()/ pion2.GetMomentum());
        else if(Kcharge*pion3.GetCharge()==-1)
            K3pi_selection->fh_odd_eop->Fill(pion3.GetEnergyLeftInEcal()/ pion3.GetMomentum());

    }

    dir1->fh_Nvtx->Fill(nvtx);
    dir1->fh_Ntracks->Fill(ntrack);
    dir1->fh_Event_Type->Fill(Event_Type);
    dir1->fh_Kaon_Charge->Fill(Kcharge);
    //Signal particle identification
    if(imu == -1 || iel1 == -1 || iel2 == -1){return 0;}

    //------------------------------------Cut one-------------------------------------


    if(Kcharge*electron1.GetCharge()==-1)
        dir1->fh_odd_eop->Fill(electron1.GetEnergyLeftInEcal()/ electron1.GetMomentum());
    else if(Kcharge*electron2.GetCharge()==-1)
        dir1->fh_odd_eop->Fill(electron2.GetEnergyLeftInEcal()/ electron2.GetMomentum());

    dir1->fh_COmPaCt_Z_Vertex->Fill(COmPaCt_Z_Vertex);
    //--Vertex Reconstruction --
    //Vertex reconstruction for each pair of tracks:
    //correct slopes are taken from the compact
    //three track vertex reconstruction routine
    double Vertex_mu_e1[3]= {0.};
    double cda_mu_e1 = 0;
    double Vertex_mu_e2[3]= {0.};
    double cda_mu_e2 = 0;
    double Vertex_e1_e2[3]= {0.};
    double cda_e1_e2 = 0;
    bool vtx_mue1 = false;
    bool vtx_mue2 = false;
    bool vtx_e1e2 = false;
    if(!closap_double_(muon.Position,electron1.Position,muon.Slopes,electron1.Slopes,&cda_mu_e1,Vertex_mu_e1) ||
       !closap_double_(muon.Position,electron2.Position,muon.Slopes,electron2.Slopes,&cda_mu_e2,Vertex_mu_e2) ||
       !closap_double_(electron1.Position,electron2.Position,electron1.Slopes,electron2.Slopes,&cda_e1_e2,Vertex_e1_e2)){
        //vtx_mue1 = true;
        //vtx_mue2 = true;
        //vtx_e1e2 = true;
        return 0;
    }

    dir1->fh_Z_Vertex->Fill(COmPaCt_Z_Vertex );

    dir1->fh_zvtxdiff_mue1_mue2->Fill( Vertex_mu_e1[2] - Vertex_mu_e2[2]);
    dir1->fh_zvtxdiff_mue1_e1e2->Fill( Vertex_mu_e1[2] - Vertex_e1_e2[2]);
    dir1->fh_zvtxdiff_mue2_e1e2->Fill( Vertex_mu_e2[2] - Vertex_e1_e2[2]);
    dir1->fh_yvtxdiff_mue1_mue2->Fill( Vertex_mu_e1[1] - Vertex_mu_e2[1]);
    dir1->fh_yvtxdiff_mue1_e1e2->Fill( Vertex_mu_e1[1] - Vertex_e1_e2[1]);
    dir1->fh_yvtxdiff_mue2_e1e2->Fill( Vertex_mu_e2[1] - Vertex_e1_e2[1]);
    dir1->fh_xvtxdiff_mue1_mue2->Fill( Vertex_mu_e1[0] - Vertex_mu_e2[0]);
    dir1->fh_xvtxdiff_mue1_e1e2->Fill( Vertex_mu_e1[0] - Vertex_e1_e2[0]);
    dir1->fh_xvtxdiff_mue2_e1e2->Fill( Vertex_mu_e2[0] - Vertex_e1_e2[0]);
    //--End of Vertex Reconstruction --

    //Fill the properties of single tracks (Momentum, E/p , Time variables etc ..)
    dir1->FillHist(muon,"muon");
    dir1->FillHist(electron1,"electron1");
    dir1->FillHist(electron2,"electron2");
    //Fill histograms for each pair of tracks. First calculates all necessary
    //two track variables like time differences between two tracks e+e-
    //invariant mass etc. . Then fills the histograms with the calculated
    //variables.
    dir1->FillHist(muon,electron1,"mue1");
    dir1->FillHist(muon,electron2,"mue2");
    dir1->FillHist(electron1,electron2,"e1e2");
    //Fill histograms for the three tracks. Calculates three track vector
    //and fills the histograms with the variables of interest (momentum,mass ..).
    dir1->FillHist(electron1,electron2,muon);

    //Producing cut variable in more readable way with the class
    //described in Cuts.h;
    Cuts cutting = make_cuts(dir1,muon,electron1,electron2,Vertex_mu_e1,Vertex_mu_e2,Vertex_e1_e2);

    //Defining variables that it would be cut on
    //
    ////Cuts
    //-- Momentum cut ---

    if(cutting.Mu_P < 10. || cutting.Mu_P > 50.){return 0;}
    if(cutting.E1_P < 3.  || cutting.E1_P > 50.){return 0;}
    if(cutting.E2_P < 3.  || cutting.E2_P > 50.){return 0;}
    if(cutting.muee_P < 44 || cutting.muee_P > 66){return 0;}
    //--ENDOF Momentum cut ---
    //-- Timing cut ---
    if(fabs(cutting.DCH_e1e2) > 10. ||
       fabs(cutting.DCH_mue1) > 10. ||
       fabs(cutting.DCH_mue2) > 10. ||
       fabs(cutting.Hod_e1e2) > 2. ||
       fabs(cutting.Hod_mue1) > 2. ||
       fabs(cutting.Hod_mue2) > 2.
        ){return 0;}
    //-- ENDOF Timing cut ---
    //--Vertex Cut --
    if(fabs(cutting.zvtx_e1e2) > 500 ||
       fabs(cutting.zvtx_mue1) > 500 ||
       fabs(cutting.zvtx_mue2) > 500 ||
       fabs(cutting.xvtx_e1e2) > 500 ||
       fabs(cutting.xvtx_mue1) > 500 ||
       fabs(cutting.xvtx_mue2) > 500 ||
       fabs(cutting.yvtx_e1e2) > 500 ||
       fabs(cutting.yvtx_mue2) > 500 ||
       fabs(cutting.yvtx_mue1) > 500
        ){return 0;}
    //--ENDOF Vertex Cut --

    dir2->FillHist(muon,"muon");
    dir2->FillHist(electron1,"electron1");
    dir2->FillHist(electron2,"electron2");
    dir2->FillHist(muon,electron1,"mue1");
    dir2->FillHist(muon,electron2,"mue2");
    dir2->FillHist(electron1,electron2,"e1e2");
    dir2->FillHist(electron1,electron2,muon);

    return 0;
}

//Calculates and returns blue tube corrections for slopes of tracks on x-y plane
void cleanup()    {


    return;
}

Cuts make_cuts(Hist_dir* dir1,Charged_Particle& muon, Charged_Particle& electron1, Charged_Particle& electron2,double* Vertex_mu_e1,double* Vertex_mu_e2,double* Vertex_e1_e2){
    Cuts cutting;

    cutting.DCH_e1e2 = electron2.GetDCHtime() - electron1.GetDCHtime();
    cutting.DCH_mue1 = muon.GetDCHtime() - electron1.GetDCHtime();
    cutting.DCH_mue2 = muon.GetDCHtime() - electron2.GetDCHtime();
    cutting.Hod_e1e2 = electron2.GetHodTime() - electron1.GetHodTime();
    cutting.Hod_mue1 = muon.GetHodTime() - electron1.GetHodTime();
    cutting.Hod_mue2 = muon.GetHodTime() - electron2.GetHodTime();
    cutting.Mu_P     = muon.GetMomentum();
    cutting.E1_P     = electron1.GetMomentum();
    cutting.E2_P     = electron2.GetMomentum();
    cutting.muee_P   = dir1->GetThreeTrackMomentum().P();
    cutting.muee_Pt  = dir1->GetThreeTrackMomentum().Pt();
    cutting.zvtx_e1e2= Vertex_mu_e1[2] - Vertex_mu_e2[2];
    cutting.zvtx_mue2= Vertex_mu_e1[2] - Vertex_e1_e2[2];
    cutting.zvtx_mue1= Vertex_mu_e2[2] - Vertex_e1_e2[2];
    cutting.yvtx_e1e2= Vertex_mu_e1[1] - Vertex_mu_e2[1];
    cutting.yvtx_mue2= Vertex_mu_e1[1] - Vertex_e1_e2[1];
    cutting.yvtx_mue1= Vertex_mu_e2[1] - Vertex_e1_e2[1];
    cutting.xvtx_e1e2= Vertex_mu_e1[0] - Vertex_mu_e2[0];
    cutting.xvtx_mue2= Vertex_mu_e1[0] - Vertex_e1_e2[0];
    cutting.xvtx_mue1= Vertex_mu_e2[0] - Vertex_e1_e2[0];
    return cutting;
}
