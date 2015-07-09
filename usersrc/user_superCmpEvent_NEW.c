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
Cuts make_cuts(Hist_dir* dir1, Charged_Particle& muon, Charged_Particle& electron1, Charged_Particle& electron2,double* Vertex_mu_e1,double* Vertex_mu_e2,double* Vertex_e1_e2, std::string decay_type);


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
    int    nclust = sevt->Ncluster; // number of clusters
    int    nvtx   = sevt->Nvtx;   // number of vtx
    double DCHbz = Geom->DCH.bz;         // z before magnet
    double DCHz = Geom->DCH.z;         // z after magnet
    double LKrz=Geom->Lkr.z;
    double COmPaCt_Z_Vertex = sevt->vtx[0].z;
    static double massKaonC = 0.493677;
    const double Electron_EoverP_up = 1.05;
    const double Electron_EoverP_down = 0.95;
    const double Muon_EoverP_up = 0.3;

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
    Initial_dir->FillCommonHist(sevt);

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
        double Vertex_pi1_pi2[3]= {0.};
        double cda_pi1_pi2 = 0;
        double Vertex_pi1_pi3[3]= {0.};
        double cda_pi1_pi3 = 0;
        double Vertex_pi2_pi3[3]= {0.};
        double cda_pi2_pi3 = 0;

//Checking if vertexes are true and if not go to the next event
        if(!closap_double_(pion1.Position,pion2.Position,pion1.Slopes,pion2.Slopes,&cda_pi1_pi2,Vertex_pi1_pi2) ||
           !closap_double_(pion1.Position,pion3.Position,pion1.Slopes,pion3.Slopes,&cda_pi1_pi3,Vertex_pi1_pi3) ||
           !closap_double_(pion2.Position,pion3.Position,pion2.Slopes,pion3.Slopes,&cda_pi2_pi3,Vertex_pi2_pi3)){
            return 0;
        }
        K3pi_selection->ComputeThreeTrack(pion1,pion2,pion3);
        Cuts cut_k3pi = make_cuts(K3pi_selection,pion1,pion2,pion3,Vertex_pi1_pi2,Vertex_pi1_pi3,Vertex_pi2_pi3,"K3pi");
        ////Cuts
        //73-3-10
        //-- CUT1 Momentum cut ---
        if(cut_k3pi.muee_P > 54. &&
           cut_k3pi.muee_P < 66. &&

           //--ENDOF CUT1 Momentum cut ---
           //-- CUT2 Timing cut ---
           fabs(cut_k3pi.DCH_e1e2) < 10.&&
           fabs(cut_k3pi.DCH_mue1) < 10.&&
           fabs(cut_k3pi.DCH_mue2) < 10.&&
           fabs(cut_k3pi.Hod_e1e2) < 2. &&
           fabs(cut_k3pi.Hod_mue1) < 2. &&
           fabs(cut_k3pi.Hod_mue2) < 2. &&
           //-- ENDOF CUT2 Timing cut ---
           //-- CUT3 Vertex Cut --
           fabs(cut_k3pi.zvtx_pi1pi2_pi2pi3) < 500 &&
           fabs(cut_k3pi.zvtx_pi1pi2_pi1pi3) < 500 &&
           fabs(cut_k3pi.zvtx_pi1pi3_pi2pi3) < 500 &&
           fabs(cut_k3pi.xvtx_pi1pi2_pi2pi3) < 500 &&
           fabs(cut_k3pi.xvtx_pi1pi2_pi1pi3) < 500 &&
           fabs(cut_k3pi.xvtx_pi1pi3_pi2pi3) < 500 &&
           fabs(cut_k3pi.yvtx_pi1pi2_pi2pi3) < 500 &&
           fabs(cut_k3pi.yvtx_pi1pi2_pi1pi3) < 500 &&
           fabs(cut_k3pi.yvtx_pi1pi3_pi2pi3) < 500 &&
           //--ENDOF CUT3 Vertex Cut --
           //-- CUT4 Three track invariant mass Cut --
           cut_k3pi.muee_M > 0.4837 &&
           cut_k3pi.muee_M < 0.5037 &&
           //-- ENDOF CUT4 Three track invariant mass Cut --
           //-- CUT5 DCH geometry Cut --
           cut_k3pi.DCH_Radius_pi1 > 14. &&
           cut_k3pi.DCH_Radius_pi1 < 110.&&
           cut_k3pi.DCH_Radius_pi2 > 14. &&
           cut_k3pi.DCH_Radius_pi2 < 110.&&
           cut_k3pi.DCH_Radius_pi3 > 14. &&
           cut_k3pi.DCH_Radius_pi3 < 110.
           //-- ENDOF CUT5 DCH geometry Cut --
            ){

            K3pi_selection->FillCommonHist(sevt);
            K3pi_selection->FillHist(pion1,"pion1");
            K3pi_selection->FillHist(pion2,"pion2");
            K3pi_selection->FillHist(pion3,"pion3");
            K3pi_selection->FillHist(K3pi_selection->GetThreeTrackMomentum(),K3pi_selection->GetNuMomentum());
            K3pi_selection->fh_Kaon_Charge->Fill(Kcharge);
            K3pi_selection->fh_Event_Type->Fill(K3pi_Event_Type);


            K3pi_selection->FillVertexHist(Vertex_pi1_pi2, cda_pi1_pi2 , Vertex_pi2_pi3, cda_pi2_pi3, Vertex_pi1_pi3, cda_pi1_pi3,"K3pi");

            if(Kcharge*pion1.GetCharge()==-1)
                K3pi_selection->fh_odd_eop->Fill(pion1.GetEnergyLeftInEcal()/ pion1.GetMomentum());
            else if(Kcharge*pion2.GetCharge()==-1)
                K3pi_selection->fh_odd_eop->Fill(pion2.GetEnergyLeftInEcal()/ pion2.GetMomentum());
            else if(Kcharge*pion3.GetCharge()==-1)
                K3pi_selection->fh_odd_eop->Fill(pion3.GetEnergyLeftInEcal()/ pion3.GetMomentum());
        }
    }

//-------------------------KMUNUEE SIGNAL SELECTION ----------------------------------------
    //Signal particle identification
    if(imu == -1 || iel1 == -1 || iel2 == -1){return 0;}

    dir1->FillCommonHist(sevt);
    dir1->fh_Event_Type->Fill(Event_Type);
    dir1->fh_Kaon_Charge->Fill(Kcharge);
    if(Kcharge*electron1.GetCharge()==-1)
        dir1->fh_odd_eop->Fill(electron1.GetEnergyLeftInEcal()/ electron1.GetMomentum());
    else if(Kcharge*electron2.GetCharge()==-1)
        dir1->fh_odd_eop->Fill(electron2.GetEnergyLeftInEcal()/ electron2.GetMomentum());

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

//Checking if vertexes are true and if not go to the next event
    if(!closap_double_(muon.Position,electron1.Position,muon.Slopes,electron1.Slopes,&cda_mu_e1,Vertex_mu_e1) ||
       !closap_double_(muon.Position,electron2.Position,muon.Slopes,electron2.Slopes,&cda_mu_e2,Vertex_mu_e2) ||
       !closap_double_(electron1.Position,electron2.Position,electron1.Slopes,electron2.Slopes,&cda_e1_e2,Vertex_e1_e2)){
        return 0;
    }

    dir1->FillVertexHist(Vertex_mu_e1, cda_mu_e1 , Vertex_mu_e2, cda_mu_e2, Vertex_e1_e2, cda_e1_e2,"munuee");

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
    dir1->ComputeThreeTrack(electron1,electron2,muon);

    dir1->FillHist(dir1->GetThreeTrackMomentum(),dir1->GetNuMomentum());

    //Producing cut variable in more readable way with the class
    //described in Cuts.h;
    Cuts cutting = make_cuts(dir1,muon,electron1,electron2,Vertex_mu_e1,Vertex_mu_e2,Vertex_e1_e2,"munuee");

    //Defining variables that it would be cut on
    //
    ////Cuts
    //-- CUT1 DCH Geometry Cut --
    if(cutting.DCH_Radius_mu < 14 &&
       cutting.DCH_Radius_mu > 110&&
       cutting.DCH_Radius_el1 < 14 &&
       cutting.DCH_Radius_el1 > 110&&
       cutting.DCH_Radius_el2 < 14 &&
       cutting.DCH_Radius_el2 > 110
        ){return 0;}
    //-- CUT1 DCH Geometry Cut --
    dir3->fh_Event_Type->Fill(Event_Type);
    dir3->fh_Kaon_Charge->Fill(Kcharge);
    dir3->FillCommonHist(sevt);
    dir3->FillVertexHist(Vertex_mu_e1, cda_mu_e1 , Vertex_mu_e2, cda_mu_e2, Vertex_e1_e2, cda_e1_e2,"munuee");
    dir3->FillHist(muon,"muon");
    dir3->FillHist(electron1,"electron1");
    dir3->FillHist(electron2,"electron2");
    dir3->FillHist(muon,electron1,"mue1");
    dir3->FillHist(muon,electron2,"mue2");
    dir3->FillHist(electron1,electron2,"e1e2");
    dir3->ComputeThreeTrack(electron1,electron2,muon);
    dir3->FillHist(dir3->GetThreeTrackMomentum(),dir3->GetNuMomentum());

    //-- CUT2 Momentum cut ---
    if(cutting.Mu_P < 10. || cutting.Mu_P > 50.){return 0;}
    if(cutting.E1_P < 3.  || cutting.E1_P > 50.){return 0;}
    if(cutting.E2_P < 3.  || cutting.E2_P > 50.){return 0;}
    if(cutting.muee_P > 66){return 0;}
    //--ENDOF CUT2 Momentum cut ---

    dir4->fh_Event_Type->Fill(Event_Type);
    dir4->fh_Kaon_Charge->Fill(Kcharge);
    dir4->FillCommonHist(sevt);
    dir4->FillVertexHist(Vertex_mu_e1, cda_mu_e1 , Vertex_mu_e2, cda_mu_e2, Vertex_e1_e2, cda_e1_e2,"munuee");
    dir4->FillHist(muon,"muon");
    dir4->FillHist(electron1,"electron1");
    dir4->FillHist(electron2,"electron2");
    dir4->FillHist(muon,electron1,"mue1");
    dir4->FillHist(muon,electron2,"mue2");
    dir4->FillHist(electron1,electron2,"e1e2");
    dir4->ComputeThreeTrack(electron1,electron2,muon);
    dir4->FillHist(dir4->GetThreeTrackMomentum(),dir4->GetNuMomentum());

    //-- CUT3 Vertex Cut --
    if(fabs(cutting.zvtx_mue1_mue2) > 500 ||
       fabs(cutting.zvtx_mue1_e1e2) > 500 ||
       fabs(cutting.zvtx_mue2_e1e2) > 500 ||
       fabs(cutting.xvtx_mue1_mue2) > 500 ||
       fabs(cutting.xvtx_mue1_e1e2) > 500 ||
       fabs(cutting.xvtx_mue2_e1e2) > 500 ||
       fabs(cutting.yvtx_mue1_mue2) > 500 ||
       fabs(cutting.yvtx_mue1_e1e2) > 500 ||
       fabs(cutting.yvtx_mue2_e1e2) > 500
        ){return 0;}
    //--ENDOF CUT3 Vertex Cut --
    dir5->fh_Event_Type->Fill(Event_Type);
    dir5->fh_Kaon_Charge->Fill(Kcharge);
    dir5->FillCommonHist(sevt);
    dir5->FillVertexHist(Vertex_mu_e1, cda_mu_e1 , Vertex_mu_e2, cda_mu_e2, Vertex_e1_e2, cda_e1_e2,"munuee");
    dir5->FillHist(muon,"muon");
    dir5->FillHist(electron1,"electron1");
    dir5->FillHist(electron2,"electron2");
    dir5->FillHist(muon,electron1,"mue1");
    dir5->FillHist(muon,electron2,"mue2");
    dir5->FillHist(electron1,electron2,"e1e2");
    dir5->ComputeThreeTrack(electron1,electron2,muon);
    dir5->FillHist(dir5->GetThreeTrackMomentum(),dir5->GetNuMomentum());

    //-- CUT3 Transverse Momentum Cut --
    if(cutting.muee_Pt < 0.022 ){return 0;}
    //--ENDOF CUT3 Transverse Momentum Cut --
    dir6->fh_Event_Type->Fill(Event_Type);
    dir6->fh_Kaon_Charge->Fill(Kcharge);
    dir6->FillCommonHist(sevt);
    dir6->FillVertexHist(Vertex_mu_e1, cda_mu_e1 , Vertex_mu_e2, cda_mu_e2, Vertex_e1_e2, cda_e1_e2,"munuee");
    dir6->FillHist(muon,"muon");
    dir6->FillHist(electron1,"electron1");
    dir6->FillHist(electron2,"electron2");
    dir6->FillHist(muon,electron1,"mue1");
    dir6->FillHist(muon,electron2,"mue2");
    dir6->FillHist(electron1,electron2,"e1e2");
    dir6->ComputeThreeTrack(electron1,electron2,muon);
    dir6->FillHist(dir6->GetThreeTrackMomentum(),dir6->GetNuMomentum());

    //-- CUT4 Invariant Mass Cut --
    //if(cutting.mee < 0.140){return 0;}
    //--ENDOF CUT4 Invariant Mass Cut --
    dir7->fh_Event_Type->Fill(Event_Type);
    dir7->fh_Kaon_Charge->Fill(Kcharge);
    dir7->FillCommonHist(sevt);
    dir7->FillVertexHist(Vertex_mu_e1, cda_mu_e1 , Vertex_mu_e2, cda_mu_e2, Vertex_e1_e2, cda_e1_e2,"munuee");
    dir7->FillHist(muon,"muon");
    dir7->FillHist(electron1,"electron1");
    dir7->FillHist(electron2,"electron2");
    dir7->FillHist(muon,electron1,"mue1");
    dir7->FillHist(muon,electron2,"mue2");
    dir7->FillHist(electron1,electron2,"e1e2");
    dir7->ComputeThreeTrack(electron1,electron2,muon);
    dir7->FillHist(dir7->GetThreeTrackMomentum(),dir7->GetNuMomentum());
    //-- CUT5 Timing cut ---
    if(fabs(cutting.DCH_e1e2) > 10. ||
       fabs(cutting.DCH_mue1) > 10. ||
       fabs(cutting.DCH_mue2) > 10. ||
       fabs(cutting.Hod_e1e2) > 2. ||
       fabs(cutting.Hod_mue1) > 2. ||
       fabs(cutting.Hod_mue2) > 2.
        ){return 0;}
    //-- ENDOF CUT5 Timing cut ---

    dir8->fh_Event_Type->Fill(Event_Type);
    dir8->fh_Kaon_Charge->Fill(Kcharge);
    dir8->FillCommonHist(sevt);
    dir8->FillVertexHist(Vertex_mu_e1, cda_mu_e1 , Vertex_mu_e2, cda_mu_e2, Vertex_e1_e2, cda_e1_e2,"munuee");
    dir8->FillHist(muon,"muon");
    dir8->FillHist(electron1,"electron1");
    dir8->FillHist(electron2,"electron2");
    dir8->FillHist(muon,electron1,"mue1");
    dir8->FillHist(muon,electron2,"mue2");
    dir8->FillHist(electron1,electron2,"e1e2");
    dir8->ComputeThreeTrack(electron1,electron2,muon);
    dir8->FillHist(dir8->GetThreeTrackMomentum(),dir8->GetNuMomentum());

    //K3pi wrong sign selection
    if(electron1.GetCharge()==electron2.GetCharge() &&
       electron1.GetCharge()!=muon.GetCharge()
        ){
        dir2->fh_Event_Type->Fill(Event_Type);
        dir2->fh_Kaon_Charge->Fill(Kcharge);
        dir2->FillCommonHist(sevt);
        dir2->FillVertexHist(Vertex_mu_e1, cda_mu_e1 , Vertex_mu_e2, cda_mu_e2, Vertex_e1_e2, cda_e1_e2,"munuee");
        dir2->FillHist(muon,"muon");
        dir2->FillHist(electron1,"electron1");
        dir2->FillHist(electron2,"electron2");
        dir2->FillHist(muon,electron1,"mue1");
        dir2->FillHist(muon,electron2,"mue2");
        dir2->FillHist(electron1,electron2,"e1e2");
        dir2->ComputeThreeTrack(electron1,electron2,muon);
        dir2->FillHist(dir2->GetThreeTrackMomentum(),dir2->GetNuMomentum());
    }

    if(muon.GetCharge()*Kcharge != 1){return 0;}

    dir9->fh_Event_Type->Fill(Event_Type);
    dir9->fh_Kaon_Charge->Fill(Kcharge);
    dir9->FillCommonHist(sevt);
    dir9->FillVertexHist(Vertex_mu_e1, cda_mu_e1 , Vertex_mu_e2, cda_mu_e2, Vertex_e1_e2, cda_e1_e2,"munuee");
    dir9->FillHist(muon,"muon");
    dir9->FillHist(electron1,"electron1");
    dir9->FillHist(electron2,"electron2");
    dir9->FillHist(muon,electron1,"mue1");
    dir9->FillHist(muon,electron2,"mue2");
    dir9->FillHist(electron1,electron2,"e1e2");
    dir9->ComputeThreeTrack(electron1,electron2,muon);
    dir9->FillHist(dir9->GetThreeTrackMomentum(),dir9->GetNuMomentum());

    Charged_Particle munuee_pion1(sevt,sbur,211,imu);
    Charged_Particle munuee_pion2(sevt,sbur,211,iel1);
    Charged_Particle munuee_pion3(sevt,sbur,211,iel2);
    dir10->ComputeThreeTrack(munuee_pion1,munuee_pion2,munuee_pion3);
    dir10->FillHist(dir9->GetThreeTrackMomentum(),dir9->GetNuMomentum());

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
Cuts make_cuts(Hist_dir* dir1,Charged_Particle& muon, Charged_Particle& electron1, Charged_Particle& electron2,double* Vertex_mu_e1,double* Vertex_mu_e2,double* Vertex_e1_e2,std::string decay_type){
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
    if(decay_type.compare("munuee") == 0 ){
        cutting.mee      = dir1->GetTwoTrackMomentum().M();
        cutting.muee_P   = dir1->GetThreeTrackMomentum().P();
        cutting.muee_Pt  = dir1->GetThreeTrackMomentum().Pt();
        cutting.muee_M   = dir1->GetThreeTrackMomentum().M();
        cutting.zvtx_mue1_mue2 = Vertex_mu_e1[2] - Vertex_mu_e2[2];
        cutting.zvtx_mue1_e1e2 = Vertex_mu_e1[2] - Vertex_e1_e2[2];
        cutting.zvtx_mue2_e1e2 = Vertex_mu_e2[2] - Vertex_e1_e2[2];
        cutting.yvtx_mue1_mue2 = Vertex_mu_e1[1] - Vertex_mu_e2[1];
        cutting.yvtx_mue1_e1e2 = Vertex_mu_e1[1] - Vertex_e1_e2[1];
        cutting.yvtx_mue2_e1e2 = Vertex_mu_e2[1] - Vertex_e1_e2[1];
        cutting.xvtx_mue1_mue2 = Vertex_mu_e1[0] - Vertex_mu_e2[0];
        cutting.xvtx_mue1_e1e2 = Vertex_mu_e1[0] - Vertex_e1_e2[0];
        cutting.xvtx_mue2_e1e2 = Vertex_mu_e2[0] - Vertex_e1_e2[0];
        cutting.DCH_Radius_mu  = muon.GetDCHradius();
        cutting.DCH_Radius_el1 = electron1.GetDCHradius();
        cutting.DCH_Radius_el2 = electron2.GetDCHradius();
    }
    if(decay_type.compare("K3pi") == 0 ){
        cutting.muee_P   = dir1->GetThreeTrackMomentum().P();
        cutting.muee_Pt  = dir1->GetThreeTrackMomentum().Pt();
        cutting.muee_M   = dir1->GetThreeTrackMomentum().M();
        cutting.zvtx_pi1pi2_pi2pi3= Vertex_mu_e1[2] - Vertex_mu_e2[2];
        cutting.zvtx_pi1pi2_pi1pi3= Vertex_mu_e1[2] - Vertex_e1_e2[2];
        cutting.zvtx_pi1pi3_pi2pi3= Vertex_mu_e2[2] - Vertex_e1_e2[2];
        cutting.yvtx_pi1pi2_pi2pi3= Vertex_mu_e1[1] - Vertex_mu_e2[1];
        cutting.yvtx_pi1pi2_pi1pi3= Vertex_mu_e1[1] - Vertex_e1_e2[1];
        cutting.yvtx_pi1pi3_pi2pi3= Vertex_mu_e2[1] - Vertex_e1_e2[1];
        cutting.xvtx_pi1pi2_pi2pi3= Vertex_mu_e1[0] - Vertex_mu_e2[0];
        cutting.xvtx_pi1pi2_pi1pi3= Vertex_mu_e1[0] - Vertex_e1_e2[0];
        cutting.xvtx_pi1pi3_pi2pi3= Vertex_mu_e2[0] - Vertex_e1_e2[0];
        cutting.DCH_Radius_pi1 = muon.GetDCHradius();
        cutting.DCH_Radius_pi2 = electron1.GetDCHradius();
        cutting.DCH_Radius_pi3 = electron2.GetDCHradius();
    }

    return cutting;
}
