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
//float __lda3(superCmpEvent*, superCmpEvent*, int);

void cleanup();
Cuts make_cuts(Hist_dir* dir1, Charged_Particle& muon, Charged_Particle& electron1, Charged_Particle& electron2,double* Vertex_mu_e1,double* Vertex_mu_e2,double* Vertex_e1_e2, std::string decay_type);

void FillMC(Hist_dir* dir, TLorentzVector p1, TLorentzVector p2, TLorentzVector p3, float* DKaon,float* part_production, float* part_decay);

void FillMC(Hist_dir* dir, TLorentzVector p1, TLorentzVector p2, TLorentzVector p3, TLorentzVector p4, float* DKaon,float* part_production, float* part_decay);

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
    double chi2_Z_Vertex = sevt->vtx[0].chi2;
    static double massKaonC = 0.493677;
    const double Electron_EoverP_up = 1.05;
    const double Electron_EoverP_down = 0.95;
    const double Muon_EoverP_up = 0.2;
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

    ////////////MC weights EoP*lda3//////////////////
    /* double w_0_5  = 0.940576*0.960979; */
    /* double w_5_10 = 0.962015*0.971205; */
    /* double w_10_15= 0.972047*0.974435; */
    /* double w_15_20= 0.975931*0.975205; */
    /* double w_20_25= 0.982185*0.973952; */
    /* double w_25_30= 0.986453*0.960957; */
    /* double w_30_35= 0.98318 *0.955836; */
    /* double w_35_40= 0.984925*0.953608; */
    double wmu_5_55  = 0.907;
    double wmu_55_6  = 0.94;
    double wmu_6_65  = 0.957;
    double wmu_65_7  = 0.963;
    double wmu_7_75  = 0.972;
    double wmu_75_8 = 0.98;
    double wmu_8_9 = 0.987;
    double wmu_9_10 = 0.99;
    double w_0_5  = 0.938829*0.959385;
    double w_5_10 = 0.960154*0.970754;
    double w_10_15= 0.973072*0.976971;
    double w_15_20= 0.978127*0.976439;
    double w_20_25= 0.980435*0.964764;
    double w_25_30= 0.980477*0.964764;
    double w_30_35= 0.978684*0.954012;
    double w_35_40= 0.977009*0.940663;

    double w_kaon= 1.0;
    double w_mu  = 1.0;
    double w_el1 = 1.0;
    double w_el2 = 1.0;
    double w_tot = 1.0;
    ////////////////////////////////////////

    double Track_Momentum;
    double Track_Charge;
    double Track_Quality;
    double Track_DeadCell_Distance;
    double Track_Energy;
    double Track_EoverP;
    double MuTrack_EoverP;

    /// apply non linearity corrections
    /// Correct LKr non-linearity (effective for energies < 11GeV)
    if (LKrCalCorr && IS_DATA) {
        user_lkrcalcor_SC (sbur,sevt,1);
    }

    Initial_dir->fh_Nrun->Fill(sbur->nrun);
    Initial_dir->FillCommonHist(sevt);

    //if(IS_MC){
    //    if( sbur->nrun > 16000)return 0;
    //}

    if(COmPaCt_Z_Vertex < -1800. || COmPaCt_Z_Vertex > 8000.){return 0;}

    if(IS_MC){
        FillMC(Initial_dir, True_Momentum[1], True_Momentum[2], True_Momentum[3], DKaon, Particle_production_zvtx, Particle_decay_zvtx);
        if(Npart >= 4){
            FillMC(Initial_dir, True_Momentum[1], True_Momentum[2], True_Momentum[3], True_Momentum[4], DKaon, Particle_production_zvtx, Particle_decay_zvtx);
        }
    }

    if(chi2_Z_Vertex > 20.){return 0;}

    if(nvtx != 1){return 0;}
    if(ntrack != 3 ) {return 0;}


    if(IS_DATA){
        if(sbur->nrun > 16000){
            if(sevt->LKRdownscaled) return 0;
        }
        if(sbur->nrun >= 16586 && sbur->nrun <= 16606){
            int dn = 5;
            int nev = sevt->nEvt;
            //double rejected;
            //double accepted;
            //double total;

            if( ((sevt->trigWord)>>4)&0x1){
                //cout << nev%16 << endl;
                //cout << "Bit 4 ==== " <<  (((sevt->trigWord)>>4)&0x1) << "Total word == " << sevt->trigWord << endl;
                if( (int)nev%16 < (int)dn){
                    //rejected++;
                    return 0;
                }//else {
                //accepted++;
                //}
                //total++;
            }

            //if(total==16){
            //
            //    //cout << "Rejected == " << rejected << " Accepted == " << accepted << " Total == " << total << endl;
            //    total=0;
            //    rejected=0;
            //    accepted=0;
            //}


            if( ((sevt->trigWord)>>1)&0x1 && sevt->trigWord > 2){

                if( (int)nev%16 < (int)dn){
                    return 0;
                }
            }
            if( ((sevt->trigWord))&0x1 && sevt->trigWord > 1){

                if( (int)nev%16 < (int)dn){
                    return 0;
                }
                //cout << "Bit 0 ==== " <<  (((sevt->trigWord))&0x1) << "Total word == " << sevt->trigWord << endl;
            }
        }

        if(sbur->nrun >= 16607){
            int dn = 6;
            int nev = sevt->nEvt;

            if( ((sevt->trigWord)>>4)&0x1){

                if( (int)nev%16 < (int)dn)return 0;

            }

            if( ((sevt->trigWord)>>1)&0x1 && sevt->trigWord > 3){

                if( (int)nev%16 < (int)dn) return 0;
            }

            if( ((sevt->trigWord))&0x1 && sevt->trigWord > 1){

                if( (int)nev%16 < (int)dn) return 0;
            }
        }
    }

    for (int j=0; j<ntrack;j++){
        Kcharge += sevt->track[j].q;
    }

    //cout << IS_DATA << "    " << IS_MC << endl;
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
        if (Track_DeadCell_Distance < 2.){continue;}

        ngoodtrack++;

        Initial_dir->fh_eop->Fill(Track_EoverP);
        Initial_dir->fh_Track_Momentum->Fill(Track_Momentum);

        //Particle identification for K to mu nu e e
        if(IS_DATA){
            if(Track_EoverP > Electron_EoverP_down && Track_EoverP < Electron_EoverP_up){
                if(iel1<0){
                    iel1 = i;
                    /* lda3_e1 = lda3(sevt,sevt,i+1); */
                    lda3_e1 = ldaC(sevt,sbur,i);
                } else if(iel2<0){
                    iel2 = i;
                    /* lda3_e2 = lda3(sevt,sevt,i+1); */
                    lda3_e2 = ldaC(sevt,sbur,i);

                } else {;}
                Nelectrons++;
            } else if (imu < 0 && Track_imu!= -1 && Track_EoverP < Muon_EoverP_up ){
                //cout << sevt->Nmuon << endl;
                //cout << sevt->track[i].iMuon <<  "Nmuon == " << sevt->Nmuon  << endl;
                imu = i;
                Nmuons++;
            }


        }// ENDIF IS_DATA

        //Particle identification for K to 3 pi charged
        if( Track_Momentum > 10. ){
            if(pi1<0){
                pi1 = i;
                if(IS_DATA)
                    lda3_pi1 = ldaC(sevt,sbur,pi1);
            } else if(pi2<0){
                pi2 = i;
                if(IS_DATA)
                    lda3_pi2 = ldaC(sevt,sbur,pi2);
            } else if(pi3<0){
                pi3 = i;
                if(IS_DATA)
                    lda3_pi3 = ldaC(sevt,sbur,pi3);
            }

            NK3pi_pions++;
        }//ENDIF K3pi selection


        if(IS_MC){
            if(Kcharge*Track_Charge == -1.){

                if( iel1 < 0 && Track_EoverP > 0.3){
                    iel1 = i;

                }

            } else if ( iel2 < 0 || imu < 0 ){

                if(imu < 0 && Track_EoverP < 0.3){
                    imu = i;
                    //cout << "iel2 = " << iel2 << " ipi =" << ipi << endl;
                    //cout << "Kcharge = " << Kcharge << "Track_charge = " << Track_Charge << " ipi =" << ipi << endl;
                    MuTrack_EoverP= sevt->cluster[sevt->track[imu].iClus].energy/sevt->track[imu].p;
                    //cout << " Pion E/p =  " << PiTrack_EoverP << "Electron E/p =  " << Track_EoverP << endl;
                    //cout << "NEXT _________+_+_+)_+)_+)_+)_+)_++__++" << endl;
                } else if(imu > 0 && Track_EoverP > MuTrack_EoverP && Track_EoverP > 0.3){
                    iel2 = i;
                    //lda3_e2 = ldaC(sevt,sbur,i);
                    //cout << " iel2 =  " << iel2 << "ipi =  " <<  ipi << endl;
                    //cout << " Pion E/p =  " << PiTrack_EoverP << "Electron E/p =  " << Track_EoverP << endl;
                }else {
                    if(Track_EoverP > 0.3){
                        iel2 = i;
                        lda3_e2 = ldaC(sevt,sbur,i);
                        //     ipi  = i;
                    }

                }
            }
        }
        //Ke4 MC PID
        //if (IS_MC){
        //    if(Track_EoverP > Electron_EoverP_down && Track_EoverP < Electron_EoverP_up && sevt->track[i].q ==1)    {
        //        //if(eop_track> 0.8)    {
        //        iel1=i;
        //    }
        //    else if(iel2 < 0 && Kcharge*sevt->track[i].q == -1 )    {
        //        iel2=i;
        //    }
        //    else if( imu < 0) {
        //        imu=i;
        //    }
        //    else    {
        //        ;
        //    }
        //}

    }//end for i
    //Kmunuee selection
    Charged_Particle muon(sevt,sbur,-13,imu);
    Charged_Particle electron1(sevt,sbur,-11,iel1);
    Charged_Particle electron2(sevt,sbur,-11,iel2);

    //K3pi wrong sign selection
    Charged_Particle k3pi_pion1(sevt,sbur,211,imu);
    Charged_Particle k3pi_pion2(sevt,sbur,211,iel1);
    Charged_Particle k3pi_pion3(sevt,sbur,211,iel2);




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

    //Checking for goodness of tracks
    if(ngoodtrack!= 3){return 0;}

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
        closap_double_(pion1.Position,pion2.Position,pion1.Slopes,pion2.Slopes,&cda_pi1_pi2,Vertex_pi1_pi2);
        closap_double_(pion1.Position,pion3.Position,pion1.Slopes,pion3.Slopes,&cda_pi1_pi3,Vertex_pi1_pi3);
        closap_double_(pion2.Position,pion3.Position,pion2.Slopes,pion3.Slopes,&cda_pi2_pi3,Vertex_pi2_pi3);
        K3pi_selection->ComputeThreeTrack(pion1,pion2,pion3);
        Cuts cut_k3pi = make_cuts(K3pi_selection,pion1,pion2,pion3,Vertex_pi1_pi2,Vertex_pi1_pi3,Vertex_pi2_pi3,"K3pi");
        ////Cuts
        double K3pi_Event_HoDTime;
        double K3pi_Event_DCHtime;
        K3pi_Event_HoDTime = ( (pion1.GetHodTime() - pion2.GetHodTime()) + (pion2.GetHodTime() - pion3.GetHodTime() ) + (pion1.GetHodTime() - pion3.GetHodTime() ))/ 3.;
        K3pi_Event_DCHtime = ( (pion1.GetDCHtime() - pion2.GetDCHtime()) + (pion2.GetDCHtime() - pion3.GetDCHtime() ) + (pion1.GetDCHtime() - pion3.GetDCHtime() ))/ 3.;
        //-- CUT1 Momentum cut ---
        if(cut_k3pi.muee_P > 54. &&
           cut_k3pi.muee_P < 66. &&

           //--ENDOF CUT1 Momentum cut ---
           //-- CUT2 Timing cut --- (DATA only)
           fabs(cut_k3pi.DCH_e1e2) < 10.&&
           fabs(cut_k3pi.DCH_mue1) < 10.&&
           fabs(cut_k3pi.DCH_mue2) < 10.&&
           fabs(cut_k3pi.Hod_e1e2) < 3. &&
           fabs(cut_k3pi.Hod_mue1) < 3. &&
           fabs(cut_k3pi.Hod_mue2) < 3. &&
           //-- ENDOF CUT2 Timing cut ---
           //fabs(K3pi_Event_HoDTime) < 3. &&
           //fabs(K3pi_Event_DCHtime) < 10. &&
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
           //cut_k3pi.muee_M > 0.4886 &&
           //cut_k3pi.muee_M < 0.4985 &&
           //-- ENDOF CUT4 Three track invariant mass Cut --

           //-- CUT5 DCH geometry Cut --
           cut_k3pi.DCH_Radius_pi1 > 14. &&
           cut_k3pi.DCH_Radius_pi1 < 110.&&
           cut_k3pi.DCH_Radius_pi2 > 14. &&
           cut_k3pi.DCH_Radius_pi2 < 110.&&
           cut_k3pi.DCH_Radius_pi3 > 14. &&
           cut_k3pi.DCH_Radius_pi3 < 110.&&

           //Lkr octagonal cut
           //          sqrt(x2+y2) > 15cm
           //|x| < 113 cm
           //|y| < 113 cm
           //|x| + |y| < 159.8 cm
           cut_k3pi.Lkr_cut_pi1  > 15. &&
           cut_k3pi.Lkr_cut_pi2  > 15. &&
           cut_k3pi.Lkr_cut_pi3  > 15. &&
           fabs(cut_k3pi.Lkr_x_pi1) < 113. &&
           fabs(cut_k3pi.Lkr_x_pi2) < 113. &&
           fabs(cut_k3pi.Lkr_x_pi3) < 113. &&
           fabs(cut_k3pi.Lkr_y_pi1) < 113. &&
           fabs(cut_k3pi.Lkr_y_pi2) < 113. &&
           fabs(cut_k3pi.Lkr_y_pi3) < 113. &&
           (fabs(cut_k3pi.Lkr_x_pi3) + fabs(cut_k3pi.Lkr_y_pi3) ) < 159.8 &&
           (fabs(cut_k3pi.Lkr_x_pi2) + fabs(cut_k3pi.Lkr_y_pi2) ) < 159.8 &&
           (fabs(cut_k3pi.Lkr_x_pi1) + fabs(cut_k3pi.Lkr_y_pi1) ) < 159.8 &&
           COmPaCt_Z_Vertex > -1800 &&
           COmPaCt_Z_Vertex < 8000.
           //-- ENDOF CUT5 DCH geometry Cut --
           ){

            double trigger = 1;
            // ------------------------------------------ SS0 --------------------------------------//
            /* if( 15304 < sbur->nrun < 15582){ */
            if( 15304 < sbur->nrun && sbur->nrun <= 15582){
                if(((sevt->trigWord)>>3)&0x1) {
                    // ref. trigg. ok:
                    K3pi_selection->fh_SS0_CPRE->Fill(trigger);
                    if( (sevt->pu[3].chan[5]>>4)&0x1 ||
                        (sevt->pu[4].chan[5]>>4)&0x1 ||
                        (sevt->pu[5].chan[5]>>4)&0x1 ||  //5 is the central time slice
                        (sevt->pu[6].chan[5]>>4)&0x1 ||
                        (sevt->pu[3].chan[5]>>8)&0x1 ||
                        (sevt->pu[4].chan[5]>>8)&0x1 ||
                        (sevt->pu[5].chan[5]>>8)&0x1 ||  //5 is the central time slice
                        (sevt->pu[6].chan[5]>>8)&0x1 ||
                        (sevt->pu[3].chan[5]>>13)&0x1 ||
                        (sevt->pu[4].chan[5]>>13)&0x1 ||
                        (sevt->pu[5].chan[5]>>13)&0x1 ||  //5 is the central time slice
                        (sevt->pu[6].chan[5]>>13)&0x1 ) {
                        if( ((sevt->trigWord)>>4)&0x1 ) //
                            K3pi_selection->fh_SS0_MB_1TRK_P->Fill(trigger);
                        if( ((sevt->trigWord)>>1)&0x1  )
                            K3pi_selection->fh_SS0_MB_1VTX->Fill(trigger);

                        if( ((sevt->trigWord))&0x1    || ((sevt->trigWord)>>1)&0x1 || ((sevt->trigWord)>>4)&0x1 )
                            K3pi_selection->fh_SS0_full_trig->Fill(trigger);

                    }
                }
            }
            //--------------------------------------- SS1 ---------------------------------------//
            if( 15633 <= sbur->nrun && sbur->nrun <= 15703){
                if(((sevt->trigWord)>>3)&0x1) {
                    // ref. trigg. ok:
                    K3pi_selection->fh_SS1_CPRE->Fill(trigger);
                    if( (sevt->pu[3].chan[5]>>4)&0x1 ||
                        (sevt->pu[4].chan[5]>>4)&0x1 ||
                        (sevt->pu[5].chan[5]>>4)&0x1 ||  //5 is the central time slice
                        (sevt->pu[6].chan[5]>>4)&0x1 ||
                        (sevt->pu[3].chan[5]>>8)&0x1 ||
                        (sevt->pu[4].chan[5]>>8)&0x1 ||
                        (sevt->pu[5].chan[5]>>8)&0x1 ||  //5 is the central time slice
                        (sevt->pu[6].chan[5]>>8)&0x1 ||
                        (sevt->pu[3].chan[5]>>13)&0x1 ||
                        (sevt->pu[4].chan[5]>>13)&0x1 ||
                        (sevt->pu[5].chan[5]>>13)&0x1 ||  //5 is the central time slice
                        (sevt->pu[6].chan[5]>>13)&0x1 ) {
                        if( ((sevt->trigWord)>>4)&0x1 ) //
                            K3pi_selection->fh_SS1_MB_1TRK_P->Fill(trigger);
                        if( ((sevt->trigWord)>>1)&0x1  )
                            K3pi_selection->fh_SS1_MB_1VTX->Fill(trigger);
                        if( ((sevt->trigWord))&0x1    || ((sevt->trigWord)>>1)&0x1 || ((sevt->trigWord)>>4)&0x1 )
                            K3pi_selection->fh_SS1_full_trig->Fill(trigger);
                    }
                }
            }
            //------------------------------------- SS2 -----------------------------------------------//
            if( 15717 <= sbur->nrun && sbur->nrun <= 15777){
                if(((sevt->trigWord)>>3)&0x1) {
                    // ref. trigg. ok:
                    K3pi_selection->fh_SS2_CPRE->Fill(trigger);
                    if( (sevt->pu[3].chan[5]>>4)&0x1 ||
                        (sevt->pu[4].chan[5]>>4)&0x1 ||
                        (sevt->pu[5].chan[5]>>4)&0x1 ||  //5 is the central time slice
                        (sevt->pu[6].chan[5]>>4)&0x1 ||
                        (sevt->pu[3].chan[5]>>8)&0x1 ||
                        (sevt->pu[4].chan[5]>>8)&0x1 ||
                        (sevt->pu[5].chan[5]>>8)&0x1 ||  //5 is the central time slice
                        (sevt->pu[6].chan[5]>>8)&0x1 ||
                        (sevt->pu[3].chan[5]>>13)&0x1 ||
                        (sevt->pu[4].chan[5]>>13)&0x1 ||
                        (sevt->pu[5].chan[5]>>13)&0x1 ||  //5 is the central time slice
                        (sevt->pu[6].chan[5]>>13)&0x1 ) {
                        if( ((sevt->trigWord)>>4)&0x1 ) //
                            K3pi_selection->fh_SS2_MB_1TRK_P->Fill(trigger);
                        if( ((sevt->trigWord)>>1)&0x1  )
                            K3pi_selection->fh_SS2_MB_1VTX->Fill(trigger);
                        if( ((sevt->trigWord))&0x1    || ((sevt->trigWord)>>1)&0x1 || ((sevt->trigWord)>>4)&0x1 )
                            K3pi_selection->fh_SS2_full_trig->Fill(trigger);
                    }
                }
            }
            //------------------------------------- SS3 -----------------------------------------------//
            if( 15778 <= sbur->nrun && sbur->nrun <= 15790){
                if(((sevt->trigWord)>>3)&0x1) {
                    // ref. trigg. ok:
                    K3pi_selection->fh_SS3_CPRE->Fill(trigger);
                    if( (sevt->pu[3].chan[5]>>4)&0x1 ||
                        (sevt->pu[4].chan[5]>>4)&0x1 ||
                        (sevt->pu[5].chan[5]>>4)&0x1 ||  //5 is the central time slice
                        (sevt->pu[6].chan[5]>>4)&0x1 ||
                        (sevt->pu[3].chan[5]>>8)&0x1 ||
                        (sevt->pu[4].chan[5]>>8)&0x1 ||
                        (sevt->pu[5].chan[5]>>8)&0x1 ||  //5 is the central time slice
                        (sevt->pu[6].chan[5]>>8)&0x1 ||
                        (sevt->pu[3].chan[5]>>13)&0x1 ||
                        (sevt->pu[4].chan[5]>>13)&0x1 ||
                        (sevt->pu[5].chan[5]>>13)&0x1 ||  //5 is the central time slice
                        (sevt->pu[6].chan[5]>>13)&0x1 ) {
                        if( ((sevt->trigWord)>>4)&0x1 ) //
                            K3pi_selection->fh_SS3_MB_1TRK_P->Fill(trigger);
                        if( ((sevt->trigWord)>>1)&0x1  )
                            K3pi_selection->fh_SS3_MB_1VTX->Fill(trigger);
                        if( ((sevt->trigWord))&0x1    || ((sevt->trigWord)>>1)&0x1 || ((sevt->trigWord)>>4)&0x1 )
                            K3pi_selection->fh_SS3_full_trig->Fill(trigger);
                    }
                }
            }

            //------------------------------------- SS4 -----------------------------------------------//
            if( 16121 <= sbur->nrun && sbur->nrun <= 16383){
                if(((sevt->trigWord)>>3)&0x1) {
                    // ref. trigg. ok:
                    K3pi_selection->fh_SS4_CPRE->Fill(trigger);
                    if( (sevt->pu[3].chan[5]>>4)&0x1 ||
                        (sevt->pu[4].chan[5]>>4)&0x1 ||
                        (sevt->pu[5].chan[5]>>4)&0x1 ||  //5 is the central time slice
                        (sevt->pu[6].chan[5]>>4)&0x1 ||
                        (sevt->pu[3].chan[5]>>8)&0x1 ||
                        (sevt->pu[4].chan[5]>>8)&0x1 ||
                        (sevt->pu[5].chan[5]>>8)&0x1 ||  //5 is the central time slice
                        (sevt->pu[6].chan[5]>>8)&0x1 ||
                        (sevt->pu[3].chan[5]>>13)&0x1 ||
                        (sevt->pu[4].chan[5]>>13)&0x1 ||
                        (sevt->pu[5].chan[5]>>13)&0x1 ||  //5 is the central time slice
                        (sevt->pu[6].chan[5]>>13)&0x1 ) {
                        if( ((sevt->trigWord)>>4)&0x1 ) //
                            K3pi_selection->fh_SS4_MB_1TRK_P->Fill(trigger);
                        if( ((sevt->trigWord)>>1)&0x1  )
                            K3pi_selection->fh_SS4_MB_1VTX->Fill(trigger);
                        if( ((sevt->trigWord))&0x1    || ((sevt->trigWord)>>1)&0x1 || ((sevt->trigWord)>>4)&0x1 )
                            K3pi_selection->fh_SS4_full_trig->Fill(trigger);
                    }
                }
            }


            //------------------------------------- SS5 -----------------------------------------------//
            if( 16428 <= sbur->nrun && sbur->nrun <= 16585){
                if(((sevt->trigWord)>>3)&0x1) {
                    // ref. trigg. ok:
                    K3pi_selection->fh_SS5_CPRE->Fill(trigger);
                    if( (sevt->pu[3].chan[5]>>4)&0x1 ||
                        (sevt->pu[4].chan[5]>>4)&0x1 ||
                        (sevt->pu[5].chan[5]>>4)&0x1 ||  //5 is the central time slice
                        (sevt->pu[6].chan[5]>>4)&0x1 ||
                        (sevt->pu[3].chan[5]>>8)&0x1 ||
                        (sevt->pu[4].chan[5]>>8)&0x1 ||
                        (sevt->pu[5].chan[5]>>8)&0x1 ||  //5 is the central time slice
                        (sevt->pu[6].chan[5]>>8)&0x1 ||
                        (sevt->pu[3].chan[5]>>13)&0x1 ||
                        (sevt->pu[4].chan[5]>>13)&0x1 ||
                        (sevt->pu[5].chan[5]>>13)&0x1 ||  //5 is the central time slice
                        (sevt->pu[6].chan[5]>>13)&0x1 ) {
                        if( ((sevt->trigWord)>>4)&0x1 ) //
                            K3pi_selection->fh_SS5_MB_1TRK_P->Fill(trigger);
                        if( ((sevt->trigWord)>>1)&0x1  )
                            K3pi_selection->fh_SS5_MB_1VTX->Fill(trigger);
                        if( ((sevt->trigWord))&0x1    || ((sevt->trigWord)>>1)&0x1 || ((sevt->trigWord)>>4)&0x1 )
                            K3pi_selection->fh_SS5_full_trig->Fill(trigger);
                    }
                }
            }

            //------------------------------------- SS6 -----------------------------------------------//
            if( 16586 <= sbur->nrun && sbur->nrun <= 16709){
                if(((sevt->trigWord)>>3)&0x1) {
                    // ref. trigg. ok:
                    K3pi_selection->fh_SS6_CPRE->Fill(trigger);
                    if( (sevt->pu[3].chan[5]>>4)&0x1 ||
                        (sevt->pu[4].chan[5]>>4)&0x1 ||
                        (sevt->pu[5].chan[5]>>4)&0x1 ||  //5 is the central time slice
                        (sevt->pu[6].chan[5]>>4)&0x1 ||
                        (sevt->pu[3].chan[5]>>8)&0x1 ||
                        (sevt->pu[4].chan[5]>>8)&0x1 ||
                        (sevt->pu[5].chan[5]>>8)&0x1 ||  //5 is the central time slice
                        (sevt->pu[6].chan[5]>>8)&0x1 ||
                        (sevt->pu[3].chan[5]>>13)&0x1 ||
                        (sevt->pu[4].chan[5]>>13)&0x1 ||
                        (sevt->pu[5].chan[5]>>13)&0x1 ||  //5 is the central time slice
                        (sevt->pu[6].chan[5]>>13)&0x1 ) {
                        if( ((sevt->trigWord)>>4)&0x1 ) //
                            K3pi_selection->fh_SS6_MB_1TRK_P->Fill(trigger);
                        if( ((sevt->trigWord)>>1)&0x1  )
                            K3pi_selection->fh_SS6_MB_1VTX->Fill(trigger);
                        if( ((sevt->trigWord))&0x1    || ((sevt->trigWord)>>1)&0x1 || ((sevt->trigWord)>>4)&0x1 )
                            K3pi_selection->fh_SS6_full_trig->Fill(trigger);
                    }
                }
            }

            //------------------------------------- SS7 -----------------------------------------------//
            if( 16722 <= sbur->nrun && sbur->nrun <= 16801){
                if(((sevt->trigWord)>>3)&0x1) {
                    // ref. trigg. ok:
                    K3pi_selection->fh_SS7_CPRE->Fill(trigger);
                    if( (sevt->pu[3].chan[5]>>4)&0x1 ||
                        (sevt->pu[4].chan[5]>>4)&0x1 ||
                        (sevt->pu[5].chan[5]>>4)&0x1 ||  //5 is the central time slice
                        (sevt->pu[6].chan[5]>>4)&0x1 ||
                        (sevt->pu[3].chan[5]>>8)&0x1 ||
                        (sevt->pu[4].chan[5]>>8)&0x1 ||
                        (sevt->pu[5].chan[5]>>8)&0x1 ||  //5 is the central time slice
                        (sevt->pu[6].chan[5]>>8)&0x1 ||
                        (sevt->pu[3].chan[5]>>13)&0x1 ||
                        (sevt->pu[4].chan[5]>>13)&0x1 ||
                        (sevt->pu[5].chan[5]>>13)&0x1 ||  //5 is the central time slice
                        (sevt->pu[6].chan[5]>>13)&0x1 ) {
                        if( ((sevt->trigWord)>>4)&0x1 ) //
                            K3pi_selection->fh_SS7_MB_1TRK_P->Fill(trigger);
                        if( ((sevt->trigWord)>>1)&0x1  )
                            K3pi_selection->fh_SS7_MB_1VTX->Fill(trigger);
                        if( ((sevt->trigWord))&0x1    || ((sevt->trigWord)>>1)&0x1 || ((sevt->trigWord)>>4)&0x1 )
                            K3pi_selection->fh_SS7_full_trig->Fill(trigger);
                    }
                }
            }

            //------------------------------------- SS8 -----------------------------------------------//
            if( 16802 <= sbur->nrun && sbur->nrun <= 16905){
                if(((sevt->trigWord)>>3)&0x1) {
                    // ref. trigg. ok:
                    K3pi_selection->fh_SS8_CPRE->Fill(trigger);
                    if( (sevt->pu[3].chan[5]>>4)&0x1 ||
                        (sevt->pu[4].chan[5]>>4)&0x1 ||
                        (sevt->pu[5].chan[5]>>4)&0x1 ||  //5 is the central time slice
                        (sevt->pu[6].chan[5]>>4)&0x1 ||
                        (sevt->pu[3].chan[5]>>8)&0x1 ||
                        (sevt->pu[4].chan[5]>>8)&0x1 ||
                        (sevt->pu[5].chan[5]>>8)&0x1 ||  //5 is the central time slice
                        (sevt->pu[6].chan[5]>>8)&0x1 ||
                        (sevt->pu[3].chan[5]>>13)&0x1 ||
                        (sevt->pu[4].chan[5]>>13)&0x1 ||
                        (sevt->pu[5].chan[5]>>13)&0x1 ||  //5 is the central time slice
                        (sevt->pu[6].chan[5]>>13)&0x1 ) {
                        if( ((sevt->trigWord)>>4)&0x1 ) //
                            K3pi_selection->fh_SS8_MB_1TRK_P->Fill(trigger);
                        if( ((sevt->trigWord)>>1)&0x1  )
                            K3pi_selection->fh_SS8_MB_1VTX->Fill(trigger);
                        if( ((sevt->trigWord))&0x1    || ((sevt->trigWord)>>1)&0x1 || ((sevt->trigWord)>>4)&0x1 )
                            K3pi_selection->fh_SS8_full_trig->Fill(trigger);
                    }
                }
            }

            K3pi_selection->FillCommonHist(sevt);

            K3pi_selection->FillHist(pion1,"pion1");
            K3pi_selection->FillHist(pion2,"pion2");
            K3pi_selection->FillHist(pion3,"pion3");

            K3pi_selection->FillHist(K3pi_selection->GetThreeTrackMomentum(),K3pi_selection->GetNuMomentum(), Kcharge);
            K3pi_selection->fh_Kaon_Charge->Fill(Kcharge);
            K3pi_selection->fh_Event_Type->Fill(K3pi_Event_Type);
            K3pi_selection->fh_HoDTotTime->Fill( K3pi_Event_HoDTime );
            K3pi_selection->fh_DCHTotTime->Fill( K3pi_Event_DCHtime );
            if(IS_DATA){
                if(pion1.GetEnergyLeftInEcal()/ pion1.GetMomentum() > 0.95 && pion1.GetEnergyLeftInEcal()/ pion1.GetMomentum() < 1.05 ) {

                    K3pi_selection->fh_lda3_p1->Fill(lda3_pi1);
                    //K3pi_selection->fh_lda3_p1->Fill(1.);
                    //cout << lda3_pi1 << endl;
                } else if(pion2.GetEnergyLeftInEcal()/ pion2.GetMomentum() > 0.95 && pion2.GetEnergyLeftInEcal()/ pion2.GetMomentum() < 1.05 ) {
                    K3pi_selection->fh_lda3_p2->Fill(lda3_pi2);
                } else if(pion3.GetEnergyLeftInEcal()/ pion3.GetMomentum() > 0.95 && pion3.GetEnergyLeftInEcal()/ pion3.GetMomentum() < 1.05 ) {
                    K3pi_selection->fh_lda3_p3->Fill(lda3_pi3);
                }
            }
            //K3pi_selection->fh_lda3_p3->Fill(lda3_pi3);

            K3pi_selection->FillVertexHist(Vertex_pi1_pi2, cda_pi1_pi2 , Vertex_pi2_pi3, cda_pi2_pi3, Vertex_pi1_pi3, cda_pi1_pi3,"K3pi");

            if(Kcharge*pion1.GetCharge()==-1 && sevt->track[pi1].iMuon != -1){
                K3pi_selection->fh_odd_eop->Fill(pion1.GetEnergyLeftInEcal()/ pion1.GetMomentum());
                K3pi_selection->fh_EoP_vs_p_odd_tr->Fill(pion1.GetEnergyLeftInEcal()/ pion1.GetMomentum(), pion1.GetMomentum());
            } else if(Kcharge*pion2.GetCharge()==-1 && sevt->track[pi2].iMuon != -1) {
                K3pi_selection->fh_odd_eop->Fill(pion2.GetEnergyLeftInEcal()/ pion2.GetMomentum());
                K3pi_selection->fh_EoP_vs_p_odd_tr->Fill(pion2.GetEnergyLeftInEcal()/ pion2.GetMomentum(), pion2.GetMomentum());
            } else if(Kcharge*pion3.GetCharge()==-1 && sevt->track[pi3].iMuon != -1){
                K3pi_selection->fh_odd_eop->Fill(pion3.GetEnergyLeftInEcal()/ pion3.GetMomentum());
                K3pi_selection->fh_EoP_vs_p_odd_tr->Fill(pion3.GetEnergyLeftInEcal()/ pion3.GetMomentum(), pion3.GetMomentum());

            }
        }
    }
    //-------------------------KMUNUEE SIGNAL SELECTION ----------------------------------------
    //Signal particle identification
    if(imu == -1 || iel1 == -1 || iel2 == -1){return 0;}
    //if(IS_DATA)
    if(!electron1.cluster_exists || !electron2.cluster_exists){return 0;}

    double Event_HoDTime;
    double Event_MuHoDTime;
    double Event_E1HoDTime;
    double Event_E2HoDTime;
    double Event_DCHTime;

    /* Event_HoDTime = ( (muon.GetHodTime() - electron1.GetHodTime()) + (electron1.GetHodTime() - electron2.GetHodTime() ) + (muon.GetHodTime() - electron2.GetHodTime() ))/ 3.; */
    Event_MuHoDTime = ( (muon.GetHodTime() + electron1.GetHodTime() + electron2.GetHodTime() )/3 -  muon.GetHodTime());
    Event_E1HoDTime = ( (muon.GetHodTime() + electron1.GetHodTime() + electron2.GetHodTime() )/3 -  electron1.GetHodTime());
    Event_E2HoDTime = ( (muon.GetHodTime() + electron1.GetHodTime() + electron2.GetHodTime() )/3 -  electron2.GetHodTime());

    if(fabs(Event_MuHoDTime) > fabs(Event_E1HoDTime) && fabs(Event_MuHoDTime) > fabs(Event_E2HoDTime) )
        Event_HoDTime = Event_MuHoDTime ;

    if(fabs(Event_E1HoDTime) > fabs(Event_MuHoDTime) && fabs(Event_MuHoDTime) > fabs(Event_E2HoDTime) )
        Event_HoDTime = Event_E1HoDTime ;

    if( fabs(Event_E2HoDTime) > fabs(Event_MuHoDTime) && fabs(Event_MuHoDTime) > fabs(Event_E1HoDTime) )
        Event_HoDTime = Event_E2HoDTime ;

    Event_DCHTime = ( muon.GetDCHtime() - (electron1.GetDCHtime() + electron2.GetDCHtime() )/3 );

    dir1->FillCommonHist(sevt);
    dir1->fh_Event_Type->Fill(Event_Type);
    dir1->fh_Kaon_Charge->Fill(Kcharge);
    //Fill MC histograms (if IS_MC == 1)
    if(IS_MC){
        FillMC(dir1, True_Momentum[1], True_Momentum[2], True_Momentum[3], DKaon, Particle_production_zvtx, Particle_decay_zvtx);
        if(Npart >= 4){
            FillMC(dir1, True_Momentum[1], True_Momentum[2], True_Momentum[3], True_Momentum[4], DKaon, Particle_production_zvtx, Particle_decay_zvtx);

        }
    }

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
    dir1->fh_lda3_e1->Fill(lda3_e1);
    dir1->fh_lda3_e2->Fill(lda3_e2);
    dir1->fh_muon_status->Fill(sevt->muon[sevt->track[imu].iMuon].status);
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
    dir1->ComputeThreeTrack(k3pi_pion1,k3pi_pion2,k3pi_pion3);
    dir1->Fill3pi(dir1->GetThreeTrackMomentum());
    dir1->ComputeThreeTrack(electron1,electron2,muon);
    dir1->FillHist(dir1->GetThreeTrackMomentum(),dir1->GetNuMomentum(), Kcharge);
    dir1->FillAngle(muon.Momentum,dir1->GetTwoTrackMomentum());
    dir1->fh_HoDTotTime->Fill(Event_HoDTime);
    dir1->fh_DCHTotTime->Fill(Event_DCHTime);


    //Producing cut variable in more readable way with the class
    //described in Cuts.h;
    Cuts cutting = make_cuts(dir1,muon,electron1,electron2,Vertex_mu_e1,Vertex_mu_e2,Vertex_e1_e2,"munuee");

    //-- CUT1 DCH Geometry and Time Cut --
    if(cutting.DCH_Radius_mu < 14  ||
       cutting.DCH_Radius_mu > 110 ||
       cutting.DCH_Radius_el1 < 14 ||
       cutting.DCH_Radius_el1 > 110||
       cutting.DCH_Radius_el2 < 14 ||
       cutting.DCH_Radius_el2 > 110
       ){return 0;}

    if(cutting.Lkr_cut_el1  < 15.     ||
       cutting.Lkr_cut_el2  < 15.
       ){return 0;}
    if(   fabs(cutting.Lkr_x_el1) > 113 ||
          fabs(cutting.Lkr_x_el2) > 113 ||
          fabs(cutting.Lkr_y_el1) > 113 ||
          fabs(cutting.Lkr_y_el2) > 113
          ) {return 0;}
    if( (fabs(cutting.Lkr_x_el1) + fabs(cutting.Lkr_y_el1) ) > 159.8 ||
        (fabs(cutting.Lkr_x_el2) + fabs(cutting.Lkr_y_el2) ) > 159.8
        ){return 0;}

    // if(IS_DATA){
    if(cutting.Lkr_cut_mu  < 15.     ||
       fabs(cutting.Lkr_x_mu ) > 113 ||
       fabs(cutting.Lkr_y_mu ) > 113 ||
       (fabs(cutting.Lkr_x_mu) +  fabs(cutting.Lkr_y_mu ) ) > 159.8
       ){return 0;}
    //}

    if(fabs(cutting.MUV_y_mu) > 130. || fabs(cutting.MUV_x_mu) > 130.) {return 0;}
    if(fabs(cutting.MUV_x_mu) < 13  && fabs(cutting.MUV_y_mu) < 13  ) {return 0;}


    //-- CUT1 DCH Geometry and Time Cut --

    if(IS_MC){
        FillMC(dir3, True_Momentum[1], True_Momentum[2], True_Momentum[3], DKaon, Particle_production_zvtx, Particle_decay_zvtx);
        if(Npart >= 4){
            FillMC(dir3, True_Momentum[1], True_Momentum[2], True_Momentum[3], True_Momentum[4], DKaon, Particle_production_zvtx, Particle_decay_zvtx);

        }
    }


    dir3->fh_HoDTotTime->Fill(Event_HoDTime);
    dir3->fh_DCHTotTime->Fill(Event_DCHTime);

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
    dir3->ComputeThreeTrack(k3pi_pion1,k3pi_pion2,k3pi_pion3);
    dir3->Fill3pi(dir3->GetThreeTrackMomentum());
    dir3->ComputeThreeTrack(electron1,electron2,muon);
    dir3->FillHist(dir3->GetThreeTrackMomentum(),dir3->GetNuMomentum(), Kcharge);
    dir3->FillAngle(muon.Momentum,dir3->GetTwoTrackMomentum());
    dir3->fh_lda3_e1->Fill(lda3_e1);
    dir3->fh_lda3_e2->Fill(lda3_e2);
    dir3->fh_muon_status->Fill(sevt->muon[sevt->track[imu].iMuon].status);
    //-- CUT2 Momentum cut ---
    //if(cutting.Mu_P < 10. || cutting.Mu_P > 50.){return 0;}
    if(cutting.E1_P < 3.  || cutting.E1_P > 50.){return 0;}
    if(cutting.E2_P < 3.  || cutting.E2_P > 50.){return 0;}
    if(cutting.muee_P > 66){return 0;}
    //--ENDOF CUT2 Momentum cut ---

    if(IS_MC){
        FillMC(dir4, True_Momentum[1], True_Momentum[2], True_Momentum[3], DKaon, Particle_production_zvtx, Particle_decay_zvtx);
        if(Npart >= 4){
            FillMC(dir4, True_Momentum[1], True_Momentum[2], True_Momentum[3], True_Momentum[4], DKaon, Particle_production_zvtx, Particle_decay_zvtx);

        }
    }

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
    dir4->ComputeThreeTrack(k3pi_pion1,k3pi_pion2,k3pi_pion3);
    dir4->Fill3pi(dir4->GetThreeTrackMomentum());
    dir4->ComputeThreeTrack(electron1,electron2,muon);
    dir4->FillHist(dir4->GetThreeTrackMomentum(),dir4->GetNuMomentum(), Kcharge);
    dir4->FillAngle(muon.Momentum,dir4->GetTwoTrackMomentum());
    dir4->fh_lda3_e1->Fill(lda3_e1);
    dir4->fh_lda3_e2->Fill(lda3_e2);
    dir4->fh_muon_status->Fill(sevt->muon[sevt->track[imu].iMuon].status);


    //New 2D plots
    dir4->fh_mee_vs_MM2->Fill(dir4->GetTwoTrackMomentum().M(), dir4->GetNuMomentum().M2());
    dir4->fh_mee_vs_Pt->Fill(dir4->GetTwoTrackMomentum().M(), dir4->GetThreeTrackMomentum().Pt());
    dir4->fh_mee_vs_zvtx->Fill(dir4->GetTwoTrackMomentum().M(), COmPaCt_Z_Vertex);
    dir4->fh_mee_vs_Pel1->Fill(dir4->GetTwoTrackMomentum().M(),electron1.GetMomentum());
    dir4->fh_mee_vs_Pel2->Fill(dir4->GetTwoTrackMomentum().M(),electron2.GetMomentum());
    dir4->fh_mee_vs_Pmu->Fill(dir4->GetTwoTrackMomentum().M(),muon.GetMomentum());
    dir4->fh_P_vs_dtrkcl->Fill(electron1.GetMomentum(),electron1.GetDistanceTrackCluster() );
    dir4->fh_P_vs_dtrkcl->Fill(electron2.GetMomentum(),electron2.GetDistanceTrackCluster() );
    dir4->fh_P_vs_dtrkcl->Fill(muon.GetMomentum(),muon.GetDistanceTrackCluster() );
    dir4->fh_Ecl_vs_dtrkcl->Fill(electron1.GetEnergyLeftInEcal(),electron1.GetDistanceTrackCluster() );
    dir4->fh_Ecl_vs_dtrkcl->Fill(electron2.GetEnergyLeftInEcal(),electron2.GetDistanceTrackCluster() );
    dir4->fh_HoDTotTime->Fill(Event_HoDTime);
    dir4->fh_DCHTotTime->Fill(Event_DCHTime);



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

    //- CUT5 Invariant Mass Cut --
    if(cutting.mee < 0.140){return 0;}
    //--ENDOF CUT5 Invariant Mass Cut --

    if(IS_MC){
        FillMC(dir7, True_Momentum[1], True_Momentum[2], True_Momentum[3], DKaon, Particle_production_zvtx, Particle_decay_zvtx);
        if(Npart >= 4){
            FillMC(dir7, True_Momentum[1], True_Momentum[2], True_Momentum[3], True_Momentum[4], DKaon, Particle_production_zvtx, Particle_decay_zvtx);

        }
    }


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
    dir7->FillHist(dir7->GetThreeTrackMomentum(),dir7->GetNuMomentum(), Kcharge);
    dir7->FillAngle(muon.Momentum,dir7->GetTwoTrackMomentum());


    //-    //New 2D plots
    dir7->fh_mee_vs_MM2->Fill(dir7->GetTwoTrackMomentum().M(), dir7->GetNuMomentum().M2());
    dir7->fh_mee_vs_Pt->Fill(dir7->GetTwoTrackMomentum().M(), dir7->GetThreeTrackMomentum().Pt());
    dir7->fh_mee_vs_zvtx->Fill(dir7->GetTwoTrackMomentum().M(), COmPaCt_Z_Vertex);
    dir7->fh_mee_vs_Pel1->Fill(dir7->GetTwoTrackMomentum().M(),electron1.GetMomentum());
    dir7->fh_mee_vs_Pel2->Fill(dir7->GetTwoTrackMomentum().M(),electron2.GetMomentum());
    dir7->fh_mee_vs_Pmu->Fill(dir7->GetTwoTrackMomentum().M(),muon.GetMomentum());
    dir7->fh_P_vs_dtrkcl->Fill(electron1.GetMomentum(),electron1.GetDistanceTrackCluster() );
    dir7->fh_P_vs_dtrkcl->Fill(electron2.GetMomentum(),electron2.GetDistanceTrackCluster() );
    dir7->fh_P_vs_dtrkcl->Fill(muon.GetMomentum(),muon.GetDistanceTrackCluster() );
    dir7->fh_Ecl_vs_dtrkcl->Fill(electron1.GetEnergyLeftInEcal(),electron1.GetDistanceTrackCluster() );
    dir7->fh_Ecl_vs_dtrkcl->Fill(electron2.GetEnergyLeftInEcal(),electron2.GetDistanceTrackCluster() );
    dir7->fh_HoDTotTime->Fill(Event_HoDTime);
    dir7->fh_DCHTotTime->Fill(Event_DCHTime);



    //K3pi wrong sign selection
    //Charged_Particle k3pi_pion1(sevt,sbur,211,imu);
    //Charged_Particle k3pi_pion2(sevt,sbur,211,iel1);
    //Charged_Particle k3pi_pion3(sevt,sbur,211,iel2);

    dir7->ComputeThreeTrack(k3pi_pion1,k3pi_pion2,k3pi_pion3);
    dir7->fh_lda3_e1->Fill(lda3_e1);
    dir7->fh_lda3_e2->Fill(lda3_e2);
    dir7->Fill3pi(dir7->GetThreeTrackMomentum());
    dir7->fh_muon_status->Fill(sevt->muon[sevt->track[imu].iMuon].status);

    dir2->ComputeThreeTrack(electron1,electron2,muon);


    //- CUT5 Invariant Mass Cut --
    //if(cutting.mee < 0.100){return 0;}
    //--ENDOF CUT5 Invariant Mass Cut --

    if(electron1.GetCharge()==electron2.GetCharge() &&
       electron1.GetCharge()!=muon.GetCharge()      &&
       lda3_e1 > 0.8                                &&
       lda3_e2 > 0.8                                &&
       dir2->GetNuMomentum().M2() > -0.015          &&
       dir2->GetNuMomentum().M2() < 0.015           &&
       sevt->muon[sevt->track[imu].iMuon].status < 3 &&
       fabs(Kcharge) < 2){




        if(IS_MC){
            FillMC(dir2, True_Momentum[1], True_Momentum[2], True_Momentum[3], DKaon, Particle_production_zvtx, Particle_decay_zvtx);
            if(Npart >= 4){
                FillMC(dir2, True_Momentum[1], True_Momentum[2], True_Momentum[3], True_Momentum[4], DKaon, Particle_production_zvtx, Particle_decay_zvtx);

            }
        }



        TLorentzVector WS_e1nu = electron1.Momentum + dir2->GetNuMomentum();
        TLorentzVector WS_e2nu = electron2.Momentum + dir2->GetNuMomentum();
        TLorentzVector WS_munu = muon.Momentum + dir2->GetNuMomentum();
        TLorentzVector WS_three_pi = k3pi_pion1.Momentum + k3pi_pion2.Momentum + k3pi_pion3.Momentum;

        if(fabs(Event_HoDTime) > 3 &&
           fabs(Event_HoDTime) < 23 ){

            dir19->fh_missing_mass->Fill(dir2->GetNuMomentum().M2(),w_tot );
            dir19->fh_muee_Pt->Fill(dir2->GetThreeTrackMomentum().Pt(),w_tot);
            dir19->fh_COmPaCt_Z_Vertex->Fill(COmPaCt_Z_Vertex,w_tot);
            dir19->fh_Mu_momentum->Fill(muon.GetMomentum(),w_tot);
            dir19->fh_mee->Fill(dir2->GetTwoTrackMomentum().M(),w_tot);
            dir19->fh_mee_z_variable->Fill(dir2->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
            dir19->fh_HoDTotTime->Fill(Event_HoDTime,w_tot);

            dir19->fh_muee_M->Fill(dir2->GetThreeTrackMomentum().M(),w_tot);
            dir19->fh_Pmu_vs_z->Fill(muon.GetMomentum(), dir2->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
            dir19->fh_Kaon_Charge->Fill(Kcharge );
            dir19->fh_eop->Fill(muon.GetEnergyLeftInEcal()/muon.GetMomentum(),w_tot);
            dir19->fh_eop->Fill(electron1.GetEnergyLeftInEcal()/electron1.GetMomentum(),w_tot);
            dir19->fh_eop->Fill(electron2.GetEnergyLeftInEcal()/electron2.GetMomentum(),w_tot);
            dir19->fh_Muee_3pi->Fill( WS_three_pi.M(),w_tot);
            dir19->fh_M3pi_vs_EoP->Fill( WS_three_pi.M(), electron2.GetEnergyLeftInEcal()/electron2.GetMomentum(),w_tot);
            dir19->fh_M3pi_vs_EoP->Fill( WS_three_pi.M(), electron1.GetEnergyLeftInEcal()/electron1.GetMomentum(),w_tot);
            dir19->fh_munu_M->Fill(WS_munu.M(),w_tot);
            dir19->fh_e1nu_M->Fill(WS_e1nu.M(),w_tot);
            dir19->fh_e2nu_M->Fill(WS_e2nu.M(),w_tot);

            if(Kcharge>0){
                dir19->fh_mee_z_Kplus->Fill(dir2->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
                dir19->fh_MM2_plus->Fill(dir2->GetNuMomentum().M2(),w_tot );

            }

            if(Kcharge<0){
                dir19->fh_MM2_minus->Fill(dir2->GetNuMomentum().M2(),w_tot );
                dir19->fh_mee_z_Kminus->Fill(dir2->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
            }

            if(magnetCurrent < 0){
                dir19->fh_mee_z_magnet_minus->Fill(dir2->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
            }
            if(magnetCurrent > 0){
                dir19->fh_mee_z_magnet_plus->Fill(dir2->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
            }

        }


        if(fabs(Event_HoDTime) < 3){

            dir2->fh_missing_mass->Fill(dir2->GetNuMomentum().M2(),w_tot );
            dir2->fh_muee_Pt->Fill(dir2->GetThreeTrackMomentum().Pt(),w_tot);
            dir2->fh_muee_M->Fill(dir2->GetThreeTrackMomentum().M(),w_tot);
            dir2->fh_COmPaCt_Z_Vertex->Fill(COmPaCt_Z_Vertex,w_tot);
            dir2->fh_Mu_momentum->Fill(muon.GetMomentum(),w_tot);
            dir2->fh_mee->Fill(dir2->GetTwoTrackMomentum().M(),w_tot);
            dir2->fh_mee_z_variable->Fill(dir2->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
            dir2->fh_HoDTotTime->Fill(Event_HoDTime,w_tot);

            dir2->fh_Pmu_vs_z->Fill(muon.GetMomentum(), dir2->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
            dir2->fh_Kaon_Charge->Fill(Kcharge );
            dir2->fh_eop->Fill(muon.GetEnergyLeftInEcal()/muon.GetMomentum(),w_tot);
            dir2->fh_eop->Fill(electron1.GetEnergyLeftInEcal()/electron1.GetMomentum(),w_tot);
            dir2->fh_eop->Fill(electron2.GetEnergyLeftInEcal()/electron2.GetMomentum(),w_tot);
            dir2->fh_Muee_3pi->Fill( WS_three_pi.M(),w_tot);
            dir2->fh_M3pi_vs_EoP->Fill( WS_three_pi.M(), electron2.GetEnergyLeftInEcal()/electron2.GetMomentum(),w_tot);
            dir2->fh_M3pi_vs_EoP->Fill( WS_three_pi.M(), electron1.GetEnergyLeftInEcal()/electron1.GetMomentum(),w_tot);
            dir2->fh_munu_M->Fill(WS_munu.M(),w_tot);
            dir2->fh_e1nu_M->Fill(WS_e1nu.M(),w_tot);
            dir2->fh_e2nu_M->Fill(WS_e2nu.M(),w_tot);

            if(Kcharge>0){
                dir2->fh_mee_z_Kplus->Fill(dir2->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
                dir2->fh_MM2_plus->Fill(dir2->GetNuMomentum().M2(),w_tot );

            }

            if(Kcharge<0){
                dir2->fh_MM2_minus->Fill(dir2->GetNuMomentum().M2(),w_tot );
                dir2->fh_mee_z_Kminus->Fill(dir2->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
            }

            if(magnetCurrent < 0){
                dir2->fh_mee_z_magnet_minus->Fill(dir2->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
            }
            if(magnetCurrent > 0){
                dir2->fh_mee_z_magnet_plus->Fill(dir2->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
            }

            //Pt > 10 MeV/c
            if(dir2->GetThreeTrackMomentum().Pt() > 0.010){

                dir5->ComputeThreeTrack(electron1,electron2,muon);

                dir5->fh_missing_mass->Fill(dir2->GetNuMomentum().M2(),w_tot );
                dir5->fh_muee_Pt->Fill(dir2->GetThreeTrackMomentum().Pt(),w_tot);
                dir5->fh_COmPaCt_Z_Vertex->Fill(COmPaCt_Z_Vertex,w_tot);
                dir5->fh_Mu_momentum->Fill(muon.GetMomentum(),w_tot);
                dir5->fh_mee->Fill(dir2->GetTwoTrackMomentum().M(),w_tot);
                dir5->fh_mee_z_variable->Fill(dir2->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
                dir5->fh_HoDTotTime->Fill(Event_HoDTime,w_tot);

                dir5->fh_muee_M->Fill(dir2->GetThreeTrackMomentum().M(),w_tot);
                dir5->fh_Pmu_vs_z->Fill(muon.GetMomentum(), dir2->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
                dir5->fh_Kaon_Charge->Fill(Kcharge );
                dir5->fh_eop->Fill(muon.GetEnergyLeftInEcal()/muon.GetMomentum(),w_tot);
                dir5->fh_eop->Fill(electron1.GetEnergyLeftInEcal()/electron1.GetMomentum(),w_tot);
                dir5->fh_eop->Fill(electron2.GetEnergyLeftInEcal()/electron2.GetMomentum(),w_tot);
                dir5->fh_Muee_3pi->Fill( WS_three_pi.M(),w_tot);
                dir5->fh_M3pi_vs_EoP->Fill( WS_three_pi.M(), electron2.GetEnergyLeftInEcal()/electron2.GetMomentum(),w_tot);
                dir5->fh_M3pi_vs_EoP->Fill( WS_three_pi.M(), electron1.GetEnergyLeftInEcal()/electron1.GetMomentum(),w_tot);
                dir5->fh_munu_M->Fill(WS_munu.M(),w_tot);
                dir5->fh_e1nu_M->Fill(WS_e1nu.M(),w_tot);
                dir5->fh_e2nu_M->Fill(WS_e2nu.M(),w_tot);

                if(Kcharge>0){
                    dir5->fh_mee_z_Kplus->Fill(dir2->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
                    dir5->fh_MM2_plus->Fill(dir2->GetNuMomentum().M2(),w_tot );

                }

                if(Kcharge<0){
                    dir5->fh_MM2_minus->Fill(dir2->GetNuMomentum().M2(),w_tot );
                    dir5->fh_mee_z_Kminus->Fill(dir2->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
                }

                if(magnetCurrent < 0){
                    dir5->fh_mee_z_magnet_minus->Fill(dir2->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
                }
                if(magnetCurrent > 0){
                    dir5->fh_mee_z_magnet_plus->Fill(dir2->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
                }
            }

            dir6->ComputeThreeTrack(electron1,electron2,muon);
            //Pt > 20 MeV/c
            if(dir2->GetThreeTrackMomentum().Pt() > 0.020){


                dir6->fh_missing_mass->Fill(dir2->GetNuMomentum().M2(),w_tot );
                dir6->fh_muee_Pt->Fill(dir2->GetThreeTrackMomentum().Pt(),w_tot);
                dir6->fh_COmPaCt_Z_Vertex->Fill(COmPaCt_Z_Vertex,w_tot);
                dir6->fh_Mu_momentum->Fill(muon.GetMomentum(),w_tot);
                dir6->fh_mee->Fill(dir2->GetTwoTrackMomentum().M(),w_tot);
                dir6->fh_mee_z_variable->Fill(dir2->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
                dir6->fh_HoDTotTime->Fill(Event_HoDTime,w_tot);

                dir6->fh_muee_M->Fill(dir2->GetThreeTrackMomentum().M(),w_tot);
                dir6->fh_Pmu_vs_z->Fill(muon.GetMomentum(), dir2->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
                dir6->fh_Kaon_Charge->Fill(Kcharge );
                dir6->fh_eop->Fill(muon.GetEnergyLeftInEcal()/muon.GetMomentum(),w_tot);
                dir6->fh_eop->Fill(electron1.GetEnergyLeftInEcal()/electron1.GetMomentum(),w_tot);
                dir6->fh_eop->Fill(electron2.GetEnergyLeftInEcal()/electron2.GetMomentum(),w_tot);
                dir6->fh_Muee_3pi->Fill( WS_three_pi.M(),w_tot);
                dir6->fh_M3pi_vs_EoP->Fill( WS_three_pi.M(), electron2.GetEnergyLeftInEcal()/electron2.GetMomentum(),w_tot);
                dir6->fh_M3pi_vs_EoP->Fill( WS_three_pi.M(), electron1.GetEnergyLeftInEcal()/electron1.GetMomentum(),w_tot);
                dir6->fh_munu_M->Fill(WS_munu.M(),w_tot);
                dir6->fh_e1nu_M->Fill(WS_e1nu.M(),w_tot);
                dir6->fh_e2nu_M->Fill(WS_e2nu.M(),w_tot);

                if(Kcharge>0){
                    dir6->fh_mee_z_Kplus->Fill(dir2->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
                    dir6->fh_MM2_plus->Fill(dir2->GetNuMomentum().M2(),w_tot );

                }

                if(Kcharge<0){
                    dir6->fh_MM2_minus->Fill(dir2->GetNuMomentum().M2(),w_tot );
                    dir6->fh_mee_z_Kminus->Fill(dir2->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
                }

                if(magnetCurrent < 0){
                    dir6->fh_mee_z_magnet_minus->Fill(dir2->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
                }
                if(magnetCurrent > 0){
                    dir6->fh_mee_z_magnet_plus->Fill(dir2->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
                }
            }
            //Pt > 30 MeV/c
            if(dir2->GetThreeTrackMomentum().Pt() > 0.030){
                dir9->ComputeThreeTrack(electron1,electron2,muon);

                dir9->fh_missing_mass->Fill(dir2->GetNuMomentum().M2(),w_tot );
                dir9->fh_muee_Pt->Fill(dir2->GetThreeTrackMomentum().Pt(),w_tot);
                dir9->fh_COmPaCt_Z_Vertex->Fill(COmPaCt_Z_Vertex,w_tot);
                dir9->fh_Mu_momentum->Fill(muon.GetMomentum(),w_tot);
                dir9->fh_mee->Fill(dir2->GetTwoTrackMomentum().M(),w_tot);
                dir9->fh_mee_z_variable->Fill(dir2->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
                dir9->fh_HoDTotTime->Fill(Event_HoDTime,w_tot);

                dir9->fh_muee_M->Fill(dir2->GetThreeTrackMomentum().M(),w_tot);
                dir9->fh_Pmu_vs_z->Fill(muon.GetMomentum(), dir2->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
                dir9->fh_Kaon_Charge->Fill(Kcharge );
                dir9->fh_eop->Fill(muon.GetEnergyLeftInEcal()/muon.GetMomentum(),w_tot);
                dir9->fh_eop->Fill(electron1.GetEnergyLeftInEcal()/electron1.GetMomentum(),w_tot);
                dir9->fh_eop->Fill(electron2.GetEnergyLeftInEcal()/electron2.GetMomentum(),w_tot);
                dir9->fh_Muee_3pi->Fill( WS_three_pi.M(),w_tot);
                dir9->fh_M3pi_vs_EoP->Fill( WS_three_pi.M(), electron2.GetEnergyLeftInEcal()/electron2.GetMomentum(),w_tot);
                dir9->fh_M3pi_vs_EoP->Fill( WS_three_pi.M(), electron1.GetEnergyLeftInEcal()/electron1.GetMomentum(),w_tot);
                dir9->fh_munu_M->Fill(WS_munu.M(),w_tot);
                dir9->fh_e1nu_M->Fill(WS_e1nu.M(),w_tot);
                dir9->fh_e2nu_M->Fill(WS_e2nu.M(),w_tot);

                if(Kcharge>0){
                    dir9->fh_mee_z_Kplus->Fill(dir2->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
                    dir9->fh_MM2_plus->Fill(dir2->GetNuMomentum().M2(),w_tot );

                }

                if(Kcharge<0){
                    dir9->fh_MM2_minus->Fill(dir2->GetNuMomentum().M2(),w_tot );
                    dir9->fh_mee_z_Kminus->Fill(dir2->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
                }

                if(magnetCurrent < 0){
                    dir9->fh_mee_z_magnet_minus->Fill(dir2->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
                }
                if(magnetCurrent > 0){
                    dir9->fh_mee_z_magnet_plus->Fill(dir2->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
                }
            }


        }
    }



    //-- CUT7 K3pi invariant mass cut ---
    if(IS_DATA)
        if(lda3_e1 < 0.8 || lda3_e2 < 0.8){return 0;}
    if(COmPaCt_Z_Vertex < -1800. || COmPaCt_Z_Vertex > 8000.){return 0;}
    if(sevt->muon[sevt->track[imu].iMuon].status > 2 ) {return 0;}
    //if(dir10->GetThreeTrackMomentum().M() <= 0.51 ){return 0;}
    dir10->ComputeThreeTrack(k3pi_pion1,k3pi_pion2,k3pi_pion3);
    dir10->Fill3pi(dir10->GetThreeTrackMomentum());
    //--END OF CUT7 K3pi invariant mass cut ---


    dir10->fh_Event_Type->Fill(Event_Type);
    dir10->fh_Kaon_Charge->Fill(Kcharge);
    dir10->FillCommonHist(sevt);

    dir10->FillVertexHist(Vertex_mu_e1, cda_mu_e1 , Vertex_mu_e2, cda_mu_e2, Vertex_e1_e2, cda_e1_e2,"munuee");
    dir10->FillHist(muon,"muon");
    dir10->FillHist(electron1,"electron1");
    dir10->FillHist(electron2,"electron2");
    dir10->FillHist(muon,electron1,"mue1");
    dir10->FillHist(muon,electron2,"mue2");
    dir10->FillHist(electron1,electron2,"e1e2");
    dir10->ComputeThreeTrack(electron1,electron2,muon);
    dir10->FillHist(dir10->GetThreeTrackMomentum(),dir10->GetNuMomentum(), Kcharge);
    dir10->FillAngle(muon.Momentum,dir10->GetTwoTrackMomentum());
    dir10->fh_lda3_e1->Fill(lda3_e1);
    dir10->fh_lda3_e2->Fill(lda3_e2);
    dir10->fh_muon_status->Fill(sevt->muon[sevt->track[imu].iMuon].status);


    //-    //New 2D plots
    dir10->fh_mee_vs_MM2->Fill(dir10->GetTwoTrackMomentum().M(), dir10->GetNuMomentum().M2());
    dir10->fh_mee_vs_Pt->Fill(dir10->GetTwoTrackMomentum().M(), dir10->GetThreeTrackMomentum().Pt());
    dir10->fh_mee_vs_zvtx->Fill(dir10->GetTwoTrackMomentum().M(), COmPaCt_Z_Vertex);
    dir10->fh_mee_vs_Pel1->Fill(dir10->GetTwoTrackMomentum().M(),electron1.GetMomentum());
    dir10->fh_mee_vs_Pel2->Fill(dir10->GetTwoTrackMomentum().M(),electron2.GetMomentum());
    dir10->fh_mee_vs_Pmu->Fill(dir10->GetTwoTrackMomentum().M(),muon.GetMomentum());
    dir10->fh_P_vs_dtrkcl->Fill(electron1.GetMomentum(),electron1.GetDistanceTrackCluster() );
    dir10->fh_P_vs_dtrkcl->Fill(electron2.GetMomentum(),electron2.GetDistanceTrackCluster() );
    dir10->fh_P_vs_dtrkcl->Fill(muon.GetMomentum(),muon.GetDistanceTrackCluster() );
    dir10->fh_Ecl_vs_dtrkcl->Fill(electron1.GetEnergyLeftInEcal(),electron1.GetDistanceTrackCluster() );
    dir10->fh_Ecl_vs_dtrkcl->Fill(electron2.GetEnergyLeftInEcal(),electron2.GetDistanceTrackCluster() );

    dir10->fh_HoDTotTime->Fill(Event_HoDTime);
    dir10->fh_DCHTotTime->Fill(Event_DCHTime);


    //LAST CUT z----- vtx -------------
    if(dir10->GetNuMomentum().M2() < -0.015 || dir10->GetNuMomentum().M2() > 0.015){return 0;}

    dir11->ComputeThreeTrack(k3pi_pion1,k3pi_pion2,k3pi_pion3);
    dir11->Fill3pi(dir11->GetThreeTrackMomentum());
    //--END OF CUT z --------- vtx ---
    double trigger = 1;
    // ------------------------------------------ SS0 --------------------------------------//
    /* if( 15304 < sbur->nrun < 15582){ */
    if( 15304 <= sbur->nrun && sbur->nrun <= 15582){
        if(((sevt->trigWord)>>3)&0x1) {
            // ref. trigg. ok:
            dir11->fh_SS0_CPRE->Fill(trigger);
            if( (sevt->pu[3].chan[5]>>4)&0x1 ||
                (sevt->pu[4].chan[5]>>4)&0x1 ||
                (sevt->pu[5].chan[5]>>4)&0x1 ||  //5 is the central time slice
                (sevt->pu[6].chan[5]>>4)&0x1 ||
                (sevt->pu[3].chan[5]>>8)&0x1 ||
                (sevt->pu[4].chan[5]>>8)&0x1 ||
                (sevt->pu[5].chan[5]>>8)&0x1 ||  //5 is the central time slice
                (sevt->pu[6].chan[5]>>8)&0x1 ||
                (sevt->pu[3].chan[5]>>13)&0x1 ||
                (sevt->pu[4].chan[5]>>13)&0x1 ||
                (sevt->pu[5].chan[5]>>13)&0x1 ||  //5 is the central time slice
                (sevt->pu[6].chan[5]>>13)&0x1 ) {
                if( ((sevt->trigWord)>>4)&0x1 ) //
                    dir11->fh_SS0_MB_1TRK_P->Fill(trigger);
                if( ((sevt->trigWord)>>1)&0x1  )
                    dir11->fh_SS0_MB_1VTX->Fill(trigger);

                if( ((sevt->trigWord))&0x1    || ((sevt->trigWord)>>1)&0x1 || ((sevt->trigWord)>>4)&0x1 )
                    dir11->fh_SS0_full_trig->Fill(trigger);

            }
        }
    }
    //--------------------------------------- SS1 ---------------------------------------//
    if( 15633 <= sbur->nrun && sbur->nrun <= 15703){
        if(((sevt->trigWord)>>3)&0x1) {
            // ref. trigg. ok:
            dir11->fh_SS1_CPRE->Fill(trigger);
            if( (sevt->pu[3].chan[5]>>4)&0x1 ||
                (sevt->pu[4].chan[5]>>4)&0x1 ||
                (sevt->pu[5].chan[5]>>4)&0x1 ||  //5 is the central time slice
                (sevt->pu[6].chan[5]>>4)&0x1 ||
                (sevt->pu[3].chan[5]>>8)&0x1 ||
                (sevt->pu[4].chan[5]>>8)&0x1 ||
                (sevt->pu[5].chan[5]>>8)&0x1 ||  //5 is the central time slice
                (sevt->pu[6].chan[5]>>8)&0x1 ||
                (sevt->pu[3].chan[5]>>13)&0x1 ||
                (sevt->pu[4].chan[5]>>13)&0x1 ||
                (sevt->pu[5].chan[5]>>13)&0x1 ||  //5 is the central time slice
                (sevt->pu[6].chan[5]>>13)&0x1 ) {
                if( ((sevt->trigWord)>>4)&0x1 ) //
                    dir11->fh_SS1_MB_1TRK_P->Fill(trigger);
                if( ((sevt->trigWord)>>1)&0x1  )
                    dir11->fh_SS1_MB_1VTX->Fill(trigger);
                if( ((sevt->trigWord))&0x1    || ((sevt->trigWord)>>1)&0x1 || ((sevt->trigWord)>>4)&0x1 )
                    dir11->fh_SS1_full_trig->Fill(trigger);
            }
        }
    }
    //------------------------------------- SS2 -----------------------------------------------//
    if( 15717 <= sbur->nrun && sbur->nrun <= 15777){
        if(((sevt->trigWord)>>3)&0x1) {
            // ref. trigg. ok:
            dir11->fh_SS2_CPRE->Fill(trigger);
            if( (sevt->pu[3].chan[5]>>4)&0x1 ||
                (sevt->pu[4].chan[5]>>4)&0x1 ||
                (sevt->pu[5].chan[5]>>4)&0x1 ||  //5 is the central time slice
                (sevt->pu[6].chan[5]>>4)&0x1 ||
                (sevt->pu[3].chan[5]>>8)&0x1 ||
                (sevt->pu[4].chan[5]>>8)&0x1 ||
                (sevt->pu[5].chan[5]>>8)&0x1 ||  //5 is the central time slice
                (sevt->pu[6].chan[5]>>8)&0x1 ||
                (sevt->pu[3].chan[5]>>13)&0x1 ||
                (sevt->pu[4].chan[5]>>13)&0x1 ||
                (sevt->pu[5].chan[5]>>13)&0x1 ||  //5 is the central time slice
                (sevt->pu[6].chan[5]>>13)&0x1 ) {
                if( ((sevt->trigWord)>>4)&0x1 ) //
                    dir11->fh_SS2_MB_1TRK_P->Fill(trigger);
                if( ((sevt->trigWord)>>1)&0x1  )
                    dir11->fh_SS2_MB_1VTX->Fill(trigger);
                if( ((sevt->trigWord))&0x1    || ((sevt->trigWord)>>1)&0x1 || ((sevt->trigWord)>>4)&0x1 )
                    dir11->fh_SS2_full_trig->Fill(trigger);
            }
        }
    }
    //------------------------------------- SS3 -----------------------------------------------//
    if( 15778 <= sbur->nrun && sbur->nrun <= 15790){
        if(((sevt->trigWord)>>3)&0x1) {
            // ref. trigg. ok:
            dir11->fh_SS3_CPRE->Fill(trigger);
            if( (sevt->pu[3].chan[5]>>4)&0x1 ||
                (sevt->pu[4].chan[5]>>4)&0x1 ||
                (sevt->pu[5].chan[5]>>4)&0x1 ||  //5 is the central time slice
                (sevt->pu[6].chan[5]>>4)&0x1 ||
                (sevt->pu[3].chan[5]>>8)&0x1 ||
                (sevt->pu[4].chan[5]>>8)&0x1 ||
                (sevt->pu[5].chan[5]>>8)&0x1 ||  //5 is the central time slice
                (sevt->pu[6].chan[5]>>8)&0x1 ||
                (sevt->pu[3].chan[5]>>13)&0x1 ||
                (sevt->pu[4].chan[5]>>13)&0x1 ||
                (sevt->pu[5].chan[5]>>13)&0x1 ||  //5 is the central time slice
                (sevt->pu[6].chan[5]>>13)&0x1 ) {
                if( ((sevt->trigWord)>>4)&0x1 ) //
                    dir11->fh_SS3_MB_1TRK_P->Fill(trigger);
                if( ((sevt->trigWord)>>1)&0x1  )
                    dir11->fh_SS3_MB_1VTX->Fill(trigger);
                if( ((sevt->trigWord))&0x1    || ((sevt->trigWord)>>1)&0x1 || ((sevt->trigWord)>>4)&0x1 )
                    dir11->fh_SS3_full_trig->Fill(trigger);
            }
        }
    }

    //------------------------------------- SS4 -----------------------------------------------//
    if( 16121 <= sbur->nrun && sbur->nrun <= 16383){
        if(((sevt->trigWord)>>3)&0x1) {
            // ref. trigg. ok:
            dir11->fh_SS4_CPRE->Fill(trigger);
            if( (sevt->pu[3].chan[5]>>4)&0x1 ||
                (sevt->pu[4].chan[5]>>4)&0x1 ||
                (sevt->pu[5].chan[5]>>4)&0x1 ||  //5 is the central time slice
                (sevt->pu[6].chan[5]>>4)&0x1 ||
                (sevt->pu[3].chan[5]>>8)&0x1 ||
                (sevt->pu[4].chan[5]>>8)&0x1 ||
                (sevt->pu[5].chan[5]>>8)&0x1 ||  //5 is the central time slice
                (sevt->pu[6].chan[5]>>8)&0x1 ||
                (sevt->pu[3].chan[5]>>13)&0x1 ||
                (sevt->pu[4].chan[5]>>13)&0x1 ||
                (sevt->pu[5].chan[5]>>13)&0x1 ||  //5 is the central time slice
                (sevt->pu[6].chan[5]>>13)&0x1 ) {
                if( ((sevt->trigWord)>>4)&0x1 ) //
                    dir11->fh_SS4_MB_1TRK_P->Fill(trigger);
                if( ((sevt->trigWord)>>1)&0x1  )
                    dir11->fh_SS4_MB_1VTX->Fill(trigger);
                if( ((sevt->trigWord))&0x1    || ((sevt->trigWord)>>1)&0x1 || ((sevt->trigWord)>>4)&0x1 )
                    dir11->fh_SS4_full_trig->Fill(trigger);
            }
        }
    }


    //------------------------------------- SS5 -----------------------------------------------//
    if( 16428 <= sbur->nrun && sbur->nrun <= 16585){
        if(((sevt->trigWord)>>3)&0x1) {
            // ref. trigg. ok:
            dir11->fh_SS5_CPRE->Fill(trigger);
            if( (sevt->pu[3].chan[5]>>4)&0x1 ||
                (sevt->pu[4].chan[5]>>4)&0x1 ||
                (sevt->pu[5].chan[5]>>4)&0x1 ||  //5 is the central time slice
                (sevt->pu[6].chan[5]>>4)&0x1 ||
                (sevt->pu[3].chan[5]>>8)&0x1 ||
                (sevt->pu[4].chan[5]>>8)&0x1 ||
                (sevt->pu[5].chan[5]>>8)&0x1 ||  //5 is the central time slice
                (sevt->pu[6].chan[5]>>8)&0x1 ||
                (sevt->pu[3].chan[5]>>13)&0x1 ||
                (sevt->pu[4].chan[5]>>13)&0x1 ||
                (sevt->pu[5].chan[5]>>13)&0x1 ||  //5 is the central time slice
                (sevt->pu[6].chan[5]>>13)&0x1 ) {
                if( ((sevt->trigWord)>>4)&0x1 ) //
                    dir11->fh_SS5_MB_1TRK_P->Fill(trigger);
                if( ((sevt->trigWord)>>1)&0x1  )
                    dir11->fh_SS5_MB_1VTX->Fill(trigger);
                if( ((sevt->trigWord))&0x1    || ((sevt->trigWord)>>1)&0x1 || ((sevt->trigWord)>>4)&0x1 )
                    dir11->fh_SS5_full_trig->Fill(trigger);
            }
        }
    }

    //------------------------------------- SS6 -----------------------------------------------//
    if( 16586 <= sbur->nrun && sbur->nrun <= 16709){
        if(((sevt->trigWord)>>3)&0x1) {
            // ref. trigg. ok:
            dir11->fh_SS6_CPRE->Fill(trigger);
            if( (sevt->pu[3].chan[5]>>4)&0x1 ||
                (sevt->pu[4].chan[5]>>4)&0x1 ||
                (sevt->pu[5].chan[5]>>4)&0x1 ||  //5 is the central time slice
                (sevt->pu[6].chan[5]>>4)&0x1 ||
                (sevt->pu[3].chan[5]>>8)&0x1 ||
                (sevt->pu[4].chan[5]>>8)&0x1 ||
                (sevt->pu[5].chan[5]>>8)&0x1 ||  //5 is the central time slice
                (sevt->pu[6].chan[5]>>8)&0x1 ||
                (sevt->pu[3].chan[5]>>13)&0x1 ||
                (sevt->pu[4].chan[5]>>13)&0x1 ||
                (sevt->pu[5].chan[5]>>13)&0x1 ||  //5 is the central time slice
                (sevt->pu[6].chan[5]>>13)&0x1 ) {
                if( ((sevt->trigWord)>>4)&0x1 ) //
                    dir11->fh_SS6_MB_1TRK_P->Fill(trigger);
                if( ((sevt->trigWord)>>1)&0x1  )
                    dir11->fh_SS6_MB_1VTX->Fill(trigger);
                if( ((sevt->trigWord))&0x1    || ((sevt->trigWord)>>1)&0x1 || ((sevt->trigWord)>>4)&0x1 )
                    dir11->fh_SS6_full_trig->Fill(trigger);
            }
        }
    }

    //------------------------------------- SS7 -----------------------------------------------//
    if( 16722 <= sbur->nrun && sbur->nrun <= 16801){
        if(((sevt->trigWord)>>3)&0x1) {
            // ref. trigg. ok:
            dir11->fh_SS7_CPRE->Fill(trigger);
            if( (sevt->pu[3].chan[5]>>4)&0x1 ||
                (sevt->pu[4].chan[5]>>4)&0x1 ||
                (sevt->pu[5].chan[5]>>4)&0x1 ||  //5 is the central time slice
                (sevt->pu[6].chan[5]>>4)&0x1 ||
                (sevt->pu[3].chan[5]>>8)&0x1 ||
                (sevt->pu[4].chan[5]>>8)&0x1 ||
                (sevt->pu[5].chan[5]>>8)&0x1 ||  //5 is the central time slice
                (sevt->pu[6].chan[5]>>8)&0x1 ||
                (sevt->pu[3].chan[5]>>13)&0x1 ||
                (sevt->pu[4].chan[5]>>13)&0x1 ||
                (sevt->pu[5].chan[5]>>13)&0x1 ||  //5 is the central time slice
                (sevt->pu[6].chan[5]>>13)&0x1 ) {
                if( ((sevt->trigWord)>>4)&0x1 ) //
                    dir11->fh_SS7_MB_1TRK_P->Fill(trigger);
                if( ((sevt->trigWord)>>1)&0x1  )
                    dir11->fh_SS7_MB_1VTX->Fill(trigger);
                if( ((sevt->trigWord))&0x1    || ((sevt->trigWord)>>1)&0x1 || ((sevt->trigWord)>>4)&0x1 )
                    dir11->fh_SS7_full_trig->Fill(trigger);
            }
        }
    }

    //------------------------------------- SS8 -----------------------------------------------//
    if( 16802 <= sbur->nrun && sbur->nrun <= 16905){
        if(((sevt->trigWord)>>3)&0x1) {
            // ref. trigg. ok:
            dir11->fh_SS8_CPRE->Fill(trigger);
            if( (sevt->pu[3].chan[5]>>4)&0x1 ||
                (sevt->pu[4].chan[5]>>4)&0x1 ||
                (sevt->pu[5].chan[5]>>4)&0x1 ||  //5 is the central time slice
                (sevt->pu[6].chan[5]>>4)&0x1 ||
                (sevt->pu[3].chan[5]>>8)&0x1 ||
                (sevt->pu[4].chan[5]>>8)&0x1 ||
                (sevt->pu[5].chan[5]>>8)&0x1 ||  //5 is the central time slice
                (sevt->pu[6].chan[5]>>8)&0x1 ||
                (sevt->pu[3].chan[5]>>13)&0x1 ||
                (sevt->pu[4].chan[5]>>13)&0x1 ||
                (sevt->pu[5].chan[5]>>13)&0x1 ||  //5 is the central time slice
                (sevt->pu[6].chan[5]>>13)&0x1 ) {
                if( ((sevt->trigWord)>>4)&0x1 ) //
                    dir11->fh_SS8_MB_1TRK_P->Fill(trigger);
                if( ((sevt->trigWord)>>1)&0x1  )
                    dir11->fh_SS8_MB_1VTX->Fill(trigger);
                if( ((sevt->trigWord))&0x1    || ((sevt->trigWord)>>1)&0x1 || ((sevt->trigWord)>>4)&0x1 )
                    dir11->fh_SS8_full_trig->Fill(trigger);
            }
        }
    }


    dir11->fh_Event_Type->Fill(Event_Type);
    dir11->fh_Kaon_Charge->Fill(Kcharge);
    dir11->FillCommonHist(sevt);
    dir11->FillVertexHist(Vertex_mu_e1, cda_mu_e1 , Vertex_mu_e2, cda_mu_e2, Vertex_e1_e2, cda_e1_e2,"munuee");
    dir11->FillHist(muon,"muon");
    dir11->FillHist(electron1,"electron1");
    dir11->FillHist(electron2,"electron2");
    dir11->FillHist(muon,electron1,"mue1");
    dir11->FillHist(muon,electron2,"mue2");
    dir11->FillHist(electron1,electron2,"e1e2");
    dir11->ComputeThreeTrack(electron1,electron2,muon);
    dir11->FillHist(dir11->GetThreeTrackMomentum(),dir11->GetNuMomentum(), Kcharge);
    dir11->FillAngle(muon.Momentum,dir11->GetTwoTrackMomentum());
    dir11->fh_lda3_e1->Fill(lda3_e1);
    dir11->fh_lda3_e2->Fill(lda3_e2);
    dir11->fh_muon_status->Fill(sevt->muon[sevt->track[imu].iMuon].status);


    //-    //New 2D plots
    dir11->fh_mee_vs_MM2->Fill(dir11->GetTwoTrackMomentum().M(), dir11->GetNuMomentum().M2());
    dir11->fh_mee_vs_Pt->Fill(dir11->GetTwoTrackMomentum().M(), dir11->GetThreeTrackMomentum().Pt());
    dir11->fh_mee_vs_zvtx->Fill(dir11->GetTwoTrackMomentum().M(), COmPaCt_Z_Vertex);
    dir11->fh_mee_vs_Pel1->Fill(dir11->GetTwoTrackMomentum().M(),electron1.GetMomentum());
    dir11->fh_mee_vs_Pel2->Fill(dir11->GetTwoTrackMomentum().M(),electron2.GetMomentum());
    dir11->fh_mee_vs_Pmu->Fill(dir11->GetTwoTrackMomentum().M(),muon.GetMomentum());
    dir11->fh_P_vs_dtrkcl->Fill(electron1.GetMomentum(),electron1.GetDistanceTrackCluster() );
    dir11->fh_P_vs_dtrkcl->Fill(electron2.GetMomentum(),electron2.GetDistanceTrackCluster() );
    dir11->fh_P_vs_dtrkcl->Fill(muon.GetMomentum(),muon.GetDistanceTrackCluster() );
    dir11->fh_Ecl_vs_dtrkcl->Fill(electron1.GetEnergyLeftInEcal(),electron1.GetDistanceTrackCluster() );
    dir11->fh_Ecl_vs_dtrkcl->Fill(electron2.GetEnergyLeftInEcal(),electron2.GetDistanceTrackCluster() );
    dir11->fh_HoDTotTime->Fill(Event_HoDTime);
    dir11->fh_DCHTotTime->Fill(Event_DCHTime);


    if(IS_MC){
        FillMC(dir10, True_Momentum[1], True_Momentum[2], True_Momentum[3], DKaon, Particle_production_zvtx, Particle_decay_zvtx);
        if(Npart >= 4){
            FillMC(dir10, True_Momentum[1], True_Momentum[2], True_Momentum[3], True_Momentum[4], DKaon, Particle_production_zvtx, Particle_decay_zvtx);

        }
    }


    //noSelectedEvents++;
    //TEST
    if(IS_MC){
        //electron1 reweighting
        if(electron1.Momentum.P() < 5.){
            //cout << "Pel1 - 0-5 GeV " << endl;
            w_el1 = w_0_5;
        }  else if(electron1.Momentum.P() >= 5. && electron1.Momentum.P() < 10.){
            //cout << "Pel1 - 5-10 GeV " << endl;
            w_el1 = w_5_10;
        } else if(electron1.Momentum.P() > 10. && electron1.Momentum.P() < 15.){
            //cout << "Pel1 - 10-15 GeV " << endl;
            w_el1 = w_10_15;
        } else if(electron1.Momentum.P() >= 15. && electron1.Momentum.P() < 20.){
            //cout << "Pel1 - 15-20 GeV " << endl;
            w_el1 = w_15_20;
        } else if(electron1.Momentum.P() >= 20. && electron1.Momentum.P() < 25.){
            //cout << "Pel1 - 20-25 GeV " << endl;
            w_el1 = w_20_25;
        } else if(electron1.Momentum.P() >= 25. && electron1.Momentum.P() < 30.){
            //cout << "Pel1 - 25-30 GeV " << endl;
            w_el1 = w_25_30;
        } else if(electron1.Momentum.P() >= 30. && electron1.Momentum.P() < 35.){
            //cout << "Pel1 - 30-35 GeV " << endl;
            w_el1 = w_30_35;
        } else if(electron1.Momentum.P() >= 35. && electron1.Momentum.P() < 40.){
            //cout << "Pel1 - 35-40 GeV " << endl;
            w_el1 = w_35_40;
        }

        //electron2 reweighting
        if(electron2.Momentum.P() < 5.){
            //cout << "Pel2 - 5-10 GeV " << endl;
            w_el2 = w_0_5;
        } else if(electron2.Momentum.P() >= 5. && electron2.Momentum.P() < 10.){
            //cout << "Pel2 - 5-10 GeV " << endl;
            w_el2 = w_5_10;
        } else if(electron2.Momentum.P() >= 10. && electron2.Momentum.P() < 15.){
            //cout << "Pel2 - 10-15 GeV " << endl;
            w_el2 = w_10_15;
        } else if(electron2.Momentum.P() >= 15. && electron2.Momentum.P() < 20.){
            //cout << "Pel2 - 15-20 GeV " << endl;
            w_el2 = w_15_20;
        } else if(electron2.Momentum.P() >= 20. && electron2.Momentum.P() < 25.){
            //cout << "Pel2 - 20-25 GeV " << endl;
            w_el2 = w_20_25;
        } else if(electron2.Momentum.P() >= 25. && electron2.Momentum.P() < 30.){
            //cout << "Pel2 - 25-30 GeV " << endl;
            w_el2 = w_25_30;
        } else if(electron2.Momentum.P() >= 30. && electron2.Momentum.P() < 35.){
            //cout << "Pel2 - 30-35 GeV " << endl;
            w_el2 = w_30_35;
        } else if(electron2.Momentum.P() >= 35. && electron2.Momentum.P() < 40.){
            //cout << "Pel2 - 30-35 GeV " << endl;
            w_el2 = w_35_40;
        }

        if(muon.Momentum.P() >= 10.){
            w_mu = 0.9957;
        }  else if(muon.Momentum.P() >= 9. && muon.Momentum.P() < 10.){
            w_mu = wmu_9_10;
        } else if(muon.Momentum.P() >= 8. && muon.Momentum.P() < 9.){
            w_mu = wmu_8_9;
        } else if(muon.Momentum.P() >= 7.5 && muon.Momentum.P() < 8.){
            w_mu = wmu_75_8;
        } else if(muon.Momentum.P() >= 7. && muon.Momentum.P() < 7.5){
            //cout << "Pmu - 30-35 GeV " << endl;
            w_mu = wmu_7_75;
        } else if(muon.Momentum.P() >= 6.5 && muon.Momentum.P() < 7.){
            w_mu = wmu_65_7;
        } else if(muon.Momentum.P() >= 6 && muon.Momentum.P() < 6.5){
            w_mu = wmu_6_65;
        } else if(muon.Momentum.P() >= 5.5 && muon.Momentum.P() < 6){
            w_mu = wmu_55_6;
        } else if(muon.Momentum.P() >= 5 && muon.Momentum.P() < 5.5){
            w_mu = wmu_5_55;
        }

        w_tot = w_el1*w_el2*w_mu;

        //if(Kcharge ==  -1) w_tot = w_tot*0.556;

        MC_reweight->fh_mee_z_variable->Fill(dir11->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
        MC_reweight->fh_mee->Fill(dir11->GetTwoTrackMomentum().M(),w_tot);
        MC_reweight->fh_missing_mass->Fill(dir11->GetNuMomentum().M2(),w_tot );
        MC_reweight->fh_muee_Pt->Fill(dir11->GetThreeTrackMomentum().Pt(),w_tot);
        MC_reweight->fh_COmPaCt_Z_Vertex->Fill(COmPaCt_Z_Vertex,w_tot);
        MC_reweight->fh_Mu_momentum->Fill(muon.GetMomentum(),w_tot);
    }
    //END OF TEST

    if(IS_MC){
        FillMC(dir11, True_Momentum[1], True_Momentum[2], True_Momentum[3], DKaon, Particle_production_zvtx, Particle_decay_zvtx);
        if(Npart >= 4){
            FillMC(dir11, True_Momentum[1], True_Momentum[2], True_Momentum[3], True_Momentum[4], DKaon, Particle_production_zvtx, Particle_decay_zvtx);

        }
    }

    if(fabs(Event_HoDTime) > 3 && fabs(Event_HoDTime) < 23 ){
        dir12->ComputeThreeTrack(k3pi_pion1,k3pi_pion2,k3pi_pion3);
        dir12->Fill3pi(dir12->GetThreeTrackMomentum());

        dir12->fh_Event_Type->Fill(Event_Type);
        dir12->fh_Kaon_Charge->Fill(Kcharge);
        dir12->FillCommonHist(sevt);

        dir12->FillVertexHist(Vertex_mu_e1, cda_mu_e1 , Vertex_mu_e2, cda_mu_e2, Vertex_e1_e2, cda_e1_e2,"munuee");
        dir12->FillHist(muon,"muon");
        dir12->FillHist(electron1,"electron1");
        dir12->FillHist(electron2,"electron2");
        dir12->FillHist(muon,electron1,"mue1");
        dir12->FillHist(muon,electron2,"mue2");
        dir12->FillHist(electron1,electron2,"e1e2");
        dir12->ComputeThreeTrack(electron1,electron2,muon);
        dir12->FillHist(dir12->GetThreeTrackMomentum(),dir12->GetNuMomentum(), Kcharge);
        dir12->FillAngle(muon.Momentum,dir12->GetTwoTrackMomentum());
        dir12->fh_lda3_e1->Fill(lda3_e1);
        dir12->fh_lda3_e2->Fill(lda3_e2);
        dir12->fh_muon_status->Fill(sevt->muon[sevt->track[imu].iMuon].status);


        //-    //New 2D plots
        dir12->fh_mee_vs_MM2->Fill(dir12->GetTwoTrackMomentum().M(), dir12->GetNuMomentum().M2());
        dir12->fh_mee_vs_Pt->Fill(dir12->GetTwoTrackMomentum().M(), dir12->GetThreeTrackMomentum().Pt());
        dir12->fh_mee_vs_zvtx->Fill(dir12->GetTwoTrackMomentum().M(), COmPaCt_Z_Vertex);
        dir12->fh_mee_vs_Pel1->Fill(dir12->GetTwoTrackMomentum().M(),electron1.GetMomentum());
        dir12->fh_mee_vs_Pel2->Fill(dir12->GetTwoTrackMomentum().M(),electron2.GetMomentum());
        dir12->fh_mee_vs_Pmu->Fill(dir12->GetTwoTrackMomentum().M(),muon.GetMomentum());
        dir12->fh_P_vs_dtrkcl->Fill(electron1.GetMomentum(),electron1.GetDistanceTrackCluster() );
        dir12->fh_P_vs_dtrkcl->Fill(electron2.GetMomentum(),electron2.GetDistanceTrackCluster() );
        dir12->fh_P_vs_dtrkcl->Fill(muon.GetMomentum(),muon.GetDistanceTrackCluster() );
        dir12->fh_Ecl_vs_dtrkcl->Fill(electron1.GetEnergyLeftInEcal(),electron1.GetDistanceTrackCluster() );
        dir12->fh_Ecl_vs_dtrkcl->Fill(electron2.GetEnergyLeftInEcal(),electron2.GetDistanceTrackCluster() );


        dir12->fh_HoDTotTime->Fill(Event_HoDTime);
        dir12->fh_DCHTotTime->Fill(Event_DCHTime);

    }


    if(IS_DATA)
        if(fabs(Event_HoDTime) > 3 ){return 0;}

    TLorentzVector mc_ee =  True_Momentum[2] + True_Momentum[3];
    TLorentzVector e1nu = electron1.Momentum + dir10->GetNuMomentum();
    TLorentzVector e2nu = electron2.Momentum + dir10->GetNuMomentum();
    TLorentzVector munu = muon.Momentum + dir10->GetNuMomentum();
    TLorentzVector three_pi = k3pi_pion1.Momentum + k3pi_pion2.Momentum + k3pi_pion3.Momentum;

    //Q vertex = +-3
    if(fabs(Kcharge) > 1 ){



        dir18->fh_missing_mass->Fill(dir10->GetNuMomentum().M2(),w_tot );
        dir18->fh_muee_Pt->Fill(dir10->GetThreeTrackMomentum().Pt(),w_tot);
        dir18->fh_COmPaCt_Z_Vertex->Fill(COmPaCt_Z_Vertex,w_tot);
        dir18->fh_Mu_momentum->Fill(muon.GetMomentum(),w_tot);
        dir18->fh_mee->Fill(dir10->GetTwoTrackMomentum().M(),w_tot);
        dir18->fh_mee_z_variable->Fill(dir10->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);

        dir18->fh_muee_M->Fill(dir10->GetThreeTrackMomentum().M(),w_tot);
        dir18->fh_Pmu_vs_z->Fill(muon.GetMomentum(), dir10->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
        dir18->fh_Kaon_Charge->Fill(Kcharge );
        dir18->fh_eop->Fill(muon.GetEnergyLeftInEcal()/muon.GetMomentum(),w_tot);
        dir18->fh_eop->Fill(electron1.GetEnergyLeftInEcal()/electron1.GetMomentum(),w_tot);
        dir18->fh_eop->Fill(electron2.GetEnergyLeftInEcal()/electron2.GetMomentum(),w_tot);
        dir18->fh_Muee_3pi->Fill( three_pi.M(),w_tot);
        dir18->fh_M3pi_vs_EoP->Fill( three_pi.M(), electron2.GetEnergyLeftInEcal()/electron2.GetMomentum(),w_tot);
        dir18->fh_M3pi_vs_EoP->Fill( three_pi.M(), electron1.GetEnergyLeftInEcal()/electron1.GetMomentum(),w_tot);
        dir18->fh_munu_M->Fill(munu.M(),w_tot);
        dir18->fh_e1nu_M->Fill(e1nu.M(),w_tot);
        dir18->fh_e2nu_M->Fill(e2nu.M(),w_tot);
        dir18->fh_mc_two_track_23_mass_z_variable->Fill(mc_ee.M()*mc_ee.M()/(0.493677*0.493677));

        dir18->fh_HoDTotTime->Fill(Event_HoDTime,w_tot);

        if(Kcharge>0){
            dir18->fh_mee_z_Kplus->Fill(dir10->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
            dir18->fh_MM2_plus->Fill(dir10->GetNuMomentum().M2(),w_tot );

        }

        if(Kcharge<0){
            dir18->fh_MM2_minus->Fill(dir10->GetNuMomentum().M2(),w_tot );
            dir18->fh_mee_z_Kminus->Fill(dir10->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
        }

        if(magnetCurrent < 0){
            dir18->fh_mee_z_magnet_minus->Fill(dir10->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
        }
        if(magnetCurrent > 0){
            dir18->fh_mee_z_magnet_plus->Fill(dir10->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
        }

    }


    if(fabs(Kcharge) > 1 ) return 0;

    dir13->ComputeThreeTrack(electron1,electron2,muon);

    dir13->fh_missing_mass->Fill(dir10->GetNuMomentum().M2(),w_tot );
    dir13->fh_muee_Pt->Fill(dir10->GetThreeTrackMomentum().Pt(),w_tot);
    dir13->fh_COmPaCt_Z_Vertex->Fill(COmPaCt_Z_Vertex,w_tot);
    dir13->fh_Mu_momentum->Fill(muon.GetMomentum(),w_tot);
    dir13->fh_mee->Fill(dir10->GetTwoTrackMomentum().M(),w_tot);
    dir13->fh_mee_z_variable->Fill(dir10->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
    dir13->fh_HoDTotTime->Fill(Event_HoDTime,w_tot);


    dir13->fh_muee_M->Fill(dir10->GetThreeTrackMomentum().M(),w_tot);
    dir13->fh_Pmu_vs_z->Fill(muon.GetMomentum(), dir10->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
    dir13->fh_Kaon_Charge->Fill(Kcharge );
    dir13->fh_eop->Fill(muon.GetEnergyLeftInEcal()/muon.GetMomentum(),w_tot);
    dir13->fh_eop->Fill(electron1.GetEnergyLeftInEcal()/electron1.GetMomentum(),w_tot);
    dir13->fh_eop->Fill(electron2.GetEnergyLeftInEcal()/electron2.GetMomentum(),w_tot);
    dir13->fh_Muee_3pi->Fill( three_pi.M(),w_tot);
    dir13->fh_M3pi_vs_EoP->Fill( three_pi.M(), electron2.GetEnergyLeftInEcal()/electron2.GetMomentum(),w_tot);
    dir13->fh_M3pi_vs_EoP->Fill( three_pi.M(), electron1.GetEnergyLeftInEcal()/electron1.GetMomentum(),w_tot);
    dir13->fh_munu_M->Fill(munu.M(),w_tot);
    dir13->fh_e1nu_M->Fill(e1nu.M(),w_tot);
    dir13->fh_e2nu_M->Fill(e2nu.M(),w_tot);
    dir13->fh_mc_two_track_23_mass_z_variable->Fill(mc_ee.M()*mc_ee.M()/(0.493677*0.493677));

    if(Kcharge>0){
        dir13->fh_mee_z_Kplus->Fill(dir10->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
        dir13->fh_MM2_plus->Fill(dir10->GetNuMomentum().M2(),w_tot );

    }

    if(Kcharge<0){
        dir13->fh_MM2_minus->Fill(dir10->GetNuMomentum().M2(),w_tot );
        dir13->fh_mee_z_Kminus->Fill(dir10->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
    }

    if(magnetCurrent < 0){
        dir13->fh_mee_z_magnet_minus->Fill(dir10->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
    }
    if(magnetCurrent > 0){
        dir13->fh_mee_z_magnet_plus->Fill(dir10->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
    }


    // Pt > 10 MeV/c
    if(dir10->GetThreeTrackMomentum().Pt() > 0.010){

        dir15->ComputeThreeTrack(electron1,electron2,muon);
        dir15->fh_missing_mass->Fill(dir10->GetNuMomentum().M2() ,w_tot);
        dir15->fh_muee_Pt->Fill(dir10->GetThreeTrackMomentum().Pt(),w_tot);
        dir15->fh_COmPaCt_Z_Vertex->Fill(COmPaCt_Z_Vertex,w_tot);
        dir15->fh_Mu_momentum->Fill(muon.GetMomentum(),w_tot);
        dir15->fh_mee->Fill(dir10->GetTwoTrackMomentum().M(),w_tot);
        dir15->fh_mee_z_variable->Fill(dir10->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
        dir15->fh_HoDTotTime->Fill(Event_HoDTime,w_tot);

        dir15->fh_muee_M->Fill(dir10->GetThreeTrackMomentum().M(),w_tot);
        dir15->fh_Pmu_vs_z->Fill(muon.GetMomentum(), dir10->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
        dir15->fh_Kaon_Charge->Fill(Kcharge );
        dir15->fh_eop->Fill(muon.GetEnergyLeftInEcal()/muon.GetMomentum(),w_tot);
        dir15->fh_eop->Fill(electron1.GetEnergyLeftInEcal()/electron1.GetMomentum(),w_tot);
        dir15->fh_eop->Fill(electron2.GetEnergyLeftInEcal()/electron2.GetMomentum(),w_tot);
        dir15->fh_Muee_3pi->Fill( three_pi.M(),w_tot);
        dir15->fh_M3pi_vs_EoP->Fill( three_pi.M(), electron2.GetEnergyLeftInEcal()/electron2.GetMomentum(),w_tot);
        dir15->fh_M3pi_vs_EoP->Fill( three_pi.M(), electron1.GetEnergyLeftInEcal()/electron1.GetMomentum(),w_tot);
        dir15->fh_munu_M->Fill(munu.M(),w_tot);
        dir15->fh_e1nu_M->Fill(e1nu.M(),w_tot);
        dir15->fh_e2nu_M->Fill(e2nu.M(),w_tot);
        dir15->fh_mc_two_track_23_mass_z_variable->Fill(mc_ee.M()*mc_ee.M()/(0.493677*0.493677));

        if(Kcharge>0){
            dir15->fh_mee_z_Kplus->Fill(dir10->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
            dir15->fh_MM2_plus->Fill(dir10->GetNuMomentum().M2(),w_tot );
        }

        if(Kcharge<0){
            dir15->fh_mee_z_Kminus->Fill(dir10->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
            dir15->fh_MM2_minus->Fill(dir10->GetNuMomentum().M2(),w_tot );
        }

        if(magnetCurrent < 0){
            dir15->fh_mee_z_magnet_minus->Fill(dir10->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
        }
        if(magnetCurrent > 0){
            dir15->fh_mee_z_magnet_plus->Fill(dir10->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
        }

    }

    // Pt > 20 MeV/c
    if(dir10->GetThreeTrackMomentum().Pt() > 0.020){

        dir16->ComputeThreeTrack(electron1,electron2,muon);
        dir16->fh_missing_mass->Fill(dir10->GetNuMomentum().M2() ,w_tot);
        dir16->fh_muee_Pt->Fill(dir10->GetThreeTrackMomentum().Pt(),w_tot);
        dir16->fh_COmPaCt_Z_Vertex->Fill(COmPaCt_Z_Vertex,w_tot);
        dir16->fh_Mu_momentum->Fill(muon.GetMomentum(),w_tot);
        dir16->fh_mee->Fill(dir10->GetTwoTrackMomentum().M(),w_tot);
        dir16->fh_mee_z_variable->Fill(dir10->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
        dir16->fh_HoDTotTime->Fill(Event_HoDTime,w_tot);

        dir16->fh_muee_M->Fill(dir10->GetThreeTrackMomentum().M(),w_tot);
        dir16->fh_Pmu_vs_z->Fill(muon.GetMomentum(), dir10->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
        dir16->fh_Kaon_Charge->Fill(Kcharge );
        dir16->fh_eop->Fill(muon.GetEnergyLeftInEcal()/muon.GetMomentum(),w_tot);
        dir16->fh_eop->Fill(electron1.GetEnergyLeftInEcal()/electron1.GetMomentum(),w_tot);
        dir16->fh_eop->Fill(electron2.GetEnergyLeftInEcal()/electron2.GetMomentum(),w_tot);
        dir16->fh_Muee_3pi->Fill( three_pi.M(),w_tot);
        dir16->fh_M3pi_vs_EoP->Fill( three_pi.M(), electron2.GetEnergyLeftInEcal()/electron2.GetMomentum(),w_tot);
        dir16->fh_M3pi_vs_EoP->Fill( three_pi.M(), electron1.GetEnergyLeftInEcal()/electron1.GetMomentum(),w_tot);
        dir16->fh_munu_M->Fill(munu.M(),w_tot);
        dir16->fh_e1nu_M->Fill(e1nu.M(),w_tot);
        dir16->fh_e2nu_M->Fill(e2nu.M(),w_tot);
        dir16->fh_mc_two_track_23_mass_z_variable->Fill(mc_ee.M()*mc_ee.M()/(0.493677*0.493677));

        if(Kcharge>0){
            dir16->fh_mee_z_Kplus->Fill(dir10->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
            dir16->fh_MM2_plus->Fill(dir10->GetNuMomentum().M2(),w_tot );
        }

        if(Kcharge<0){
            dir16->fh_mee_z_Kminus->Fill(dir10->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
            dir16->fh_MM2_minus->Fill(dir10->GetNuMomentum().M2(),w_tot );
        }

        if(magnetCurrent < 0){
            dir16->fh_mee_z_magnet_minus->Fill(dir10->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
        }
        if(magnetCurrent > 0){
            dir16->fh_mee_z_magnet_plus->Fill(dir10->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
        }

    }

    // Pt > 30 MeV/c
    if(dir10->GetThreeTrackMomentum().Pt() > 0.030){

        dir17->ComputeThreeTrack(electron1,electron2,muon);
        dir17->fh_missing_mass->Fill(dir10->GetNuMomentum().M2() ,w_tot);
        dir17->fh_muee_Pt->Fill(dir10->GetThreeTrackMomentum().Pt(),w_tot);
        dir17->fh_COmPaCt_Z_Vertex->Fill(COmPaCt_Z_Vertex,w_tot);
        dir17->fh_Mu_momentum->Fill(muon.GetMomentum(),w_tot);
        dir17->fh_mee->Fill(dir10->GetTwoTrackMomentum().M(),w_tot);
        dir17->fh_mee_z_variable->Fill(dir10->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
        dir17->fh_HoDTotTime->Fill(Event_HoDTime,w_tot);

        dir17->fh_muee_M->Fill(dir10->GetThreeTrackMomentum().M(),w_tot);
        dir17->fh_Pmu_vs_z->Fill(muon.GetMomentum(), dir10->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
        dir17->fh_Kaon_Charge->Fill(Kcharge );
        dir17->fh_eop->Fill(muon.GetEnergyLeftInEcal()/muon.GetMomentum(),w_tot);
        dir17->fh_eop->Fill(electron1.GetEnergyLeftInEcal()/electron1.GetMomentum(),w_tot);
        dir17->fh_eop->Fill(electron2.GetEnergyLeftInEcal()/electron2.GetMomentum(),w_tot);
        dir17->fh_Muee_3pi->Fill( three_pi.M(),w_tot);
        dir17->fh_M3pi_vs_EoP->Fill( three_pi.M(), electron2.GetEnergyLeftInEcal()/electron2.GetMomentum(),w_tot);
        dir17->fh_M3pi_vs_EoP->Fill( three_pi.M(), electron1.GetEnergyLeftInEcal()/electron1.GetMomentum(),w_tot);
        dir17->fh_munu_M->Fill(munu.M(),w_tot);
        dir17->fh_e1nu_M->Fill(e1nu.M(),w_tot);
        dir17->fh_e2nu_M->Fill(e2nu.M(),w_tot);
        dir17->fh_mc_two_track_23_mass_z_variable->Fill(mc_ee.M()*mc_ee.M()/(0.493677*0.493677));

        if(Kcharge>0){
            dir17->fh_mee_z_Kplus->Fill(dir10->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
            dir17->fh_MM2_plus->Fill(dir10->GetNuMomentum().M2(),w_tot );
        }

        if(Kcharge<0){
            dir17->fh_mee_z_Kminus->Fill(dir10->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
            dir17->fh_MM2_minus->Fill(dir10->GetNuMomentum().M2(),w_tot );
        }

        if(magnetCurrent < 0){
            dir17->fh_mee_z_magnet_minus->Fill(dir10->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
        }
        if(magnetCurrent > 0){
            dir17->fh_mee_z_magnet_plus->Fill(dir10->GetTwoTrackMomentum().M2()/(massKaonC*massKaonC),w_tot);
        }

    }







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
    cutting.cluster_e1e2 = electron1.GetClusterTime() - electron2.GetClusterTime();
    //cutting.cluster_mue1 = muon.GetClusterTime() - electron1.GetClusterTime();
    //cutting.cluster_mue2 = muon.GetClusterTime() - electron2.GetClusterTime();
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

        cutting.Lkr_x_el1 = electron1.extrapolated_track_Lkr[0];
        cutting.Lkr_x_el2 = electron2.extrapolated_track_Lkr[0];
        cutting.Lkr_x_mu  = muon.extrapolated_track_Lkr[0];
        cutting.Lkr_y_el1 = electron1.extrapolated_track_Lkr[1];
        cutting.Lkr_y_el2 = electron2.extrapolated_track_Lkr[1];
        cutting.Lkr_y_mu  = muon.extrapolated_track_Lkr[1];
        cutting.Lkr_cut_el1 = electron1.GetLkrRadius();
        cutting.Lkr_cut_el2 = electron2.GetLkrRadius();
        cutting.Lkr_cut_mu  = muon.GetLkrRadius();

        cutting.MUV_x_el1  = electron1.extrapolated_track_MUV2[0];
        cutting.MUV_x_el2  = electron2.extrapolated_track_MUV2[0];
        cutting.MUV_x_mu   = muon.extrapolated_track_MUV2[0]     ;
        cutting.MUV_y_el1  = electron1.extrapolated_track_MUV2[1];
        cutting.MUV_y_el2  = electron1.extrapolated_track_MUV2[1];
        cutting.MUV_y_mu   = muon.extrapolated_track_MUV2[1]     ;


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

        cutting.Lkr_x_pi1 = muon.extrapolated_track_Lkr[0];
        cutting.Lkr_x_pi2 = electron1.extrapolated_track_Lkr[0];
        cutting.Lkr_x_pi3 = electron2.extrapolated_track_Lkr[0];
        cutting.Lkr_y_pi1 = muon.extrapolated_track_Lkr[1];
        cutting.Lkr_y_pi2 = electron1.extrapolated_track_Lkr[1];
        cutting.Lkr_y_pi3 = electron2.extrapolated_track_Lkr[1];
        cutting.Lkr_cut_pi1 = muon.GetLkrRadius();
        cutting.Lkr_cut_pi2 = electron1.GetLkrRadius();
        cutting.Lkr_cut_pi3 = electron2.GetLkrRadius();


        cutting.MUV_x_pi1  = muon.extrapolated_track_MUV2[0]     ;
        cutting.MUV_x_pi2  = electron1.extrapolated_track_MUV2[0];
        cutting.MUV_x_pi3  = electron2.extrapolated_track_MUV2[0];
        cutting.MUV_y_pi1  = muon.extrapolated_track_MUV2[1]     ;
        cutting.MUV_y_pi2  = electron1.extrapolated_track_MUV2[1];
        cutting.MUV_y_pi3  = electron2.extrapolated_track_MUV2[1];

    }

    return cutting;
}

void FillMC(Hist_dir* dir, TLorentzVector p1, TLorentzVector p2, TLorentzVector p3, float* DKaon, float* part_production, float* part_decay){
    TLorentzVector mc_Two_Track_Momentum;
    TLorentzVector mc_Three_Track_Momentum;
    mc_Three_Track_Momentum = p1 + p2 + p3;
    mc_Two_Track_Momentum   = p2 + p3;

    //NOT wORKING ///////////////////////////////////////////
    double P1_distance_prod_dec = DKaon[2] - part_production[1];
    dir->fh_mc_P1_dist_prod_dec->Fill(P1_distance_prod_dec);
    double P2_distance_prod_dec = DKaon[2] - part_production[2];
    dir->fh_mc_P2_dist_prod_dec->Fill(P2_distance_prod_dec);
    double P3_distance_prod_dec = DKaon[2] - part_production[3];
    dir->fh_mc_P3_dist_prod_dec->Fill(P3_distance_prod_dec);
    ///////////////////////////////////////////////////////
    dir->fh_mc_KDzvtx->Fill(DKaon[2]);
    dir->fh_mc_P1_mass->Fill(p1.M());
    dir->fh_mc_P1_momentum->Fill(p1.P());
    dir->fh_mc_P1_Pzvtx->Fill(part_production[1]);
    dir->fh_mc_P2_mass->Fill(p2.M());
    dir->fh_mc_P2_momentum->Fill(p2.P());
    dir->fh_mc_P2_Pzvtx->Fill(part_production[2]);
    dir->fh_mc_P3_mass->Fill(p3.M());
    dir->fh_mc_P3_momentum->Fill(p3.P());
    dir->fh_mc_P3_Pzvtx->Fill(part_production[3]);
    dir->fh_mc_three_track_123_momentum->Fill(mc_Three_Track_Momentum.P());
    dir->fh_mc_three_track_123_mass->Fill(mc_Three_Track_Momentum.M());
    dir->fh_mc_two_track_23_momentum->Fill(mc_Two_Track_Momentum.P());
    dir->fh_mc_two_track_23_mass->Fill(mc_Two_Track_Momentum.M());
    dir->fh_mc_two_track_23_mass_z_variable->Fill(mc_Two_Track_Momentum.M()*mc_Two_Track_Momentum.M()/(0.493677*0.493677));
    dir->fh_Pmu_vs_z->Fill( p1.P(),mc_Two_Track_Momentum.M()*mc_Two_Track_Momentum.M()/(0.493677*0.493677) );

}
void FillMC(Hist_dir* dir, TLorentzVector p1, TLorentzVector p2, TLorentzVector p3, TLorentzVector p4, float* DKaon, float* part_production, float* part_decay){
    TLorentzVector mc_Four_Track_Momentum;

    mc_Four_Track_Momentum = p1 + p2 + p3 + p4;
    dir->fh_mc_P4_Pzvtx->Fill(part_production[4]);
    dir->fh_mc_P4_mass->Fill(p4.M());
    dir->fh_mc_P4_momentum->Fill(p4.P());
    dir->fh_mc_four_track_1234_mass->Fill(mc_Four_Track_Momentum.M());
    dir->fh_mc_four_track_1234_momentum->Fill(mc_Four_Track_Momentum.P());


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
