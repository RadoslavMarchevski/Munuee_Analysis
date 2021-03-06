/* COmPACT user routine: user_superCmpEvent(superCmpEvent *sevt) */

#include <math.h>
#include "cmpio.h"
#include "user.h"
#include "reader.h"
#include <constants.h>

#include <TAxis.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TDirectory.h>
#include <TVector3.h>

#include <iostream>
#include <vector>
#include <math.h>
using namespace std;


vector<TVector3*> trkCorrSlopesV;
vector<TVector3*> trkCorrMidPointsV;
TVector3* trkCorrSlope     = NULL;
TVector3* trkCorrMidPoint     = NULL;

vector<TVector3*> clcorrpos;
TVector3* clpos = NULL;

TLorentzVector *el1lv=NULL, *el2lv=NULL, *eelv=NULL;
TLorentzVector *pi1lv=NULL, *pi2lv=NULL,*pi3lv=NULL;
TLorentzVector *pchlv=NULL, *muchlv=NULL , *muvee=NULL;
TLorentzVector *pK=NULL, *neutrino=NULL,*k3pi_miss_mass=NULL ;
TLorentzVector *pieelv=NULL,*k3pic=NULL;
TLorentzVector *glv =NULL;
TLorentzVector *pi0lv=NULL;
TLorentzVector *el1lvcm=NULL, *el2lvcm=NULL,  *glvcm=NULL, *pi0lvcm=NULL;
TLorentzVector *Klv=NULL, *pi0lvKcm=NULL, *pchlvKcm=NULL, *muchlvKcm=NULL, *KlvKcm=NULL,*muveeKcm=NULL;

TLorentzVector *v_first_e1;
TLorentzVector *v_second_e1;
TLorentzVector *v_first_e2;
TLorentzVector *v_second_e2;
TLorentzVector *v_first_mu;
TLorentzVector *v_second_mu;


void cleanup();

//int thetaFunc(int t){
//if(t>0) return 1;
//return 0;
//}

void applyBlueTubeCorrection(superCmpEvent *sevt, int itrk, double pKaon[3], double vKaon[3], double* trkSlopes, double* trkMidPoints);

int user_superCmpEvent(superBurst *sbur,superCmpEvent *sevt) {
  /* WARNING: do not alter things before this line */
  /*---------- Add user C code here ----------*/
  //cout << "Start" << endl;
  static int nuserevt=0;
  int  a, i, j, k;
  int  l, m, n;
  //cout << "Before" << endl;  
  ////thetaFunc(sevt->ntrack);
  ////return 0;
  //cout << "After" << endl;
 
  if(nuserevt<1)
    {
      //gkh      printSuperCmpEvent(sevt,fprt);
    }

  nuserevt++;
  //printf ("testVariable = %f\n", testVariable);
  //testVariable++;

  // definitions of particle masses (from PDG 2006)
  static double massElec = 0.000510998918;
  static double massMuonC = 0.105658369;
  static double massPionC = 0.13957018;
  static double massKaonC = 0.493677;
  //##########################################################
  // End of Cut definitions
  //#########################################################

  float  beamCorrX, beamOffsetX, beamCorrY, beamOffsetY;
  int    Kcharge=0;
  double vKaon[3] = {0.};
  double pKaon[3] = {0.};

  double clenergy=0, clx=0, cly=0, clz=0;
  //    vector<double> trkP;
  //    double trkcorrP=0;
  double eop_track;             // E/p for the track
  double p_track, e_track;      // momentum and energy of the track
  double eop_mu_track;             //E/p of the selected muon
  double eop_e1_track,eop_e2_track;//E/p of the selected electrons
  double eop_pi1_track,eop_pi2_track,eop_pi3_track;
  double ddcell_track;             // distance to the dead cell
  int    icl_track;            // cluster index associated to the track
  int    imu_track;            // muon index associated to the track
  int    clst_track=999;            // cluster status
  float  quality_track;
  float  trkxcl, trkycl, dtrkcl;
  int    nelectrack=0;
  int    npitrack=0;
  int    nmutrack=0;
  int    nphcluster=0;
  int    evttypebytrks=-1;
  int    k3pi_evttypebytrks=-1;
  double el1trkhodtime=-1e+3;
  double el2trkhodtime=-1e+3;
  double el1trktime=-1e+3;
  double el2trktime=-1e+3;
  double pi1trkhodtime=-1e+3;
  double pi2trkhodtime=-1e+3;
  double pi3trkhodtime=-1e+3;
  double pi1trktime=-1e+3;
  double pi2trktime=-1e+3;
  double pi3trktime=-1e+3;

  double norm1, norm2,normmu;
  double p1x,p1y, p1z, p2x,p2y, p2z, el1E, el2E, pK1x,pK1y, pK1z,pKE,pmux,pmuy,pmuz;
  double normpi1, normpi2,normpi3;
  double pi1x,pi1y, pi1z, pi2x,pi2y, pi2z, pi1E, pi2E, pi3x,pi3y, pi3z,pi3E;
   
  int    gclindx=-1;
  double gclenergy=0.;
  double gcltime=0.;
  double gclx=0.;
  double gcly=0.;
  double gclz=0.;
  double gtrkDHC1x, gtrkDHC1y, gtrkRadDHC1 ;
  
  
  double pchtrkhodtime, pchE , muchE , muchtrkhodtime;
  
  double eldistDCH1plane, el1xLKrplane, el1yLKrplane, el2xLKrplane, el2yLKrplane, eldistLKrplane;
  //  double pidistDCH1plane, pi1xLKrplane, pi1yLKrplane, pi2xLKrplane, pi2yLKrplane, pidistLKrplane;
  
  double gel1distLKrplane, gel2distLKrplane;
  double pche1distDHC1plane, pche2distDHC1plane, pche1distLKrplane, pche2distLKrplane;
  double muche1distDHC1plane, muche2distDHC1plane, muche1distLKrplane, muche2distLKrplane;
    
  double pchxLKrplane, pchyLKrplane, gpchdistLKrplane;
  double muchxLKrplane, muchyLKrplane, gmuchdistLKrplane;
  double Dgvtx_x,Dgvtx_y, Dgvtx_z, gnorm,  gpx, gpy, gpz;

  unsigned int vertex_OK ,k3pic12_vertex_OK , k3pic23_vertex_OK, k3pic13_vertex_OK ;
  double vertex[3], vertex_k3pic12[3], vertex_k3pic23[3],vertex_k3pic13[3];
  double cda, cda_k3pic12, cda_k3pic23, cda_k3pic13;

  unsigned int muchee_vertex_OK,muchee2_vertex_OK;

  double vertex_e2much[3], vertex_e1much[3];
  double cda_e2much, cda_e1much;

  double elp1_uncorr[3], elp1[3], elv1_uncorr[3], elv1[3], elp2_uncorr[3],elp2[3], elv2_uncorr[3],elv2[3], pchp_uncorr[3],pchp[3], pchv_uncorr[3], pchv[3] ,muchp_uncorr[3],muchp[3], muchv_uncorr[3], muchv[3];
  double pilp1_uncorr[3], pilp1[3], pilv1_uncorr[3], pilv1[3], pilp2_uncorr[3],pilp2[3], pilv2_uncorr[3],pilv2[3],pilp3_uncorr[3], pilp3[3], pilv3_uncorr[3], pilv3[3];
  
  int    pchtrkindx=-1;
  double pchtrktime=-1e+3;
  int    pchtrkq=0;
  double pchcltime=-1e+3;
  double pcheopopt=0;
  double pcheopposition=0;

  int    muchtrkindx=-1;
  double muchtrktime=-1e+3;
  int    muchtrkq=0;
  double muchcltime=-1e+3;
  double mucheopopt=0;
  double mucheopposition=0;
    
  int    eltrkindx1=-1;
  int    eltrkindx2=-1;
  int    eltrkq1=0;
  int    eltrkq2=0;
  int    pi1indx=-1;
  int    pi2indx=-1;
  int    pi3indx=-1;
  int    pitrkq1=0;
  int    pitrkq2=0;
  int    pitrkq3=0;
  double pi1cltime=-1e+3;
  double pi2cltime=-1e+3;
  double pi3cltime=-1e+3;
  
  double eleop1=0;
  double eleop2=0;
  double el1cltime=-1e+3;
  double el2cltime=-1e+3;
  double el1P3;
  double el2P3;

  double eelvmass;

  double c_eop_e=0.9;//0.95;
  double c_eop_e_up=1.05;
  double c_eop_mu_up=0.1;
  double c_eop_pi_min=0.2;
  double c_eop_pi_max=0.85;
  
  int flux_id =0;
  int f_cout=0;
  
  bool   LKrCalCorr=true;
  bool   CorrAlphaBeta=true;
  bool Ke2L3trig=false;
  bool Km2L3trig=false;
  
 
    
  int       cutcounter=1;
  int    ntrack = sevt->Ntrack; // number of tracks
  int    nclust = sevt->Ncluster; // number of tracks

  double beta=0.;
  double DCHbz = Geom->DCH.bz;         // z before magnet
  double DCHz = Geom->DCH.z;         // z after magnet
  double LKrz=Geom->Lkr.z;

  int       ngoodtrack=0;


  // --------------------------------------------------------------------------------------------->
  // Kaon beam coordinates
  // Kaon flight direction along the beam axis
  //
  for (j=0; j<ntrack; j++) {
    Kcharge += sevt->track[j].q;
  }

  //        if (SQL_Database)    {
  if (Kcharge >= 0)                // K+ beam
    {
      beamCorrX     = abcog_params.pkdxdzp;
      beamOffsetX = abcog_params.pkxoffp;
      beamCorrY     = abcog_params.pkdydzp;
      beamOffsetY = abcog_params.pkyoffp;
    }
  else
    {
      beamCorrX     = abcog_params.pkdxdzm;
      beamOffsetX = abcog_params.pkxoffm;
      beamCorrY     = abcog_params.pkdydzm;
      beamOffsetY = abcog_params.pkyoffm;
    }
  //        }

  vKaon[0] = beamCorrX;
  vKaon[1] = beamCorrY;
  vKaon[2] = 1;
  pKaon[0] = beamOffsetX + beamCorrX*DCHbz;
  pKaon[1] = beamOffsetY + beamCorrY*DCHbz;
  pKaon[2] = DCHbz;
  // <----------------------------------------------------------------------------------------------



  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------
  //                            General corrections for objects
  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------


  /// apply non linearity corrections
  /// Correct LKr non-linearity (effective for energies < 11GeV)
  if (LKrCalCorr && IS_DATA) {
    user_lkrcalcor_SC (sbur,sevt,1);
  }


  //apply projectivity corrections
  for (int i = 0; i < nclust; ++i) {

    clenergy = sevt->cluster[i].energy;
    clx= sevt->cluster[i].x;
    cly= sevt->cluster[i].y;
    double x0 = clx;
    double y0 = cly;

    clz = LKrz + 16.5 + 4.3*log(clenergy);

    if(IS_DATA)    {
      clx = (x0 + 0.0136 + y0 * 0.87e-3 ) * ( 1 + (clz-LKrz)/10998.);
      cly = (y0 + 0.300  + x0 * 0.87e-3 ) * ( 1 + (clz-LKrz)/10998.);
    }
    else    {
      clx = (x0 - 0.013 ) * ( 1 + (clz-LKrz)/10998.);
      cly = y0              * ( 1 + (clz-LKrz)/10998.);
    }

    clpos = new TVector3(clx, cly, clz);
    clcorrpos.push_back(clpos);
  }
 

  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------
  //                                        object selection
  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------
  

  // ~~~ NO CUT ~~~ NO CUT ~~~ NO CUT ~~~
  // ************* Initial cuts ****************** ///
  dirs[cutcounter-1]->cd();
  
  if (IS_DATA){
    ((TH1D*)gDirectory->FindObject( h_vz_initial->GetName()))->Fill( sevt->vtx[0].z);
  }
    
  
 
  if(sevt->vtx[0].z < -2000 || sevt->vtx[0].z > 8000) {cleanup();return 0;}
     
  h_ntrack->Fill(ntrack);
  h_nclust->Fill(nclust);
  
  
  if (IS_MC){
    ////if(true_4.P()  < 3. || true_5.P() < 3.) {cleanup();return 0;}
    
    ((TH1D*)gDirectory->FindObject( h_mc_muvee_M->GetName()))->Fill( Muee_three_track.M());
    ((TH1D*)gDirectory->FindObject( h_mc_muvee_P->GetName()))->Fill( Muee_three_track.P());
    ((TH1D*)gDirectory->FindObject( h_mc_muvee_Pt->GetName()))->Fill( Muee_three_track.Pt());
    ((TH1D*)gDirectory->FindObject( h_mc_elsP->GetName()))->Fill( mc_ee.P());
    ((TH1D*)gDirectory->FindObject( h_mc_elsM->GetName()))->Fill( mc_ee.M());
  
    ((TH1D*)gDirectory->FindObject( h_mc_nu_M->GetName()))->Fill( Nu_true.M());
    ((TH1D*)gDirectory->FindObject( h_mc_nu_M2->GetName()))->Fill( Nu_true.M2());
    ((TH1D*)gDirectory->FindObject( h_mc_nu_P->GetName()))->Fill( Nu_true.P());
    ((TH1D*)gDirectory->FindObject( h_mc_nu_Pt->GetName()))->Fill(Nu_true.Pt());
    
    
    
    ((TH1D*)gDirectory->FindObject(h_true_M_0->GetName()))->Fill( true_0.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_0->GetName()))->Fill( true_0.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_0->GetName()))->Fill(true_0.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_0->GetName()))->Fill( true_0.E());
  
    ((TH1D*)gDirectory->FindObject(h_true_M_1->GetName()))->Fill( true_1.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_1->GetName()))->Fill( true_1.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_1->GetName()))->Fill( true_1.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_1->GetName()))->Fill( true_1.E());
  
    ((TH1D*)gDirectory->FindObject(h_true_M_2->GetName()))->Fill( true_2.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_2->GetName()))->Fill( true_2.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_2->GetName()))->Fill(true_2.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_2->GetName()))->Fill( true_2.E());
  
    ((TH1D*)gDirectory->FindObject(h_true_M_3->GetName()))->Fill( true_3.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_3->GetName()))->Fill( true_3.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_3->GetName()))->Fill(true_3.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_3->GetName()))->Fill( true_3.E());
  
    ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_4.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_4.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_4.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_4.E());
  
    ((TH1D*)gDirectory->FindObject(h_true_M_5->GetName()))->Fill( true_5.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_5->GetName()))->Fill( true_5.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_5->GetName()))->Fill(true_5.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_5->GetName()))->Fill( true_5.E());
  
    ((TH1D*)gDirectory->FindObject(h_true_M_6->GetName()))->Fill( true_6.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_6->GetName()))->Fill( true_6.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_6->GetName()))->Fill(true_6.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_6->GetName()))->Fill( true_6.E());
  
    ((TH1D*)gDirectory->FindObject(h_true_M_7->GetName()))->Fill( true_7.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_7->GetName()))->Fill( true_7.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_7->GetName()))->Fill(true_7.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_7->GetName()))->Fill( true_7.E());
  
  }
  
  /// ************* tracks  &  clusters  *******************************
  
  // ~~~ CUT ~~~ CUT ~~~ CUT ~~~
  //if (flux_id ==1 ){
  //Ke2L3trig=sevt->DETstatus[0].LV3ABTrig & 1 ;
  //Km2L3trig=sevt->DETstatus[0].LV3ABTrig & 2 ;
  //if (Km2L3trig)
  //  cout << Km2L3trig << endl;
  
  if(ntrack!=3 ) {cleanup();return 0;}
    
  //if(nclust< 2  ) {cleanup();return 0;}
  //}
  //"3tracks"
  
  
  for (int i=0; i<ntrack; i++) {
    
    // Alpha + Beta Correction --> call p_corr_ab() routine
    if (CorrAlphaBeta)
      {
	//            cout<<"Not corr. p="<<sevt->track[i].p;
	//            beta = abcog_params.beta; //for Kaons
	sevt->track[i].p = p_corr_ab(sevt->track[i].p, sevt->track[i].q);
	//            cout<<"        Corr. p="<<sevt->track[i].p<<endl;
      }
    
    p_track            = sevt->track[i].p;
    quality_track     = sevt->track[i].quality;
    ddcell_track    = sevt->track[i].dDeadCell;
    icl_track         = sevt->track[i].iClus;
    imu_track         = sevt->track[i].iMuon;
    e_track  = sevt->cluster[icl_track].energy;
    eop_track  = e_track / p_track;

    if(icl_track>-1)
      clst_track     = sevt->cluster[icl_track].status;

    //        if(imu_track>-1)
    //            cout<<i<<"    imu_track="<<imu_track<<"    "<<sevt->muon[imu_track].chi2<<endl;
   
    if (quality_track<0.8)        { trkCorrSlopesV.push_back(NULL); trkCorrMidPointsV.push_back(NULL); continue; } //0.8 is in Manuel's thesis
    if (icl_track<0)             { trkCorrSlopesV.push_back(NULL); trkCorrMidPointsV.push_back(NULL); continue; }
    //if (imu_track>-1)             { trkCorrSlopesV.push_back(NULL); trkCorrMidPointsV.push_back(NULL); continue; }
    if(IS_DATA)
      if (ddcell_track<2.)        { trkCorrSlopesV.push_back(NULL); trkCorrMidPointsV.push_back(NULL); continue; }
    // if(flux_id==1){
    //  if (p_track<5.)                { trkCorrSlopesV.push_back(NULL); trkCorrMidPointsV.push_back(NULL); continue; }
    //  if (p_track>50.)            { trkCorrSlopesV.push_back(NULL); trkCorrMidPointsV.push_back(NULL); continue; }
    // }
    //if (flux_id==0){
    if (p_track<3.)                { trkCorrSlopesV.push_back(NULL); trkCorrMidPointsV.push_back(NULL); continue; }
    if (p_track>75.)            { trkCorrSlopesV.push_back(NULL); trkCorrMidPointsV.push_back(NULL); continue;// }
    }
    if(IS_DATA)
      if (clst_track>3)            { trkCorrSlopesV.push_back(NULL); trkCorrMidPointsV.push_back(NULL); continue; }
    
        
    //calculate distance between the track and its matching cluster in LKr
    gclx=clcorrpos[icl_track]->x();
    gcly=clcorrpos[icl_track]->y();
    gclz=clcorrpos[icl_track]->z();
    trkxcl = sevt->track[i].x + (gclz-DCHz) * sevt->track[i].dxdz ;
    trkycl = sevt->track[i].y + (gclz-DCHz) * sevt->track[i].dydz ;
    dtrkcl = sqrt( pow(trkxcl-gclx, 2.) + pow(trkycl-gcly, 2.) );

    //((TH1F*)gDirectory->FindObject(h_dtrkcl->GetName()))        ->Fill( dtrkcl );
    
    //if (IS_DATA)
    //  if(dtrkcl>6.)                 { trkCorrSlopesV.push_back(NULL); trkCorrMidPointsV.push_back(NULL); continue; }

    //number of good tracks , cut variable
    ngoodtrack++;
    //cout << "you are stupid !!!!!!!!!!!!!"  << IS_DATA << endl;
    //
    //calculate Blue tube corrected slopes and track middle points before these corrections
    double trkSlopes[3]     = {0.};
    double trkMidPoints[3]     = {0.};
   
    applyBlueTubeCorrection(sevt, i, pKaon, vKaon, trkSlopes, trkMidPoints);
    
    //cout << "trkSlopes =" << trkSlopes[0] << "=1=  "<< trkSlopes[1]<< "=2=  " << trkSlopes[2] << "=3=  " << endl;
    //cout << "trkmidpoints =" << trkMidPoints[0] << "=1=  "<< trkMidPoints[1]<< "=2=  " << trkMidPoints[2] << "=3=  " << endl;
    trkCorrSlope    = new TVector3(trkSlopes[0],     trkSlopes[1],     trkSlopes[2]);
    trkCorrMidPoint = new TVector3(trkMidPoints[0], trkMidPoints[1],trkMidPoints[2]);
    trkCorrSlopesV.push_back(trkCorrSlope);
    trkCorrMidPointsV.push_back(trkCorrMidPoint);
    // cout << "kakvo stava tuka =" << trkCorrSlopesV[1]->x() << "=1=  "<< trkCorrSlopesV[2]<< "=2=  " << trkCorrSlopesV[0] << "=3=  " << trkCorrSlopesV[3] << "=4=  "  << "give me i =" << i << endl;
    //cout << "kakvo stava tuka =" << trkCorrMidPointsV[0] << "=1=  "<< trkCorrMidPointsV[1]<< "=2=  " << trkCorrMidPointsV[2] << "=3=  " << endl;

    if( p_track > 10. && eop_track < c_eop_e_up )    {
      if(pi1indx<0)
	//	 && sevt->track[i].q==-1)
	//	//for Piee only
	//	//p_track < 20
	{
	  pi1indx =i;
	  pitrkq1 = sevt->track[i].q;
     	  
	}
      else if(pi2indx<0 /*&& imu_track==-1*/)    {
	pi2indx=i;
	pitrkq2 = sevt->track[i].q;
      }
      //            cout<< "electron ";
      //WORK ON  making e+ , e-!!!!!!
      else  if( pi3indx   < 0 )    {
	pi3indx=i;
	pitrkq3=sevt->track[i].q;
	//nmutrack++;
      } else {
	;
      }
      //nelectrack++;
	
      npitrack++;
    
    }
    else    {
      ;//            cout<< "?? ";;
    }
     
 
    if(eop_track>c_eop_e  &&  eop_track<c_eop_e_up)    {


      if(eltrkindx1<0)
       	//	 && sevt->track[i].q==-1)
       	//	//for Piee only
       	//	//p_track < 20
       	{
       	  eltrkindx1=i;
       	  eltrkq1 = sevt->track[i].q;
       	  
       	}
      else if(eltrkindx2<0 /*&& imu_track==-1*/)    {
       	eltrkindx2=i;
       	eltrkq2 = sevt->track[i].q;
      }
      else { ; }
	
    

      //if (eltrkindx1<0 || eltrkindx2<0){
      //
      //
      //	if(sevt->track[i].q > 0){
      //	  if (eltrkindx2<0 ){
      //	    eltrkindx2=i;
      //	    eltrkq2 = sevt->track[i].q;
      //	  }
      //	  
      //	} else if (eltrkindx1<0 ){
      //	   
      //	  eltrkindx1=i;
      //	  eltrkq1 = sevt->track[i].q;
      //	}
      //	
      //
      //}else { ; }

      
      nelectrack++;
    }
    // muon track
    //Mahame za Piee MC 
    else if(eop_track < c_eop_mu_up)
      //else if(eop_track < c_eop_pi_max)
    
      
    
      //Piee
      //sevt->track[i].q*sevt->vtx[0].charge > 0 && 

      //eop_track < 0.95//Piee
      {
	if( muchtrkindx < 0 && imu_track !=-1)    {
	  muchtrkindx=i;
	  muchtrkq=sevt->track[i].q;
	}
    
	nmutrack++;
    
      }
    //track with too high e/p. skip it.
    else    {
      ;//            cout<< "?? ";;
    }

  }

  //        cout<<"1"<<endl;
  
  //file1->cd();
  //h_cutflow->Fill(cutcounter);
  //h_cutflow->GetXaxis()->SetBinLabel( cutcounter, cuts[cutcounter-1].c_str() );// "3tracks_4clusters");
  //cutcounter++;

  
  
  // ~~~ CUT ~~~ CUT ~~~ CUT ~~~
  //CHANGE
  //if(IS_DATA){
   
  if(ngoodtrack!=3) {cleanup();return 0;}
 
  // cout << "WTFFFF" << nmutrack << nelectrack << endl;
  //Munue+e- track assignment selection
  if(nelectrack+nmutrack==3)    {
    if(muchtrkq>0.)    {
      if(eltrkq1 != eltrkq2)    {// mu+, e+, e-
	evttypebytrks=0;
      }
      else if(eltrkq1>0)    {    // mu+, e+, e+
	evttypebytrks=1;
      }
      else    {                // mu+, e-, e-
	evttypebytrks=2;
      }
    }
    else    {
      if(eltrkq1 != eltrkq2)    {// mu-, e+, e-
	evttypebytrks=3;
      }
      else if(eltrkq1>0)    {    // mu-, e+, e+
	evttypebytrks=4;
      }
      else    {                // mu-, e-, e-
	evttypebytrks=5;
      }
    }
  }
     
  if(npitrack==3)    {
    if((pitrkq1+pitrkq2+pitrkq3)==3.)    {// pi+, pi+, pi+
      k3pi_evttypebytrks=0;
    } else if ((pitrkq1+pitrkq2+pitrkq3)== -3.){// pi-, pi-, pi-
      k3pi_evttypebytrks=1;
    } else if ((pitrkq1+pitrkq2+pitrkq3)== 1.){// pi+, pi+, pi- , pi-, pi+, pi+ ,  pi+, pi-, pi+
      k3pi_evttypebytrks=2;
    } else if ((pitrkq1+pitrkq2+pitrkq3)== -1.){// pi+, pi-, pi- , pi-, pi-, pi+ ,  pi-, pi+, pi-
      k3pi_evttypebytrks=3;
    }



  }
  ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_track);
  ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_track, p_track);
  ((TH1F*)gDirectory->FindObject(h_dtrkcl->GetName()))        ->Fill( dtrkcl );
    
  
  //Pion initial histograms
  ((TH1I*)gDirectory->FindObject(h_k3pi_evttypebytrks->GetName()))->Fill(k3pi_evttypebytrks);
  ((TH1I*)gDirectory->FindObject(h_k3pi_pitrack->GetName()))->Fill(npitrack);
  ((TH1I*)gDirectory->FindObject(h_k3pi_pi1trackq->GetName()))->Fill(pitrkq1);
  ((TH1I*)gDirectory->FindObject(h_k3pi_pi2trackq->GetName()))->Fill(pitrkq2);
  ((TH1I*)gDirectory->FindObject(h_k3pi_pi3trackq->GetName()))->Fill(pitrkq3);
  ///Pion initial histograms
  ((TH1I*)gDirectory->FindObject(h_evttypebytrks->GetName()))->Fill(evttypebytrks);
  ((TH1I*)gDirectory->FindObject(h_nelectrack->GetName()))->Fill(nelectrack);
  ((TH1I*)gDirectory->FindObject(h_mutrack->GetName()))->Fill(nmutrack);
  ((TH2I*)gDirectory->FindObject(h_neltrk_mutrk->GetName()))->Fill(nelectrack, nmutrack);
  ((TH1I*)gDirectory->FindObject(h_mutrackq->GetName()))->Fill(muchtrkq);
  ((TH1I*)gDirectory->FindObject(h_e1trackq->GetName()))->Fill(eltrkq1);
  ((TH1I*)gDirectory->FindObject(h_e2trackq->GetName()))->Fill(eltrkq2);

 
  // ~~~ CUT ~~~ CUT ~~~ CUT ~~~
  //if(IS_DATA)
  
  //if(dtrkcl>6.)        {cleanup();return 0;}

  //if(dtrkcl<200.)        {cleanup();return 0;}
  
 
  //"cluster_isolation"
  //printf("abcog_params.cogX1p=%f, abcog_params.cogX1n=%f\n",abcog_params.cogX1p,abcog_params.cogX1n);
  //printf("abcog_params.cogX4p=%f,abcog_params.cogX4n=%f\n",abcog_params.cogX4p,abcog_params.cogX4n);
  //printf("abcog_params.cogY1p=%f, abcog_params.cogY1n=%f\n",abcog_params.cogY1p,abcog_params.cogY1n);
  //printf("abcog_params.cogY4p=%f,abcog_params.cogY4n=%f\n",abcog_params.cogY4p,abcog_params.cogY4n);
  
  //cout << "---------BEGINNING OF EVENT ----------------"<< endl;
  //cout << "Event number =" << f_cout << "Mu index =" << muchtrkindx << "E1 index =" << eltrkindx1 << "E2 index =" << eltrkindx2 << endl;
  //cout << "Time stamp =" << sevt->timeStamp << endl;
  //f_cout++;
  //cout << "---------ENDING OF EVENT ----------------"<< endl;

  
  //  ------------------------------ K3PIC selection part ------------------------------------
  if( pi1indx != -1 &&  pi2indx != -1 && pi3indx !=-1){
    //cout << pi1indx << pi2indx << pi3indx << endl;
    //return 0;
    pilp1_uncorr[0] = sevt->track[pi1indx].bx;
    pilp1_uncorr[1] = sevt->track[pi1indx].by;
    pilp1_uncorr[2] = DCHbz;
    pilp2_uncorr[0] = sevt->track[pi2indx].bx;
    pilp2_uncorr[1] = sevt->track[pi2indx].by;
    pilp2_uncorr[2] = DCHbz;
    pilp3_uncorr[0] = sevt->track[pi3indx].bx;
    pilp3_uncorr[1] = sevt->track[pi3indx].by;
    pilp3_uncorr[2] = DCHbz;

    pilv1_uncorr[0] = sevt->track[pi1indx].bdxdz;
    pilv1_uncorr[1] = sevt->track[pi1indx].bdydz;
    pilv1_uncorr[2] = 1.;
    pilv2_uncorr[0] = sevt->track[pi2indx].bdxdz;
    pilv2_uncorr[1] = sevt->track[pi2indx].bdydz;
    pilv2_uncorr[2] = 1.;
    pilv3_uncorr[0] = sevt->track[pi3indx].bdxdz;
    pilv3_uncorr[1] = sevt->track[pi3indx].bdydz;
    pilv3_uncorr[2] = 1.;


    // if (IS_DATA){
    //cout << "Index breaks when eltrkindx1 =" << eltrkindx1 << " \t and eltrkindx2 = " << eltrkindx2 << endl;
    //blue field corrected values
    pilv1[0] = trkCorrSlopesV[pi1indx]->x();
    pilv1[1] = trkCorrSlopesV[pi1indx]->y();
    pilv1[2] = 1.;
    pilp1[0] = trkCorrMidPointsV[pi1indx]->x();
    pilp1[1] = trkCorrMidPointsV[pi1indx]->y();
    pilp1[2] = trkCorrMidPointsV[pi1indx]->z();

  
    pilv2[0] = trkCorrSlopesV[pi2indx]->x();
    pilv2[1] = trkCorrSlopesV[pi2indx]->y();
    pilv2[2] = 1.;
    pilp2[0] = trkCorrMidPointsV[pi2indx]->x();
    pilp2[1] = trkCorrMidPointsV[pi2indx]->y();
    pilp2[2] = trkCorrMidPointsV[pi2indx]->z();
   
    pilv3[0] = trkCorrSlopesV[pi3indx]->x();
    pilv3[1] = trkCorrSlopesV[pi3indx]->y();
    pilv3[2] = 1.;
    pilp3[0] = trkCorrMidPointsV[pi3indx]->x();
    pilp3[1] = trkCorrMidPointsV[pi3indx]->y();
    pilp3[2] = trkCorrMidPointsV[pi3indx]->z();

    //}
    normpi1 = 1./sqrt(pilv1_uncorr[0]*pilv1_uncorr[0] + pilv1_uncorr[1]*pilv1_uncorr[1] + pilv1_uncorr[2]*pilv1_uncorr[2] );
    normpi2 = 1./sqrt(pilv2_uncorr[0]*pilv2_uncorr[0] + pilv2_uncorr[1]*pilv2_uncorr[1] + pilv2_uncorr[2]*pilv2_uncorr[2] );
    normpi3 = 1./sqrt(pilv3_uncorr[0]*pilv3_uncorr[0] + pilv3_uncorr[1]*pilv3_uncorr[1] + pilv3_uncorr[2]*pilv3_uncorr[2] );
    pi1x = normpi1 * pilv1_uncorr[0] * sevt->track[pi1indx].p ;
    pi1y = normpi1 * pilv1_uncorr[1] * sevt->track[pi1indx].p ;
    pi1z = normpi1 *                   sevt->track[pi1indx].p ;
    pi2x = normpi2 * pilv2_uncorr[0] * sevt->track[pi2indx].p ;
    pi2y = normpi2 * pilv2_uncorr[1] * sevt->track[pi2indx].p ;
    pi2z = normpi2 *                   sevt->track[pi2indx].p ;
    pi3x = normpi3 * pilv3_uncorr[0] * sevt->track[pi3indx].p ;
    pi3y = normpi3 * pilv3_uncorr[1] * sevt->track[pi3indx].p ;
    pi3z = normpi3 *                   sevt->track[pi3indx].p ;
    pi1E = sqrt( pi1x*pi1x+pi1y*pi1y+pi1z*pi1z + massPionC*massPionC);
    pi2E = sqrt( pi2x*pi2x+pi2y*pi2y+pi2z*pi2z + massPionC*massPionC);
    pi3E = sqrt( pi3x*pi3x+pi3y*pi3y+pi3z*pi3z + massPionC*massPionC);


    pi1lv = new TLorentzVector(pi1x, pi1y, pi1z, pi1E);
    pi2lv = new TLorentzVector(pi2x, pi2y, pi2z, pi2E);
    pi3lv = new TLorentzVector(pi3x, pi3y, pi3z, pi3E);
    //eelv  = new TLorentzVector( (*el1lv)+(*el2lv) );
    //eelvmass = (*eelv).M();

    pi1trktime      = sevt->track[pi1indx].time ;
    pi1trkhodtime    = sevt->track[pi1indx].hodTime ;
    pi1cltime        = sevt->cluster[sevt->track[pi1indx].iClus].time;

    pi2trktime      = sevt->track[pi2indx].time ;
    pi2trkhodtime    = sevt->track[pi2indx].hodTime ;
    // icl_track       = sevt->track[eltrkindx2].iClus;
    pi2cltime        = sevt->cluster[sevt->track[pi2indx].iClus].time;

    pi3trktime       = sevt->track[pi3indx].time ;
    pi3trkhodtime    = sevt->track[pi3indx].hodTime ;
    pi3cltime        = sevt->cluster[sevt->track[pi3indx].iClus].time;

    if ( closap_double_(pilp1, pilp2, pilv1, pilv2, &cda_k3pic12, vertex_k3pic12) )
      k3pic12_vertex_OK  = TRUE;
    else
      k3pic12_vertex_OK  = FALSE;


    if ( closap_double_(pilp2, pilp3, pilv2, pilv3, &cda_k3pic23, vertex_k3pic23) )
      k3pic23_vertex_OK  = TRUE;
    else
      k3pic23_vertex_OK  = FALSE;

    if ( closap_double_(pilp1, pilp3, pilv1, pilv3, &cda_k3pic13, vertex_k3pic13) )
      k3pic13_vertex_OK  = TRUE;
    else
      k3pic13_vertex_OK  = FALSE;


    k3pic = new TLorentzVector( (*pi1lv) + (*pi2lv) + (*pi3lv) );

    pK1x = 0;
    pK1y = 0;                                   
    pK1z = 60.;                                     
    pKE =  sqrt( pK1z*pK1z + massKaonC*massKaonC);                             


    pK = new TLorentzVector(pK1x , pK1y , pK1z , pKE);


    k3pi_miss_mass = new TLorentzVector((*pK) - (*pi1lv) - (*pi2lv) - (*pi3lv) );
    // find vertex of mu and ee

    //eldistDCH1plane = sqrt( pow((elp1_uncorr[0]-elp2_uncorr[0]),2.) + pow((elp1_uncorr[1]-elp2_uncorr[1]),2.) ) ;

    //el1xLKrplane = sevt->track[eltrkindx1].x + (LKrz-DCHz) * sevt->track[eltrkindx1].dxdz ;
    //el1yLKrplane = sevt->track[eltrkindx1].y + (LKrz-DCHz) * sevt->track[eltrkindx1].dydz ;
    //el2xLKrplane = sevt->track[eltrkindx2].x + (LKrz-DCHz) * sevt->track[eltrkindx2].dxdz ;
    //el2yLKrplane = sevt->track[eltrkindx2].y + (LKrz-DCHz) * sevt->track[eltrkindx2].dydz ;
    //eldistLKrplane = sqrt( pow((el1xLKrplane-el2xLKrplane),2.) + pow((el1yLKrplane-el2yLKrplane),2.) ) ;
    //
    //muche1distDHC1plane = sqrt( pow((elp1_uncorr[0] - sevt->track[muchtrkindx].x),2.) + pow((elp1_uncorr[1] - sevt->track[muchtrkindx].y),2.) ) ;
    //muche2distDHC1plane = sqrt( pow((elp2_uncorr[0] - sevt->track[muchtrkindx].x),2.) + pow((elp2_uncorr[1] - sevt->track[muchtrkindx].y),2.) ) ;
    //muchxLKrplane         = sevt->track[muchtrkindx].x + (LKrz-DCHz) * sevt->track[muchtrkindx].dxdz ;
    //muchyLKrplane         = sevt->track[muchtrkindx].y + (LKrz-DCHz) * sevt->track[muchtrkindx].dydz ;
    //muche1distLKrplane     = sqrt( pow((el1xLKrplane-muchxLKrplane),2.) + pow((el1yLKrplane-muchxLKrplane),2.) ) ;
    //muche2distLKrplane     = sqrt( pow((el2xLKrplane-muchxLKrplane),2.) + pow((el2yLKrplane-muchxLKrplane),2.) ) ;

    eop_pi1_track=sevt->cluster[sevt->track[pi1indx].iClus].energy/sevt->track[pi1indx].p;
    eop_pi2_track=sevt->cluster[sevt->track[pi2indx].iClus].energy/sevt->track[pi2indx].p;
    eop_pi3_track=sevt->cluster[sevt->track[pi3indx].iClus].energy/sevt->track[pi3indx].p;
    
 
    //DCH INNER OUTER RADIUS CUTS FOR k3pic
    if (pow(pilp1_uncorr[0],2.) + pow(pilp1_uncorr[1],2.) < 12100. &&   
	pow(pilp2_uncorr[0],2.) + pow(pilp2_uncorr[1],2.) < 12100. &&
	pow(pilp3_uncorr[0],2.) + pow(pilp3_uncorr[1],2.) < 12100. &&
	pow(pilp1_uncorr[0],2.) + pow(pilp1_uncorr[1],2.) > 196.   &&
	pow(pilp2_uncorr[0],2.) + pow(pilp2_uncorr[1],2.) > 196.   &&
	pow(pilp3_uncorr[0],2.) + pow(pilp3_uncorr[1],2.) > 196.   &&
	(vertex_k3pic12[2] > -1500.  || vertex_k3pic12[2] < 7000.) &&
	(vertex_k3pic13[2] > -1500.  || vertex_k3pic13[2] < 7000.) &&
	(vertex_k3pic23[2] > -1500.  || vertex_k3pic23[2] < 7000.) &&
	fabs(vertex_k3pic12[2] - vertex_k3pic23[2] ) < 100         &&
	fabs(vertex_k3pic12[2] - vertex_k3pic13[2] ) < 100         &&
	fabs(vertex_k3pic23[2] - vertex_k3pic13[2] ) < 100         &&
	k3pic23_vertex_OK==1                                       &&
	k3pic13_vertex_OK==1                                       &&
	k3pic12_vertex_OK==1                                       &&
	(*k3pic).M() > 0.4837                                      &&
	(*k3pic).M() < 0.5037         
	){
      
      ((TH1F*)gDirectory->FindObject(h_eop_pi3_track->GetName()))->Fill(eop_pi3_track);
      ((TH1F*)gDirectory->FindObject(h_eop_pi1_track->GetName()))->Fill(eop_pi1_track);
      ((TH1F*)gDirectory->FindObject(h_eop_pi2_track->GetName()))->Fill(eop_pi2_track);
      ((TH1F*)gDirectory->FindObject(h_eop_pi_track->GetName()))->Fill(eop_pi3_track);
      ((TH1F*)gDirectory->FindObject(h_eop_pi_track->GetName()))->Fill(eop_pi1_track);
      ((TH1F*)gDirectory->FindObject(h_eop_pi_track->GetName()))->Fill(eop_pi2_track);


      //((TH1F*)gDirectory->FindObject(h_k3pi_nu_m->GetName()))->Fill( (*k3pi_miss_mass).M2() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_3track_E->GetName()))->Fill( (*k3pic).E());
      ((TH1F*)gDirectory->FindObject(h_k3pi_3track_P->GetName()))->Fill( (*k3pic).P());
      ((TH1F*)gDirectory->FindObject(h_k3pi_3track_Pt->GetName()))->Fill((*k3pic).Pt());
      ((TH1F*)gDirectory->FindObject(h_k3pi_3track_m->GetName()))->Fill ((*k3pic).M());
      ((TH1F*)gDirectory->FindObject(h_k3pi_nu_E->GetName()))->Fill( (*k3pi_miss_mass).E());
      ((TH1F*)gDirectory->FindObject(h_k3pi_nu_P->GetName()))->Fill( (*k3pi_miss_mass).P());
      ((TH1F*)gDirectory->FindObject(h_k3pi_nu_Pt->GetName()))->Fill((*k3pi_miss_mass).Pt());
      ((TH1F*)gDirectory->FindObject(h_k3pi_nu_m->GetName()))->Fill ((*k3pi_miss_mass).M());
 
      //((TH1F*)gDirectory->FindObject(h_3trk_invm->GetName()))->Fill ((*k3pic).M2());
 
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_P->GetName()))->Fill( (*pi1lv).P() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_P->GetName()))->Fill( (*pi2lv).P() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_P->GetName()))->Fill( (*pi3lv).P() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_E->GetName()))->Fill( (*pi1lv).E() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_E->GetName()))->Fill( (*pi2lv).E() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_E->GetName()))->Fill( (*pi3lv).E() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_Pt->GetName()))->Fill((*pi1lv).Pt() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_Pt->GetName()))->Fill((*pi2lv).Pt() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_Pt->GetName()))->Fill((*pi3lv).Pt() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_m->GetName()))->Fill( (*pi1lv).M() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_m->GetName()))->Fill( (*pi2lv).M() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_m->GetName()))->Fill( (*pi3lv).M() );
    
      if(sevt->track[pi1indx].q * (sevt->track[pi1indx].q + sevt->track[pi2indx].q + sevt->track[pi3indx].q) == -1){                                                                           
	((TH1F*)gDirectory->FindObject(h_piodd_eovp->GetName()))->Fill(eop_pi1_track);                  
      } else if (sevt->track[pi2indx].q * (sevt->track[pi1indx].q + sevt->track[pi2indx].q + sevt->track[pi3indx].q) == -1){                                                                    
	((TH1F*)gDirectory->FindObject(h_piodd_eovp->GetName()))->Fill( eop_pi2_track); 
      }  else if ( sevt->track[pi3indx].q * (sevt->track[pi1indx].q + sevt->track[pi2indx].q + sevt->track[pi3indx].q)== -1){                                                                  
	((TH1F*)gDirectory->FindObject(h_piodd_eovp->GetName()))->Fill( eop_pi3_track);                  
      }   
    
      ((TH1F*)gDirectory->FindObject(h_pi1trktime->GetName()))           ->Fill( pi1trktime );
      ((TH1F*)gDirectory->FindObject(h_pi1trkhodtime->GetName()))        ->Fill( pi1trkhodtime );
      ((TH1F*)gDirectory->FindObject(h_pi1trktimediff->GetName()))       ->Fill( pi1trktime - pi1trkhodtime );
      ((TH1F*)gDirectory->FindObject(h_ecltimepi1->GetName()))           ->Fill( pi1cltime );
      ((TH1F*)gDirectory->FindObject(h_pi12hodtimediff->GetName()))      ->Fill( pi1trkhodtime - pi2trkhodtime );
      ((TH1F*)gDirectory->FindObject(h_pi13hodtimediff->GetName()))      ->Fill( pi1trkhodtime - pi3trkhodtime );
      ((TH1F*)gDirectory->FindObject(h_pi23hodtimediff->GetName()))      ->Fill( pi2trkhodtime - pi3trkhodtime );
      
      ((TH1F*)gDirectory->FindObject(h_cda_k3pic12->GetName()))->Fill(cda_k3pic12);
      ((TH1F*)gDirectory->FindObject(h_k3pi_vx_12->GetName()))->Fill(vertex_k3pic12[0]);
      ((TH1F*)gDirectory->FindObject(h_k3pi_vy_12->GetName()))->Fill(vertex_k3pic12[1]);
      ((TH1F*)gDirectory->FindObject(h_k3pi_vz_12->GetName()))->Fill(vertex_k3pic12[2]);
     
      ((TH1F*)gDirectory->FindObject(h_cda_k3pic13->GetName()))->Fill(cda_k3pic13);
      ((TH1F*)gDirectory->FindObject(h_k3pi_vx_13->GetName()))->Fill(vertex_k3pic13[0]);
      ((TH1F*)gDirectory->FindObject(h_k3pi_vy_13->GetName()))->Fill(vertex_k3pic13[1]);
      ((TH1F*)gDirectory->FindObject(h_k3pi_vz_13->GetName()))->Fill(vertex_k3pic13[2]);
 
      ((TH1F*)gDirectory->FindObject(h_cda_k3pic23->GetName()))->Fill(cda_k3pic23);
      ((TH1F*)gDirectory->FindObject(h_k3pi_vx_23->GetName()))->Fill(vertex_k3pic23[0]);
      ((TH1F*)gDirectory->FindObject(h_k3pi_vy_23->GetName()))->Fill(vertex_k3pic23[1]);
      ((TH1F*)gDirectory->FindObject(h_k3pi_vz_23->GetName()))->Fill(vertex_k3pic23[2]);

      ((TH1F*)gDirectory->FindObject(h_vz_pi12_13_diff->GetName()))->Fill( vertex_k3pic12[2] -vertex_k3pic13[2] );
      ((TH1F*)gDirectory->FindObject(h_vz_pi13_23_diff->GetName()))->Fill( vertex_k3pic13[2] -vertex_k3pic23[2] );  
      ((TH1F*)gDirectory->FindObject(h_vz_pi12_23_diff->GetName()))->Fill( vertex_k3pic12[2] -vertex_k3pic23[2] );
 
    }//DCH inner outer radius and all other cuts
  }
  //  ------------------------------ END OF K3PIC selection part ------------------------------------

  if (eltrkindx1 == -1 || eltrkindx2 == -1 || muchtrkindx==-1)    {cleanup();return 0;}

 


  if (nelectrack!=2)        {cleanup();return 0;}

  
  
  //Making four vectors and time variables of the tracks
  //0000
  //Changes
  
  elp1_uncorr[0] = sevt->track[eltrkindx1].bx;
  elp1_uncorr[1] = sevt->track[eltrkindx1].by;
  elp1_uncorr[2] = DCHbz;
  elp2_uncorr[0] = sevt->track[eltrkindx2].bx;
  elp2_uncorr[1] = sevt->track[eltrkindx2].by;
  elp2_uncorr[2] = DCHbz;

  elv1_uncorr[0] = sevt->track[eltrkindx1].bdxdz;
  elv1_uncorr[1] = sevt->track[eltrkindx1].bdydz;
  elv1_uncorr[2] = 1.;
  elv2_uncorr[0] = sevt->track[eltrkindx2].bdxdz;
  elv2_uncorr[1] = sevt->track[eltrkindx2].bdydz;
  elv2_uncorr[2] = 1.;

  // if (IS_DATA){
  //cout << "Index breaks when eltrkindx1 =" << eltrkindx1 << " \t and eltrkindx2 = " << eltrkindx2 << endl;
  //blue field corrected values
  elv1[0] = trkCorrSlopesV[eltrkindx1]->x();
  elv1[1] = trkCorrSlopesV[eltrkindx1]->y();
  elv1[2] = 1.;
  elp1[0] = trkCorrMidPointsV[eltrkindx1]->x();
  elp1[1] = trkCorrMidPointsV[eltrkindx1]->y();
  elp1[2] = trkCorrMidPointsV[eltrkindx1]->z();

  //cout << "Slopes e1 =" << trkCorrSlopesV[eltrkindx1]->x() << "=1=  "<< trkCorrSlopesV[eltrkindx1]->y() << "=2=  " << trkCorrSlopesV[eltrkindx1]->z() << "=3=  " << eltrkindx1  << endl;
  //cout << "Slopes e2 =" << trkCorrSlopesV[eltrkindx2]->x() << "=1=  "<< trkCorrSlopesV[eltrkindx2]->y() << "=2=  " << trkCorrSlopesV[eltrkindx2]->z() << "=3=  " << eltrkindx2  << endl;
  
  
  //cout << "Mid points e1 =" <<  trkCorrMidPointsV[eltrkindx1]->x() << "=1=  "<< trkCorrMidPointsV[eltrkindx1]->x() << "=2=  " << trkCorrMidPointsV[eltrkindx1]->x() << "=3=  " << eltrkindx1 << endl;
  //cout << "Mid points e2 =" <<  trkCorrMidPointsV[eltrkindx2]->x() << "=1=  "<< trkCorrMidPointsV[eltrkindx2]->x() << "=2=  " << trkCorrMidPointsV[eltrkindx2]->x() << "=3=  " << eltrkindx2 << endl;

  
  
 
  
  elv2[0] = trkCorrSlopesV[eltrkindx2]->x();
  elv2[1] = trkCorrSlopesV[eltrkindx2]->y();
  elv2[2] = 1.;
  elp2[0] = trkCorrMidPointsV[eltrkindx2]->x();
  elp2[1] = trkCorrMidPointsV[eltrkindx2]->y();
  elp2[2] = trkCorrMidPointsV[eltrkindx2]->z();

  //}
  norm1 = 1./sqrt(elv1_uncorr[0]*elv1_uncorr[0] + elv1_uncorr[1]*elv1_uncorr[1] + elv1_uncorr[2]*elv1_uncorr[2] );
  norm2 = 1./sqrt(elv2_uncorr[0]*elv2_uncorr[0] + elv2_uncorr[1]*elv2_uncorr[1] + elv2_uncorr[2]*elv2_uncorr[2] );
  p1x = norm1 * elv1_uncorr[0] * sevt->track[eltrkindx1].p ;
  p1y = norm1 * elv1_uncorr[1] * sevt->track[eltrkindx1].p ;
  p1z = norm1 *                  sevt->track[eltrkindx1].p ;
  p2x = norm2 * elv2_uncorr[0] * sevt->track[eltrkindx2].p ;
  p2y = norm2 * elv2_uncorr[1] * sevt->track[eltrkindx2].p ;
  p2z = norm2 *                  sevt->track[eltrkindx2].p ;
  el1E = sqrt( p1x*p1x+p1y*p1y+p1z*p1z + massElec*massElec);
  el2E = sqrt( p2x*p2x+p2y*p2y+p2z*p2z + massElec*massElec);



  el1lv = new TLorentzVector(p1x, p1y, p1z, el1E);
  el2lv = new TLorentzVector(p2x, p2y, p2z, el2E);
  eelv  = new TLorentzVector( (*el1lv)+(*el2lv) );
  eelvmass = (*eelv).M();

  el1trktime      = sevt->track[eltrkindx1].time ;
  el1trkhodtime    = sevt->track[eltrkindx1].hodTime ;
  el1cltime        = sevt->cluster[sevt->track[eltrkindx1].iClus].time;

  el2trktime      = sevt->track[eltrkindx2].time ;
  el2trkhodtime    = sevt->track[eltrkindx2].hodTime ;
  // icl_track       = sevt->track[eltrkindx2].iClus;
  el2cltime        = sevt->cluster[sevt->track[eltrkindx2].iClus].time;

  if ( closap_double_(elp1, elp2, elv1, elv2, &cda, vertex) )
    vertex_OK = TRUE;
  else
    vertex_OK = FALSE;

 

  muchv_uncorr[0] = sevt->track[muchtrkindx].bdxdz;
  muchv_uncorr[1] = sevt->track[muchtrkindx].bdydz;
  muchv_uncorr[2] = 1.;
  muchp_uncorr[0] = sevt->track[muchtrkindx].bx;
  muchp_uncorr[1] = sevt->track[muchtrkindx].by;
  muchp_uncorr[2] = DCHbz;
  
  //cout << "mu =" << muchv_uncorr[0] << "=chv[0]=" << muchv_uncorr[1] << "=chv[1]=" << muchv_uncorr[2] << "=chv[2]=" << muchp_uncorr[0] << "=chp[0]=" << muchp_uncorr[1] << "=chp[1]=" << muchp_uncorr[2] << "=chp[2]=" << endl;
  //cout << "mu =" << trkCorrSlopesV[eltrkindx2]->x() << endl;
  //ALERT
  // if(IS_DATA){
  muchv[0] = trkCorrSlopesV[muchtrkindx]->x();
  muchv[1] = trkCorrSlopesV[muchtrkindx]->y();
  muchv[2] = 1.;

  muchp[0] = trkCorrMidPointsV[muchtrkindx]->x();
  muchp[1] = trkCorrMidPointsV[muchtrkindx]->y();
  muchp[2] = trkCorrMidPointsV[muchtrkindx]->z();
  
  // }
  
 
  normmu = 1./sqrt(muchv_uncorr[0]*muchv_uncorr[0] + muchv_uncorr[1]*muchv_uncorr[1] + muchv_uncorr[2]*muchv_uncorr[2] );
  pmux = normmu * muchv_uncorr[0] * sevt->track[muchtrkindx].p ;
  pmuy = normmu * muchv_uncorr[1] *  sevt->track[muchtrkindx].p ;
  pmuz = normmu *                   sevt->track[muchtrkindx].p ;
  //Pions
  //muchE = sqrt( pmux*pmux+pmuy*pmuy+pmuz*pmuz + massPionC*massPionC);
  muchE = sqrt( pmux*pmux+pmuy*pmuy+pmuz*pmuz + massMuonC*massMuonC);
  muchlv = new TLorentzVector(pmux, pmuy, pmuz, muchE);

  muchtrktime         = sevt->track[muchtrkindx].time;
  muchcltime         = sevt->cluster[ sevt->track[muchtrkindx].iClus ].time;
  muchtrkhodtime     = sevt->track[muchtrkindx].hodTime;



  if ( closap_double_(elp1, muchp, elv1, muchv, &cda_e1much, vertex_e1much) )
    muchee_vertex_OK = TRUE;
  else
    muchee_vertex_OK = FALSE;
  //2nd vertex
  if ( closap_double_(elp2, muchp, elv2, muchv, &cda_e2much, vertex_e2much) )
    muchee2_vertex_OK = TRUE;
  else
    muchee2_vertex_OK = FALSE;
  
  muvee = new TLorentzVector( (*el1lv) + (*el2lv) + (*muchlv) );

  pK1x = 0;
  pK1y = 0;                                   
  pK1z = 60.;                                     
  pKE =  sqrt( pK1z*pK1z + massKaonC*massKaonC);                             


  pK = new TLorentzVector(pK1x , pK1y , pK1z , pKE);


  neutrino = new TLorentzVector((*pK) - (*muchlv) - (*el1lv) - (*el2lv));
  // find vertex of mu and ee
  
  eldistDCH1plane = sqrt( pow((elp1_uncorr[0]-elp2_uncorr[0]),2.) + pow((elp1_uncorr[1]-elp2_uncorr[1]),2.) ) ;

  el1xLKrplane = sevt->track[eltrkindx1].x + (LKrz-DCHz) * sevt->track[eltrkindx1].dxdz ;
  el1yLKrplane = sevt->track[eltrkindx1].y + (LKrz-DCHz) * sevt->track[eltrkindx1].dydz ;
  el2xLKrplane = sevt->track[eltrkindx2].x + (LKrz-DCHz) * sevt->track[eltrkindx2].dxdz ;
  el2yLKrplane = sevt->track[eltrkindx2].y + (LKrz-DCHz) * sevt->track[eltrkindx2].dydz ;
  eldistLKrplane = sqrt( pow((el1xLKrplane-el2xLKrplane),2.) + pow((el1yLKrplane-el2yLKrplane),2.) ) ;

  muche1distDHC1plane = sqrt( pow((elp1_uncorr[0] - sevt->track[muchtrkindx].x),2.) + pow((elp1_uncorr[1] - sevt->track[muchtrkindx].y),2.) ) ;
  muche2distDHC1plane = sqrt( pow((elp2_uncorr[0] - sevt->track[muchtrkindx].x),2.) + pow((elp2_uncorr[1] - sevt->track[muchtrkindx].y),2.) ) ;
  muchxLKrplane         = sevt->track[muchtrkindx].x + (LKrz-DCHz) * sevt->track[muchtrkindx].dxdz ;
  muchyLKrplane         = sevt->track[muchtrkindx].y + (LKrz-DCHz) * sevt->track[muchtrkindx].dydz ;
  muche1distLKrplane     = sqrt( pow((el1xLKrplane-muchxLKrplane),2.) + pow((el1yLKrplane-muchxLKrplane),2.) ) ;
  muche2distLKrplane     = sqrt( pow((el2xLKrplane-muchxLKrplane),2.) + pow((el2yLKrplane-muchxLKrplane),2.) ) ;
  
  eop_mu_track=sevt->cluster[sevt->track[muchtrkindx].iClus].energy/sevt->track[muchtrkindx].p;
  eop_e1_track=sevt->cluster[sevt->track[eltrkindx1].iClus].energy/sevt->track[eltrkindx1].p;
  eop_e2_track=sevt->cluster[sevt->track[eltrkindx2].iClus].energy/sevt->track[eltrkindx2].p;
  

  //end 0000
  //Histos kinematics + timing of the three tracks
 
  
  //v_second_e1 = new TLorentzVector;
  
  ((TH1F*)gDirectory->FindObject(h_el1trktime->GetName()))        ->Fill( el1trktime );
  ((TH1F*)gDirectory->FindObject(h_el1trkhodtime->GetName()))        ->Fill( el1trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_el1trktimediff->GetName()))    ->Fill( el1trktime - el1trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_ecltime1->GetName()))            ->Fill( el1cltime );
  ((TH1F*)gDirectory->FindObject(h_e1trkcltimediff->GetName()))    ->Fill( el1trktime- el1cltime );

  ((TH1F*)gDirectory->FindObject(h_el2trktime->GetName()))        ->Fill( el2trktime );
  ((TH1F*)gDirectory->FindObject(h_el2trkhodtime->GetName()))        ->Fill( el2trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_el2trktimediff->GetName()))    ->Fill( el2trktime - el2trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_ecltime2->GetName()))            ->Fill( el2cltime );
  ((TH1F*)gDirectory->FindObject(h_e2trkcltimediff->GetName()))    ->Fill( el2trktime-el2cltime );

  ((TH1F*)gDirectory->FindObject(h_e1e2timediff->GetName()))        ->Fill( el1trktime-el2trktime );
  ((TH1F*)gDirectory->FindObject(h_e1e2cltimediff->GetName()))    ->Fill( el1cltime-el2cltime );
  ((TH1F*)gDirectory->FindObject(h_e1e2hodtimediff->GetName()))        ->Fill( el1trkhodtime-el2trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_e1muhodtimediff->GetName()))        ->Fill( el1trkhodtime-muchtrkhodtime );
  ((TH1F*)gDirectory->FindObject(h_e2muhodtimediff->GetName()))        ->Fill( el2trkhodtime-muchtrkhodtime );


  ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_mu_track);
  ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_e1_track);
  ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_e2_track);

  
  ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_mu_track, (*muchlv).P());
  ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_e1_track, (*el1lv).P() );
  ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_e2_track, (*el2lv).P() );
 
  ((TH1F*)gDirectory->FindObject(h_eop_mu->GetName()))->Fill(eop_mu_track);
  ((TH2F*)gDirectory->FindObject(h_eop_mu_vs_p->GetName()))->Fill(eop_mu_track, (*muchlv).P());
  ((TH1F*)gDirectory->FindObject(h_eop_el->GetName()))->Fill(eop_e1_track);
  ((TH2F*)gDirectory->FindObject(h_eop_el_vs_p->GetName()))->Fill(eop_e1_track, (*el1lv).P());
  ((TH1F*)gDirectory->FindObject(h_eop_el->GetName()))->Fill(eop_e2_track);
  ((TH2F*)gDirectory->FindObject(h_eop_el_vs_p->GetName()))->Fill(eop_e2_track, (*el2lv).P());


  ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill( elp1_uncorr[0] );
  ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill( elp1_uncorr[1] );
  ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill( elv1_uncorr[0] );
  ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill( elv1_uncorr[1] );
  ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill( elp1_uncorr[0], elp1_uncorr[1] );
  ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill( elv1_uncorr[0], elv1_uncorr[1] );

  ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( el1lv->P() );
  ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( el2lv->P() );
  ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el1lv->E() );
  ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el2lv->E() );
  ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el1lv->P() );
  ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el2lv->P() );
  ((TH1F*)gDirectory->FindObject(h_e1m->GetName()))->Fill( el1lv->M() );
  ((TH1F*)gDirectory->FindObject(h_e2m->GetName()))->Fill( el2lv->M() );
  ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el1lv->Pt() );
  ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el2lv->Pt() );
  ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
  ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
  ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );
  ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);
  ((TH1F*)gDirectory->FindObject(h_angle_e1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()) );
  ((TH1F*)gDirectory->FindObject(h_angle_e1e2_mu->GetName()))->Fill( (*eelv).Vect().Angle(muchlv->Vect()) );

  ((TH1F*)gDirectory->FindObject(h_muchtrktime->GetName()))        ->Fill(muchtrktime);
  ((TH1F*)gDirectory->FindObject(h_muchtrkhodtime->GetName()))        ->Fill( muchtrkhodtime );
  ((TH1F*)gDirectory->FindObject(h_muchtrktimediff->GetName()))    ->Fill( muchtrktime-muchtrkhodtime );
  //((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( muchtrktime-muchcltime );

  ((TH1F*)gDirectory->FindObject(h_muchlvE->GetName()))->Fill( (*muchlv).E() );
  ((TH1F*)gDirectory->FindObject(h_muchlvP->GetName()))->Fill( (*muchlv).P() );
  ((TH1F*)gDirectory->FindObject(h_muchlvPt->GetName()))->Fill((*muchlv).Pt() );

  ((TH1F*)gDirectory->FindObject(h_KE->GetName()))->Fill((*pK).E());
  ((TH1F*)gDirectory->FindObject(h_KP->GetName()))->Fill((*pK).P());

  ((TH1F*)gDirectory->FindObject(h_nu_P->GetName()))->Fill(  (*neutrino).P() );
  ((TH1F*)gDirectory->FindObject(h_nu_E->GetName()))->Fill(  (*neutrino).E() );
  ((TH1F*)gDirectory->FindObject(h_nu_Pt->GetName()))->Fill( (*neutrino).Pt() );
  ((TH1F*)gDirectory->FindObject(h_nu_m->GetName()))->Fill(  (*neutrino).M() );
  ((TH1F*)gDirectory->FindObject(h_missing_mass->GetName()))->Fill( (*neutrino).M2() );


  ((TH1F*)gDirectory->FindObject(h_mueeE->GetName()))->Fill( (*muvee).E());
  ((TH1F*)gDirectory->FindObject(h_mueeP->GetName()))->Fill( (*muvee).P());
  ((TH1F*)gDirectory->FindObject(h_mueePt->GetName()))->Fill((*muvee).Pt());
  ((TH1F*)gDirectory->FindObject(h_mueeM->GetName()))->Fill ((*muvee).M());
  ((TH1F*)gDirectory->FindObject(h_3trk_invm->GetName()))->Fill ((*muvee).M2());


  ((TH1F*)gDirectory->FindObject(h_bx_much->GetName()))            ->Fill(muchp_uncorr[0]);
  ((TH1F*)gDirectory->FindObject(h_by_much->GetName()))            ->Fill(muchp_uncorr[1]);
  ((TH1F*)gDirectory->FindObject(h_bdxdz_much->GetName()))            ->Fill(muchv_uncorr[0]);
  ((TH1F*)gDirectory->FindObject(h_bdydz_much->GetName()))            ->Fill(muchv_uncorr[1]) ;
  ((TH1F*)gDirectory->FindObject(h_bx_vs_by_much->GetName()))        ->Fill(muchp_uncorr[0], muchp_uncorr[1]);
  ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_much->GetName()))->Fill(muchv_uncorr[0], muchv_uncorr[1]);

  ((TH1F*)gDirectory->FindObject(h_cda_ee->GetName()))->Fill(cda);
  ((TH1F*)gDirectory->FindObject(h_vx->GetName()))->Fill(vertex[0]);
  ((TH1F*)gDirectory->FindObject(h_vy->GetName()))->Fill(vertex[1]);
  ((TH1F*)gDirectory->FindObject(h_vz->GetName()))->Fill(vertex[2]);
  ((TH2F*)gDirectory->FindObject(h_vx_vs_vy->GetName()))->Fill(vertex[0],vertex[1]);

  ((TH1F*)gDirectory->FindObject(h_vx_e1much->GetName()))->Fill(vertex_e1much[0]);
  ((TH1F*)gDirectory->FindObject(h_vy_e1much->GetName()))->Fill(vertex_e1much[1]);
  ((TH1F*)gDirectory->FindObject(h_vz_e1much->GetName()))->Fill(vertex_e1much[2]);
  ((TH1F*)gDirectory->FindObject(h_cda_e1much->GetName()))->Fill(cda_e1much);
  ((TH1F*)gDirectory->FindObject(h_vz_e1much_diff->GetName()))->Fill( vertex[2] - vertex_e1much[2] );
  
  ((TH1F*)gDirectory->FindObject(h_vx_e2much->GetName()))->Fill(vertex_e2much[0]);
  ((TH1F*)gDirectory->FindObject(h_vy_e2much->GetName()))->Fill(vertex_e2much[1]);
  ((TH1F*)gDirectory->FindObject(h_vz_e2much->GetName()))->Fill(vertex_e2much[2]);
  ((TH1F*)gDirectory->FindObject(h_cda_e2much->GetName()))->Fill(cda_e2much);
  ((TH1F*)gDirectory->FindObject(h_vz_e2much_diff->GetName()))->Fill( vertex[2] - vertex_e2much[2] );
  ((TH1F*)gDirectory->FindObject(h_vz_mu1mu2_diff->GetName()))->Fill( vertex_e1much[2] - vertex_e2much[2] );

  ((TH1F*)gDirectory->FindObject(h_cda_e1much->GetName()))->Fill(cda_e1much);
  ((TH2F*)gDirectory->FindObject(h_vx_vs_vy_e1much->GetName()))->Fill(vertex_e1much[0],vertex_e1much[1]);

  ((TH1F*)gDirectory->FindObject(h_muche1distDHC1plane->GetName()))->Fill( muche1distDHC1plane );
  ((TH1F*)gDirectory->FindObject(h_muche2distDHC1plane->GetName()))->Fill( muche2distDHC1plane );
  ((TH1F*)gDirectory->FindObject(h_muche1distLKrplane->GetName()))->Fill(  muche1distLKrplane );
  ((TH1F*)gDirectory->FindObject(h_muche2distLKrplane->GetName()))->Fill(  muche2distLKrplane );
  
  ((TH1F*)gDirectory->FindObject(h_eldistDCH1plane->GetName()))->Fill( eldistDCH1plane );
  ((TH1F*)gDirectory->FindObject(h_eldistLKrplane->GetName()))->Fill( eldistLKrplane );
  file1->cd();
  h_cutflow->Fill(cutcounter);
  h_cutflow->GetXaxis()->SetBinLabel( cutcounter, cuts[cutcounter-1].c_str() );// "Initial cuts" // "cluster_isolation");
  cutcounter++;
  ///THE END OF ************* Initial cuts ****************** ///


  /// ************* Geometry cuts ****************** ///

  //"e1e2_DCHouterR"
  if( pow(elp1_uncorr[0],2.) + pow(elp1_uncorr[1],2.) > 12100. )    {cleanup();return 0;} //110*110
  if( pow(elp2_uncorr[0],2.) + pow(elp2_uncorr[1],2.) > 12100. )    {cleanup();return 0;}
  //"e1e2_DCHinnerR"
  if( pow(elp1_uncorr[0],2.) + pow(elp1_uncorr[1],2.) < 196. )    {cleanup();return 0;} //14*14 /it was 12*12/
  if( pow(elp2_uncorr[0],2.) + pow(elp2_uncorr[1],2.) < 196. )    {cleanup();return 0;}
  //"e1e2_dist_LKr"
  //if (flux_id ==1 ){
  //if(eldistDCH1plane<2.)       {cleanup();return 0;}//return 0;
  //if( eldistLKrplane<15. )        {cleanup();return 0;}
  // }
  // if (flux_id==0){
  if(eldistDCH1plane<1.)       {cleanup();return 0;}//return 0;
  if( eldistLKrplane<20. )        {cleanup();return 0;}
  // }
  //"e1e2_dist_DCH1"

  //if(muche1distLKrplane<30.)    {cleanup();return 0;}
  //if(muche2distLKrplane<30.)    {cleanup();return 0;}

  dirs[cutcounter-1]->cd();

  if (IS_MC){
    //if(true_4.P()  < 3. ||true_5.P() < 3.) {cleanup();return 0;}
    ((TH1D*)gDirectory->FindObject( h_mc_muvee_M->GetName()))->Fill( Muee_three_track.M());
    ((TH1D*)gDirectory->FindObject( h_mc_muvee_P->GetName()))->Fill( Muee_three_track.P());
    ((TH1D*)gDirectory->FindObject( h_mc_muvee_Pt->GetName()))->Fill( Muee_three_track.Pt());
    ((TH1D*)gDirectory->FindObject( h_mc_elsP->GetName()))->Fill( mc_ee.P());
    ((TH1D*)gDirectory->FindObject( h_mc_elsM->GetName()))->Fill( mc_ee.M());
    ((TH1D*)gDirectory->FindObject( h_mc_nu_M->GetName()))->Fill( Nu_true.M());
    ((TH1D*)gDirectory->FindObject( h_mc_nu_M2->GetName()))->Fill( Nu_true.M2());
    ((TH1D*)gDirectory->FindObject( h_mc_nu_P->GetName()))->Fill( Nu_true.P());
    ((TH1D*)gDirectory->FindObject( h_mc_nu_Pt->GetName()))->Fill(Nu_true.Pt());
  

    ((TH1D*)gDirectory->FindObject(h_true_M_0->GetName()))->Fill( true_0.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_0->GetName()))->Fill( true_0.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_0->GetName()))->Fill(true_0.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_0->GetName()))->Fill( true_0.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_1->GetName()))->Fill( true_1.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_1->GetName()))->Fill( true_1.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_1->GetName()))->Fill(true_1.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_1->GetName()))->Fill( true_1.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_2->GetName()))->Fill( true_2.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_2->GetName()))->Fill( true_2.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_2->GetName()))->Fill(true_2.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_2->GetName()))->Fill( true_2.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_3->GetName()))->Fill( true_3.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_3->GetName()))->Fill( true_3.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_3->GetName()))->Fill(true_3.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_3->GetName()))->Fill( true_3.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_4.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_4.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_4.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_4.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_5->GetName()))->Fill( true_5.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_5->GetName()))->Fill( true_5.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_5->GetName()))->Fill(true_5.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_5->GetName()))->Fill( true_5.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_6->GetName()))->Fill( true_6.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_6->GetName()))->Fill( true_6.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_6->GetName()))->Fill(true_6.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_6->GetName()))->Fill( true_6.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_7->GetName()))->Fill( true_7.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_7->GetName()))->Fill( true_7.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_7->GetName()))->Fill(true_7.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_7->GetName()))->Fill( true_7.E());
  }
  
  
  
  ((TH1D*)gDirectory->FindObject( h_vz_initial->GetName()))->Fill( sevt->vtx[0].z);
  ((TH1I*)gDirectory->FindObject(h_ntrack->GetName()))->Fill(ntrack);
  ((TH1I*)gDirectory->FindObject(h_nclust->GetName()))->Fill(nclust);
  ((TH1I*)gDirectory->FindObject(h_evttypebytrks->GetName()))->Fill(evttypebytrks);
  ((TH1I*)gDirectory->FindObject(h_nelectrack->GetName()))->Fill(nelectrack);
  ((TH1I*)gDirectory->FindObject(h_mutrack->GetName()))->Fill(nmutrack);
  ((TH2I*)gDirectory->FindObject(h_neltrk_mutrk->GetName()))->Fill(nelectrack, nmutrack);
  ((TH1I*)gDirectory->FindObject(h_mutrackq->GetName()))->Fill(muchtrkq);
  ((TH1I*)gDirectory->FindObject(h_e1trackq->GetName()))->Fill(eltrkq1);
  ((TH1I*)gDirectory->FindObject(h_e2trackq->GetName()))->Fill(eltrkq2);
  
  ((TH1F*)gDirectory->FindObject(h_el1trktime->GetName()))        ->Fill( el1trktime );
  ((TH1F*)gDirectory->FindObject(h_el1trkhodtime->GetName()))        ->Fill( el1trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_el1trktimediff->GetName()))    ->Fill( el1trktime - el1trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_ecltime1->GetName()))            ->Fill( el1cltime );
  ((TH1F*)gDirectory->FindObject(h_e1trkcltimediff->GetName()))    ->Fill( el1trktime- el1cltime );
  
  ((TH1F*)gDirectory->FindObject(h_el2trktime->GetName()))        ->Fill( el2trktime );
  ((TH1F*)gDirectory->FindObject(h_el2trkhodtime->GetName()))        ->Fill( el2trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_el2trktimediff->GetName()))    ->Fill( el2trktime - el2trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_ecltime2->GetName()))            ->Fill( el2cltime );
  ((TH1F*)gDirectory->FindObject(h_e2trkcltimediff->GetName()))    ->Fill( el2trktime-el2cltime );
  
  ((TH1F*)gDirectory->FindObject(h_e1e2timediff->GetName()))        ->Fill( el1trktime-el2trktime );
  ((TH1F*)gDirectory->FindObject(h_e1e2cltimediff->GetName()))    ->Fill( el1cltime-el2cltime );
  ((TH1F*)gDirectory->FindObject(h_e1e2hodtimediff->GetName()))        ->Fill( el1trkhodtime-el2trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_e1muhodtimediff->GetName()))        ->Fill( el1trkhodtime-muchtrkhodtime );
  ((TH1F*)gDirectory->FindObject(h_e2muhodtimediff->GetName()))        ->Fill( el2trkhodtime-muchtrkhodtime );
  

  ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_mu_track);
  ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_e1_track);
  ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_e2_track);

  
  ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_mu_track, (*muchlv).P());
  ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_e1_track, (*el1lv).P() );
  ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_e2_track, (*el2lv).P() );
 
  ((TH1F*)gDirectory->FindObject(h_eop_mu->GetName()))->Fill(eop_mu_track);
  ((TH2F*)gDirectory->FindObject(h_eop_mu_vs_p->GetName()))->Fill(eop_mu_track, (*muchlv).P());
  ((TH1F*)gDirectory->FindObject(h_eop_el->GetName()))->Fill(eop_e1_track);
  ((TH2F*)gDirectory->FindObject(h_eop_el_vs_p->GetName()))->Fill(eop_e1_track, (*el1lv).P());
  ((TH1F*)gDirectory->FindObject(h_eop_el->GetName()))->Fill(eop_e2_track);
  ((TH2F*)gDirectory->FindObject(h_eop_el_vs_p->GetName()))->Fill(eop_e2_track, (*el2lv).P());


  
  
   
 
  ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill( elp1_uncorr[0] );
  ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill( elp1_uncorr[1] );
  ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill( elv1_uncorr[0] );
  ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill( elv1_uncorr[1] );
  ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill( elp1_uncorr[0], elp1_uncorr[1] );
  ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill( elv1_uncorr[0], elv1_uncorr[1] );
  
  ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( el1lv->P() );
  ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( el2lv->P() );
  ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el1lv->E() );
  ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el2lv->E() );
  ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el1lv->P() );
  ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el2lv->P() );
  ((TH1F*)gDirectory->FindObject(h_e1m->GetName()))->Fill( el1lv->M() );
  ((TH1F*)gDirectory->FindObject(h_e2m->GetName()))->Fill( el2lv->M() );
  ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el1lv->Pt() );
  ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el2lv->Pt() );
  ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
  ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
  ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );
  ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);
  ((TH1F*)gDirectory->FindObject(h_angle_e1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()) );
  ((TH1F*)gDirectory->FindObject(h_angle_e1e2_mu->GetName()))->Fill( (*eelv).Vect().Angle(muchlv->Vect()) );
  
  ((TH1F*)gDirectory->FindObject(h_muchtrktime->GetName()))        ->Fill(muchtrktime);
  ((TH1F*)gDirectory->FindObject(h_muchtrkhodtime->GetName()))        ->Fill( muchtrkhodtime );
  ((TH1F*)gDirectory->FindObject(h_muchtrktimediff->GetName()))    ->Fill( muchtrktime-muchtrkhodtime );
  //((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( muchtrktime-muchcltime );
  
  ((TH1F*)gDirectory->FindObject(h_muchlvE->GetName()))->Fill( (*muchlv).E() );
  ((TH1F*)gDirectory->FindObject(h_muchlvP->GetName()))->Fill( (*muchlv).P() );
  ((TH1F*)gDirectory->FindObject(h_muchlvPt->GetName()))->Fill((*muchlv).Pt() );
  
  ((TH1F*)gDirectory->FindObject(h_KE->GetName()))->Fill((*pK).E());
  ((TH1F*)gDirectory->FindObject(h_KP->GetName()))->Fill((*pK).P());
  
  ((TH1F*)gDirectory->FindObject(h_nu_P->GetName()))->Fill(  (*neutrino).P() );
  ((TH1F*)gDirectory->FindObject(h_nu_E->GetName()))->Fill(  (*neutrino).E() );
  ((TH1F*)gDirectory->FindObject(h_nu_Pt->GetName()))->Fill( (*neutrino).Pt() );
  ((TH1F*)gDirectory->FindObject(h_nu_m->GetName()))->Fill(  (*neutrino).M() );
  ((TH1F*)gDirectory->FindObject(h_missing_mass->GetName()))->Fill( (*neutrino).M2() );
  
  
  ((TH1F*)gDirectory->FindObject(h_mueeE->GetName()))->Fill( (*muvee).E());
  ((TH1F*)gDirectory->FindObject(h_mueeP->GetName()))->Fill( (*muvee).P());
  ((TH1F*)gDirectory->FindObject(h_mueePt->GetName()))->Fill((*muvee).Pt());
  ((TH1F*)gDirectory->FindObject(h_mueeM->GetName()))->Fill ((*muvee).M());
  ((TH1F*)gDirectory->FindObject(h_3trk_invm->GetName()))->Fill ((*muvee).M2());
  
  
  ((TH1F*)gDirectory->FindObject(h_bx_much->GetName()))            ->Fill(muchp_uncorr[0]);
  ((TH1F*)gDirectory->FindObject(h_by_much->GetName()))            ->Fill(muchp_uncorr[1]);
  ((TH1F*)gDirectory->FindObject(h_bdxdz_much->GetName()))            ->Fill(muchv_uncorr[0]);
  ((TH1F*)gDirectory->FindObject(h_bdydz_much->GetName()))            ->Fill(muchv_uncorr[1]) ;
  ((TH1F*)gDirectory->FindObject(h_bx_vs_by_much->GetName()))        ->Fill(muchp_uncorr[0], muchp_uncorr[1]);
  ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_much->GetName()))->Fill(muchv_uncorr[0], muchv_uncorr[1]);
  
  ((TH1F*)gDirectory->FindObject(h_cda_ee->GetName()))->Fill(cda);
  ((TH1F*)gDirectory->FindObject(h_vx->GetName()))->Fill(vertex[0]);
  ((TH1F*)gDirectory->FindObject(h_vy->GetName()))->Fill(vertex[1]);
  ((TH1F*)gDirectory->FindObject(h_vz->GetName()))->Fill(vertex[2]);
  ((TH2F*)gDirectory->FindObject(h_vx_vs_vy->GetName()))->Fill(vertex[0],vertex[1]);
  
  ((TH1F*)gDirectory->FindObject(h_vx_e1much->GetName()))->Fill(vertex_e1much[0]);
  ((TH1F*)gDirectory->FindObject(h_vy_e1much->GetName()))->Fill(vertex_e1much[1]);
  ((TH1F*)gDirectory->FindObject(h_vz_e1much->GetName()))->Fill(vertex_e1much[2]);
  ((TH1F*)gDirectory->FindObject(h_cda_e1much->GetName()))->Fill(cda_e1much);
  ((TH1F*)gDirectory->FindObject(h_vz_e1much_diff->GetName()))->Fill( vertex[2] - vertex_e1much[2] );
  
  ((TH1F*)gDirectory->FindObject(h_vx_e2much->GetName()))->Fill(vertex_e2much[0]);
  ((TH1F*)gDirectory->FindObject(h_vy_e2much->GetName()))->Fill(vertex_e2much[1]);
  ((TH1F*)gDirectory->FindObject(h_vz_e2much->GetName()))->Fill(vertex_e2much[2]);
  ((TH1F*)gDirectory->FindObject(h_cda_e2much->GetName()))->Fill(cda_e2much);
  ((TH1F*)gDirectory->FindObject(h_vz_e2much_diff->GetName()))->Fill( vertex[2] - vertex_e2much[2] );
  ((TH1F*)gDirectory->FindObject(h_vz_mu1mu2_diff->GetName()))->Fill( vertex_e1much[2] - vertex_e2much[2] );
  
  
  ((TH1F*)gDirectory->FindObject(h_cda_e1much->GetName()))->Fill(cda_e1much);
  ((TH2F*)gDirectory->FindObject(h_vx_vs_vy_e1much->GetName()))->Fill(vertex_e1much[0],vertex_e1much[1]);
  
  ((TH1F*)gDirectory->FindObject(h_muche1distDHC1plane->GetName()))->Fill( muche1distDHC1plane );
  ((TH1F*)gDirectory->FindObject(h_muche2distDHC1plane->GetName()))->Fill( muche2distDHC1plane );
  ((TH1F*)gDirectory->FindObject(h_muche1distLKrplane->GetName()))->Fill(  muche1distLKrplane );
  ((TH1F*)gDirectory->FindObject(h_muche2distLKrplane->GetName()))->Fill(  muche2distLKrplane );
  
  ((TH1F*)gDirectory->FindObject(h_eldistDCH1plane->GetName()))->Fill( eldistDCH1plane );
  ((TH1F*)gDirectory->FindObject(h_eldistLKrplane->GetName()))->Fill( eldistLKrplane );
  
  file1->cd();
  h_cutflow->Fill(cutcounter);
  h_cutflow->GetXaxis()->SetBinLabel( cutcounter, cuts[cutcounter-1].c_str() );
  cutcounter++;
  //THE END OF ************* Geometry cuts ****************** ///

  // ************* Vertex matching ****************** ///
  
  //"e1e2_vtx_found"
  if(IS_DATA){
    if(vertex_OK!=1)    {cleanup();return 0;}//return 0;
    if( muchee_vertex_OK!=1)    {cleanup();return 0;}
    if( muchee2_vertex_OK!=1)    {cleanup();return 0;}
  
    //"e1e2_vtx_good" 
    //if (flux_id==1){
    //  if( vertex[2]< -1800. )    {cleanup();return 0;}
    //  if( vertex_e1much[2]< -1800. )    {cleanup();return 0;}
    //  if( vertex_e2much[2]< -1800. )    {cleanup();return 0;}  
    //} 
    //if (flux_id==0){
    if( vertex[2]<-1500. || vertex[2]>7000. )    {cleanup();return 0;}
    if( vertex_e1much[2]<-1500.  || vertex_e1much[2]>7000. )    {cleanup();return 0;}  
    if( vertex_e2much[2]<-1500.  || vertex_e2much[2]>7000. )    {cleanup();return 0;}  
    //}
    //"mu_vtx_good"
    if( fabs(vertex[2] - vertex_e2much[2]) > 100  )    {cleanup();return 0;}
    if( fabs(vertex[2] - vertex_e1much[2]) > 100  )    {cleanup();return 0;}
  
  }
  dirs[cutcounter-1]->cd();
  
  if (IS_MC){
    //if(true_4.P()  < 3. ||true_5.P() < 3.) {cleanup();return 0;}
    ((TH1D*)gDirectory->FindObject( h_mc_muvee_M->GetName()))->Fill( Muee_three_track.M());
    ((TH1D*)gDirectory->FindObject( h_mc_muvee_P->GetName()))->Fill( Muee_three_track.P());
    ((TH1D*)gDirectory->FindObject( h_mc_muvee_Pt->GetName()))->Fill( Muee_three_track.Pt());
    ((TH1D*)gDirectory->FindObject( h_mc_elsP->GetName()))->Fill( mc_ee.P());
    ((TH1D*)gDirectory->FindObject( h_mc_elsM->GetName()))->Fill( mc_ee.M());
    ((TH1D*)gDirectory->FindObject( h_mc_nu_M->GetName()))->Fill( Nu_true.M());
    ((TH1D*)gDirectory->FindObject( h_mc_nu_M2->GetName()))->Fill( Nu_true.M2());
    ((TH1D*)gDirectory->FindObject( h_mc_nu_P->GetName()))->Fill( Nu_true.P());
    ((TH1D*)gDirectory->FindObject( h_mc_nu_Pt->GetName()))->Fill(Nu_true.Pt());
  

    ((TH1D*)gDirectory->FindObject(h_true_M_0->GetName()))->Fill( true_0.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_0->GetName()))->Fill( true_0.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_0->GetName()))->Fill(true_0.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_0->GetName()))->Fill( true_0.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_1->GetName()))->Fill( true_1.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_1->GetName()))->Fill( true_1.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_1->GetName()))->Fill(true_1.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_1->GetName()))->Fill( true_1.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_2->GetName()))->Fill( true_2.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_2->GetName()))->Fill( true_2.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_2->GetName()))->Fill(true_2.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_2->GetName()))->Fill( true_2.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_3->GetName()))->Fill( true_3.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_3->GetName()))->Fill( true_3.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_3->GetName()))->Fill(true_3.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_3->GetName()))->Fill( true_3.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_4.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_4.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_4.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_4.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_5->GetName()))->Fill( true_5.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_5->GetName()))->Fill( true_5.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_5->GetName()))->Fill(true_5.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_5->GetName()))->Fill( true_5.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_6->GetName()))->Fill( true_6.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_6->GetName()))->Fill( true_6.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_6->GetName()))->Fill(true_6.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_6->GetName()))->Fill( true_6.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_7->GetName()))->Fill( true_7.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_7->GetName()))->Fill( true_7.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_7->GetName()))->Fill(true_7.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_7->GetName()))->Fill( true_7.E());
  }
  
  



  ((TH1D*)gDirectory->FindObject( h_vz_initial->GetName()))->Fill( sevt->vtx[0].z);
  ((TH1I*)gDirectory->FindObject(h_ntrack->GetName()))->Fill(ntrack);
  ((TH1I*)gDirectory->FindObject(h_nclust->GetName()))->Fill(nclust);
  
  ((TH1I*)gDirectory->FindObject(h_evttypebytrks->GetName()))->Fill(evttypebytrks);
  ((TH1I*)gDirectory->FindObject(h_nelectrack->GetName()))->Fill(nelectrack);
  ((TH1I*)gDirectory->FindObject(h_mutrack->GetName()))->Fill(nmutrack);
  ((TH2I*)gDirectory->FindObject(h_neltrk_mutrk->GetName()))->Fill(nelectrack, nmutrack);
  ((TH1I*)gDirectory->FindObject(h_mutrackq->GetName()))->Fill(muchtrkq);
  ((TH1I*)gDirectory->FindObject(h_e1trackq->GetName()))->Fill(eltrkq1);
  ((TH1I*)gDirectory->FindObject(h_e2trackq->GetName()))->Fill(eltrkq2);

  ((TH1F*)gDirectory->FindObject(h_el1trktime->GetName()))        ->Fill( el1trktime );
  ((TH1F*)gDirectory->FindObject(h_el1trkhodtime->GetName()))        ->Fill( el1trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_el1trktimediff->GetName()))    ->Fill( el1trktime - el1trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_ecltime1->GetName()))            ->Fill( el1cltime );
  ((TH1F*)gDirectory->FindObject(h_e1trkcltimediff->GetName()))    ->Fill( el1trktime- el1cltime );

  ((TH1F*)gDirectory->FindObject(h_el2trktime->GetName()))        ->Fill( el2trktime );
  ((TH1F*)gDirectory->FindObject(h_el2trkhodtime->GetName()))        ->Fill( el2trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_el2trktimediff->GetName()))    ->Fill( el2trktime - el2trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_ecltime2->GetName()))            ->Fill( el2cltime );
  ((TH1F*)gDirectory->FindObject(h_e2trkcltimediff->GetName()))    ->Fill( el2trktime-el2cltime );

  ((TH1F*)gDirectory->FindObject(h_e1e2timediff->GetName()))        ->Fill( el1trktime-el2trktime );
  ((TH1F*)gDirectory->FindObject(h_e1e2cltimediff->GetName()))    ->Fill( el1cltime-el2cltime );
  ((TH1F*)gDirectory->FindObject(h_e1e2hodtimediff->GetName()))        ->Fill( el1trkhodtime-el2trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_e1muhodtimediff->GetName()))        ->Fill( el1trkhodtime-muchtrkhodtime );
  ((TH1F*)gDirectory->FindObject(h_e2muhodtimediff->GetName()))        ->Fill( el2trkhodtime-muchtrkhodtime );


  ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_mu_track);
  ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_e1_track);
  ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_e2_track);

  
  ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_mu_track, (*muchlv).P());
  ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_e1_track, (*el1lv).P() );
  ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_e2_track, (*el2lv).P() );
 
  ((TH1F*)gDirectory->FindObject(h_eop_mu->GetName()))->Fill(eop_mu_track);
  ((TH2F*)gDirectory->FindObject(h_eop_mu_vs_p->GetName()))->Fill(eop_mu_track, (*muchlv).P());
  ((TH1F*)gDirectory->FindObject(h_eop_el->GetName()))->Fill(eop_e1_track);
  ((TH2F*)gDirectory->FindObject(h_eop_el_vs_p->GetName()))->Fill(eop_e1_track, (*el1lv).P());
  ((TH1F*)gDirectory->FindObject(h_eop_el->GetName()))->Fill(eop_e2_track);
  ((TH2F*)gDirectory->FindObject(h_eop_el_vs_p->GetName()))->Fill(eop_e2_track, (*el2lv).P());


  ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill( elp1_uncorr[0] );
  ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill( elp1_uncorr[1] );
  ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill( elv1_uncorr[0] );
  ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill( elv1_uncorr[1] );
  ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill( elp1_uncorr[0], elp1_uncorr[1] );
  ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill( elv1_uncorr[0], elv1_uncorr[1] );

  ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( el1lv->P() );
  ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( el2lv->P() );
  ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el1lv->E() );
  ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el2lv->E() );
  ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el1lv->P() );
  ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el2lv->P() );
  ((TH1F*)gDirectory->FindObject(h_e1m->GetName()))->Fill( el1lv->M() );
  ((TH1F*)gDirectory->FindObject(h_e2m->GetName()))->Fill( el2lv->M() );
  ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el1lv->Pt() );
  ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el2lv->Pt() );
  ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
  ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
  ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );
  ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);
  ((TH1F*)gDirectory->FindObject(h_angle_e1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()) );
  ((TH1F*)gDirectory->FindObject(h_angle_e1e2_mu->GetName()))->Fill( (*eelv).Vect().Angle(muchlv->Vect()) );

  ((TH1F*)gDirectory->FindObject(h_muchtrktime->GetName()))        ->Fill(muchtrktime);
  ((TH1F*)gDirectory->FindObject(h_muchtrkhodtime->GetName()))        ->Fill( muchtrkhodtime );
  ((TH1F*)gDirectory->FindObject(h_muchtrktimediff->GetName()))    ->Fill( muchtrktime-muchtrkhodtime );
  //((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( muchtrktime-muchcltime );

  ((TH1F*)gDirectory->FindObject(h_muchlvE->GetName()))->Fill( (*muchlv).E() );
  ((TH1F*)gDirectory->FindObject(h_muchlvP->GetName()))->Fill( (*muchlv).P() );
  ((TH1F*)gDirectory->FindObject(h_muchlvPt->GetName()))->Fill((*muchlv).Pt() );

  ((TH1F*)gDirectory->FindObject(h_KE->GetName()))->Fill((*pK).E());
  ((TH1F*)gDirectory->FindObject(h_KP->GetName()))->Fill((*pK).P());

  ((TH1F*)gDirectory->FindObject(h_nu_P->GetName()))->Fill(  (*neutrino).P() );
  ((TH1F*)gDirectory->FindObject(h_nu_E->GetName()))->Fill(  (*neutrino).E() );
  ((TH1F*)gDirectory->FindObject(h_nu_Pt->GetName()))->Fill( (*neutrino).Pt() );
  ((TH1F*)gDirectory->FindObject(h_nu_m->GetName()))->Fill(  (*neutrino).M() );
  ((TH1F*)gDirectory->FindObject(h_missing_mass->GetName()))->Fill( (*neutrino).M2() );


  ((TH1F*)gDirectory->FindObject(h_mueeE->GetName()))->Fill( (*muvee).E());
  ((TH1F*)gDirectory->FindObject(h_mueeP->GetName()))->Fill( (*muvee).P());
  ((TH1F*)gDirectory->FindObject(h_mueePt->GetName()))->Fill((*muvee).Pt());
  ((TH1F*)gDirectory->FindObject(h_mueeM->GetName()))->Fill ((*muvee).M());
  ((TH1F*)gDirectory->FindObject(h_3trk_invm->GetName()))->Fill ((*muvee).M2());


  ((TH1F*)gDirectory->FindObject(h_bx_much->GetName()))            ->Fill(muchp_uncorr[0]);
  ((TH1F*)gDirectory->FindObject(h_by_much->GetName()))            ->Fill(muchp_uncorr[1]);
  ((TH1F*)gDirectory->FindObject(h_bdxdz_much->GetName()))            ->Fill(muchv_uncorr[0]);
  ((TH1F*)gDirectory->FindObject(h_bdydz_much->GetName()))            ->Fill(muchv_uncorr[1]) ;
  ((TH1F*)gDirectory->FindObject(h_bx_vs_by_much->GetName()))        ->Fill(muchp_uncorr[0], muchp_uncorr[1]);
  ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_much->GetName()))->Fill(muchv_uncorr[0], muchv_uncorr[1]);

  ((TH1F*)gDirectory->FindObject(h_cda_ee->GetName()))->Fill(cda);
  ((TH1F*)gDirectory->FindObject(h_vx->GetName()))->Fill(vertex[0]);
  ((TH1F*)gDirectory->FindObject(h_vy->GetName()))->Fill(vertex[1]);
  ((TH1F*)gDirectory->FindObject(h_vz->GetName()))->Fill(vertex[2]);
  ((TH2F*)gDirectory->FindObject(h_vx_vs_vy->GetName()))->Fill(vertex[0],vertex[1]);

  ((TH1F*)gDirectory->FindObject(h_vx_e1much->GetName()))->Fill(vertex_e1much[0]);
  ((TH1F*)gDirectory->FindObject(h_vy_e1much->GetName()))->Fill(vertex_e1much[1]);
  ((TH1F*)gDirectory->FindObject(h_vz_e1much->GetName()))->Fill(vertex_e1much[2]);
  ((TH1F*)gDirectory->FindObject(h_cda_e1much->GetName()))->Fill(cda_e1much);
  ((TH1F*)gDirectory->FindObject(h_vz_e1much_diff->GetName()))->Fill( vertex[2] - vertex_e1much[2] );

  ((TH1F*)gDirectory->FindObject(h_vx_e2much->GetName()))->Fill(vertex_e2much[0]);
  ((TH1F*)gDirectory->FindObject(h_vy_e2much->GetName()))->Fill(vertex_e2much[1]);
  ((TH1F*)gDirectory->FindObject(h_vz_e2much->GetName()))->Fill(vertex_e2much[2]);
  ((TH1F*)gDirectory->FindObject(h_cda_e2much->GetName()))->Fill(cda_e2much);
  ((TH1F*)gDirectory->FindObject(h_vz_e2much_diff->GetName()))->Fill( vertex[2] - vertex_e2much[2] );
  ((TH1F*)gDirectory->FindObject(h_vz_mu1mu2_diff->GetName()))->Fill( vertex_e1much[2] - vertex_e2much[2] );


  ((TH1F*)gDirectory->FindObject(h_cda_e1much->GetName()))->Fill(cda_e1much);
  ((TH2F*)gDirectory->FindObject(h_vx_vs_vy_e1much->GetName()))->Fill(vertex_e1much[0],vertex_e1much[1]);

  ((TH1F*)gDirectory->FindObject(h_muche1distDHC1plane->GetName()))->Fill( muche1distDHC1plane );
  ((TH1F*)gDirectory->FindObject(h_muche2distDHC1plane->GetName()))->Fill( muche2distDHC1plane );
  ((TH1F*)gDirectory->FindObject(h_muche1distLKrplane->GetName()))->Fill(  muche1distLKrplane );
  ((TH1F*)gDirectory->FindObject(h_muche2distLKrplane->GetName()))->Fill(  muche2distLKrplane );
  ((TH1F*)gDirectory->FindObject(h_eldistDCH1plane->GetName()))->Fill( eldistDCH1plane );
  ((TH1F*)gDirectory->FindObject(h_eldistLKrplane->GetName()))->Fill( eldistLKrplane );

  file1->cd();
  h_cutflow->Fill(cutcounter);
  h_cutflow->GetXaxis()->SetBinLabel( cutcounter, cuts[cutcounter-1].c_str() );
  cutcounter++;
  ///THE END OF ************* Vertex matching ****************** ///

  /// ************* Time matching ****************** ///
  //"e1e2_hod_Tmatch"
  if(IS_DATA)    {
    if ( el1cltime-el2cltime < -3 || el1cltime-el2cltime > 3 ) {cleanup();return 0;}
    if (  el1trkhodtime-el2trkhodtime < -2 || el1trkhodtime-el2trkhodtime > 2 ) {cleanup();return 0;}
    if ( el1trkhodtime-muchtrkhodtime  < -2 || el1trkhodtime - muchtrkhodtime > 2 ) {cleanup();return 0;}
    if ( el2trkhodtime-muchtrkhodtime  < -2 || el2trkhodtime - muchtrkhodtime > 2 ) {cleanup();return 0;}
    if ( el1trktime-el1cltime > 10) {cleanup();return 0;}
    if ( el2trktime-el1cltime > 10) {cleanup();return 0;}
 
    //5ns is according to Manuel's thesis
    if( fabs(el1trktime - el1trkhodtime)>12. )     return 0;
    if( fabs(el2trktime - el2trkhodtime)>12. )     return 0;
    if( fabs(muchtrktime-muchtrkhodtime) > 12. )    return 0;  //5 is from Manuel's thesis
    if( el1trkhodtime<136. || el1trkhodtime>146. )     {cleanup();return 0;}
    if( el2trkhodtime<136. || el2trkhodtime>146. )     {cleanup();return 0;}
    if( muchtrkhodtime<135. || muchtrkhodtime>148. )    return 0;
    if( fabs((long double)(el1trktime-muchtrktime)) > 10. )    {cleanup();return 0;}
    if( fabs((long double)(el2trktime-muchtrktime)) > 10. )    {cleanup();return 0;}
  
  }
 
  //if( fabs((long double)(el1cltime-muchcltime))   > 10. )    {cleanup();return 0;}

  dirs[cutcounter-1]->cd();
    
  if (IS_MC){
    //if(true_4.P()  < 3. ||true_5.P() < 3.) {cleanup();return 0;}
    ((TH1D*)gDirectory->FindObject( h_mc_muvee_M->GetName()))->Fill( Muee_three_track.M());
    ((TH1D*)gDirectory->FindObject( h_mc_muvee_P->GetName()))->Fill( Muee_three_track.P());
    ((TH1D*)gDirectory->FindObject( h_mc_muvee_Pt->GetName()))->Fill( Muee_three_track.Pt());
    ((TH1D*)gDirectory->FindObject( h_mc_elsP->GetName()))->Fill( mc_ee.P());
    ((TH1D*)gDirectory->FindObject( h_mc_elsM->GetName()))->Fill( mc_ee.M());
    ((TH1D*)gDirectory->FindObject( h_mc_nu_M->GetName()))->Fill( Nu_true.M());
    ((TH1D*)gDirectory->FindObject( h_mc_nu_M2->GetName()))->Fill( Nu_true.M2());
    ((TH1D*)gDirectory->FindObject( h_mc_nu_P->GetName()))->Fill( Nu_true.P());
    ((TH1D*)gDirectory->FindObject( h_mc_nu_Pt->GetName()))->Fill(Nu_true.Pt());
  
    ((TH1D*)gDirectory->FindObject(h_true_M_0->GetName()))->Fill( true_0.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_0->GetName()))->Fill( true_0.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_0->GetName()))->Fill(true_0.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_0->GetName()))->Fill( true_0.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_1->GetName()))->Fill( true_1.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_1->GetName()))->Fill( true_1.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_1->GetName()))->Fill(true_1.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_1->GetName()))->Fill( true_1.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_2->GetName()))->Fill( true_2.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_2->GetName()))->Fill( true_2.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_2->GetName()))->Fill(true_2.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_2->GetName()))->Fill( true_2.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_3->GetName()))->Fill( true_3.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_3->GetName()))->Fill( true_3.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_3->GetName()))->Fill(true_3.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_3->GetName()))->Fill( true_3.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_4.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_4.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_4.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_4.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_5->GetName()))->Fill( true_5.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_5->GetName()))->Fill( true_5.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_5->GetName()))->Fill(true_5.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_5->GetName()))->Fill( true_5.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_6->GetName()))->Fill( true_6.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_6->GetName()))->Fill( true_6.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_6->GetName()))->Fill(true_6.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_6->GetName()))->Fill( true_6.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_7->GetName()))->Fill( true_7.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_7->GetName()))->Fill( true_7.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_7->GetName()))->Fill(true_7.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_7->GetName()))->Fill( true_7.E());
  
  }


((TH1D*)gDirectory->FindObject( h_vz_initial->GetName()))->Fill( sevt->vtx[0].z);
  ((TH1I*)gDirectory->FindObject(h_ntrack->GetName()))->Fill(ntrack);
  ((TH1I*)gDirectory->FindObject(h_nclust->GetName()))->Fill(nclust);

  ((TH1I*)gDirectory->FindObject(h_evttypebytrks->GetName()))->Fill(evttypebytrks);
  ((TH1I*)gDirectory->FindObject(h_nelectrack->GetName()))->Fill(nelectrack);
  ((TH1I*)gDirectory->FindObject(h_mutrack->GetName()))->Fill(nmutrack);
  ((TH2I*)gDirectory->FindObject(h_neltrk_mutrk->GetName()))->Fill(nelectrack, nmutrack);
  ((TH1I*)gDirectory->FindObject(h_mutrackq->GetName()))->Fill(muchtrkq);
  ((TH1I*)gDirectory->FindObject(h_e1trackq->GetName()))->Fill(eltrkq1);
  ((TH1I*)gDirectory->FindObject(h_e2trackq->GetName()))->Fill(eltrkq2);

  ((TH1F*)gDirectory->FindObject(h_el1trktime->GetName()))        ->Fill( el1trktime );
  ((TH1F*)gDirectory->FindObject(h_el1trkhodtime->GetName()))        ->Fill( el1trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_el1trktimediff->GetName()))    ->Fill( el1trktime - el1trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_ecltime1->GetName()))            ->Fill( el1cltime );
  ((TH1F*)gDirectory->FindObject(h_e1trkcltimediff->GetName()))    ->Fill( el1trktime- el1cltime );

  ((TH1F*)gDirectory->FindObject(h_el2trktime->GetName()))        ->Fill( el2trktime );
  ((TH1F*)gDirectory->FindObject(h_el2trkhodtime->GetName()))        ->Fill( el2trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_el2trktimediff->GetName()))    ->Fill( el2trktime - el2trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_ecltime2->GetName()))            ->Fill( el2cltime );
  ((TH1F*)gDirectory->FindObject(h_e2trkcltimediff->GetName()))    ->Fill( el2trktime-el2cltime );

  ((TH1F*)gDirectory->FindObject(h_e1e2timediff->GetName()))        ->Fill( el1trktime-el2trktime );
  ((TH1F*)gDirectory->FindObject(h_e1e2cltimediff->GetName()))    ->Fill( el1cltime-el2cltime );
  ((TH1F*)gDirectory->FindObject(h_e1e2hodtimediff->GetName()))        ->Fill( el1trkhodtime-el2trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_e1muhodtimediff->GetName()))        ->Fill( el1trkhodtime-muchtrkhodtime );
  ((TH1F*)gDirectory->FindObject(h_e2muhodtimediff->GetName()))        ->Fill( el2trkhodtime-muchtrkhodtime );



  ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_mu_track);
  ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_e1_track);
  ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_e2_track);

  
  ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_mu_track, (*muchlv).P());
  ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_e1_track, (*el1lv).P() );
  ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_e2_track, (*el2lv).P() );
 
  ((TH1F*)gDirectory->FindObject(h_eop_mu->GetName()))->Fill(eop_mu_track);
  ((TH2F*)gDirectory->FindObject(h_eop_mu_vs_p->GetName()))->Fill(eop_mu_track, (*muchlv).P());
  ((TH1F*)gDirectory->FindObject(h_eop_el->GetName()))->Fill(eop_e1_track);
  ((TH2F*)gDirectory->FindObject(h_eop_el_vs_p->GetName()))->Fill(eop_e1_track, (*el1lv).P());
  ((TH1F*)gDirectory->FindObject(h_eop_el->GetName()))->Fill(eop_e2_track);
  ((TH2F*)gDirectory->FindObject(h_eop_el_vs_p->GetName()))->Fill(eop_e2_track, (*el2lv).P());


 
  ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill( elp1_uncorr[0] );
  ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill( elp1_uncorr[1] );
  ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill( elv1_uncorr[0] );
  ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill( elv1_uncorr[1] );
  ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill( elp1_uncorr[0], elp1_uncorr[1] );
  ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill( elv1_uncorr[0], elv1_uncorr[1] );

  ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( el1lv->P() );
  ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( el2lv->P() );
  ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el1lv->E() );
  ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el2lv->E() );
  ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el1lv->P() );
  ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el2lv->P() );
  ((TH1F*)gDirectory->FindObject(h_e1m->GetName()))->Fill( el1lv->M() );
  ((TH1F*)gDirectory->FindObject(h_e2m->GetName()))->Fill( el2lv->M() );
  ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el1lv->Pt() );
  ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el2lv->Pt() );
  ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
  ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
  ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );
  ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);
  ((TH1F*)gDirectory->FindObject(h_angle_e1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()) );
  ((TH1F*)gDirectory->FindObject(h_angle_e1e2_mu->GetName()))->Fill( (*eelv).Vect().Angle(muchlv->Vect()) );

  ((TH1F*)gDirectory->FindObject(h_muchtrktime->GetName()))        ->Fill(muchtrktime);
  ((TH1F*)gDirectory->FindObject(h_muchtrkhodtime->GetName()))        ->Fill( muchtrkhodtime );
  ((TH1F*)gDirectory->FindObject(h_muchtrktimediff->GetName()))    ->Fill( muchtrktime-muchtrkhodtime );
  //((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( muchtrktime-muchcltime );

  ((TH1F*)gDirectory->FindObject(h_muchlvE->GetName()))->Fill( (*muchlv).E() );
  ((TH1F*)gDirectory->FindObject(h_muchlvP->GetName()))->Fill( (*muchlv).P() );
  ((TH1F*)gDirectory->FindObject(h_muchlvPt->GetName()))->Fill((*muchlv).Pt() );

  ((TH1F*)gDirectory->FindObject(h_KE->GetName()))->Fill((*pK).E());
  ((TH1F*)gDirectory->FindObject(h_KP->GetName()))->Fill((*pK).P());

  ((TH1F*)gDirectory->FindObject(h_nu_P->GetName()))->Fill(  (*neutrino).P() );
  ((TH1F*)gDirectory->FindObject(h_nu_E->GetName()))->Fill(  (*neutrino).E() );
  ((TH1F*)gDirectory->FindObject(h_nu_Pt->GetName()))->Fill( (*neutrino).Pt() );
  ((TH1F*)gDirectory->FindObject(h_nu_m->GetName()))->Fill(  (*neutrino).M() );
  ((TH1F*)gDirectory->FindObject(h_missing_mass->GetName()))->Fill( (*neutrino).M2() );


  ((TH1F*)gDirectory->FindObject(h_mueeE->GetName()))->Fill( (*muvee).E());
  ((TH1F*)gDirectory->FindObject(h_mueeP->GetName()))->Fill( (*muvee).P());
  ((TH1F*)gDirectory->FindObject(h_mueePt->GetName()))->Fill((*muvee).Pt());
  ((TH1F*)gDirectory->FindObject(h_mueeM->GetName()))->Fill ((*muvee).M());
  ((TH1F*)gDirectory->FindObject(h_3trk_invm->GetName()))->Fill ((*muvee).M2());


  ((TH1F*)gDirectory->FindObject(h_bx_much->GetName()))            ->Fill(muchp_uncorr[0]);
  ((TH1F*)gDirectory->FindObject(h_by_much->GetName()))            ->Fill(muchp_uncorr[1]);
  ((TH1F*)gDirectory->FindObject(h_bdxdz_much->GetName()))            ->Fill(muchv_uncorr[0]);
  ((TH1F*)gDirectory->FindObject(h_bdydz_much->GetName()))            ->Fill(muchv_uncorr[1]) ;
  ((TH1F*)gDirectory->FindObject(h_bx_vs_by_much->GetName()))        ->Fill(muchp_uncorr[0], muchp_uncorr[1]);
  ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_much->GetName()))->Fill(muchv_uncorr[0], muchv_uncorr[1]);

  ((TH1F*)gDirectory->FindObject(h_cda_ee->GetName()))->Fill(cda);
  ((TH1F*)gDirectory->FindObject(h_vx->GetName()))->Fill(vertex[0]);
  ((TH1F*)gDirectory->FindObject(h_vy->GetName()))->Fill(vertex[1]);
  ((TH1F*)gDirectory->FindObject(h_vz->GetName()))->Fill(vertex[2]);
  ((TH2F*)gDirectory->FindObject(h_vx_vs_vy->GetName()))->Fill(vertex[0],vertex[1]);

  ((TH1F*)gDirectory->FindObject(h_vx_e1much->GetName()))->Fill(vertex_e1much[0]);
  ((TH1F*)gDirectory->FindObject(h_vy_e1much->GetName()))->Fill(vertex_e1much[1]);
  ((TH1F*)gDirectory->FindObject(h_vz_e1much->GetName()))->Fill(vertex_e1much[2]);
  ((TH1F*)gDirectory->FindObject(h_cda_e1much->GetName()))->Fill(cda_e1much);
  ((TH1F*)gDirectory->FindObject(h_vz_e1much_diff->GetName()))->Fill( vertex[2] - vertex_e1much[2] );

  ((TH1F*)gDirectory->FindObject(h_vx_e2much->GetName()))->Fill(vertex_e2much[0]);
  ((TH1F*)gDirectory->FindObject(h_vy_e2much->GetName()))->Fill(vertex_e2much[1]);
  ((TH1F*)gDirectory->FindObject(h_vz_e2much->GetName()))->Fill(vertex_e2much[2]);
  ((TH1F*)gDirectory->FindObject(h_cda_e2much->GetName()))->Fill(cda_e2much);
  ((TH1F*)gDirectory->FindObject(h_vz_e2much_diff->GetName()))->Fill( vertex[2] - vertex_e2much[2] );
  ((TH1F*)gDirectory->FindObject(h_vz_mu1mu2_diff->GetName()))->Fill( vertex_e1much[2] - vertex_e2much[2] );


  ((TH1F*)gDirectory->FindObject(h_cda_e1much->GetName()))->Fill(cda_e1much);
  ((TH2F*)gDirectory->FindObject(h_vx_vs_vy_e1much->GetName()))->Fill(vertex_e1much[0],vertex_e1much[1]);



  file1->cd();
  h_cutflow->Fill(cutcounter);
  h_cutflow->GetXaxis()->SetBinLabel( cutcounter, cuts[cutcounter-1].c_str() );
  cutcounter++;
  ///THE END OF ************* Time matching ****************** ///

  /// ************* Momentum cuts ****************** ///
  //if(flux_id==1){
  //  ;
  //}
  //if (flux_id==0){
  if ( (*el1lv).P() < 3. || (*el2lv).P() < 3. || (*muchlv).P() < 10.) {cleanup();return 0;}
  //}
  dirs[cutcounter-1]->cd();
  
  if (IS_MC){
    //if(true_4.P()  < 3. ||true_5.P() < 3.) {cleanup();return 0;}
    ((TH1D*)gDirectory->FindObject( h_mc_muvee_M->GetName()))->Fill( Muee_three_track.M());
    ((TH1D*)gDirectory->FindObject( h_mc_muvee_P->GetName()))->Fill( Muee_three_track.P());
    ((TH1D*)gDirectory->FindObject( h_mc_muvee_Pt->GetName()))->Fill( Muee_three_track.Pt());
    ((TH1D*)gDirectory->FindObject( h_mc_elsP->GetName()))->Fill( mc_ee.P());
    ((TH1D*)gDirectory->FindObject( h_mc_elsM->GetName()))->Fill( mc_ee.M());
    ((TH1D*)gDirectory->FindObject( h_mc_nu_M->GetName()))->Fill( Nu_true.M());
    ((TH1D*)gDirectory->FindObject( h_mc_nu_M2->GetName()))->Fill( Nu_true.M2());
    ((TH1D*)gDirectory->FindObject( h_mc_nu_P->GetName()))->Fill( Nu_true.P());
    ((TH1D*)gDirectory->FindObject( h_mc_nu_Pt->GetName()))->Fill(Nu_true.Pt());
  
    ((TH1D*)gDirectory->FindObject(h_true_M_0->GetName()))->Fill( true_0.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_0->GetName()))->Fill( true_0.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_0->GetName()))->Fill(true_0.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_0->GetName()))->Fill( true_0.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_1->GetName()))->Fill( true_1.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_1->GetName()))->Fill( true_1.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_1->GetName()))->Fill(true_1.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_1->GetName()))->Fill( true_1.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_2->GetName()))->Fill( true_2.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_2->GetName()))->Fill( true_2.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_2->GetName()))->Fill(true_2.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_2->GetName()))->Fill( true_2.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_3->GetName()))->Fill( true_3.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_3->GetName()))->Fill( true_3.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_3->GetName()))->Fill(true_3.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_3->GetName()))->Fill( true_3.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_4.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_4.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_4.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_4.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_5->GetName()))->Fill( true_5.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_5->GetName()))->Fill( true_5.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_5->GetName()))->Fill(true_5.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_5->GetName()))->Fill( true_5.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_6->GetName()))->Fill( true_6.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_6->GetName()))->Fill( true_6.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_6->GetName()))->Fill(true_6.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_6->GetName()))->Fill( true_6.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_7->GetName()))->Fill( true_7.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_7->GetName()))->Fill( true_7.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_7->GetName()))->Fill(true_7.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_7->GetName()))->Fill( true_7.E());
  
  }
  ((TH1D*)gDirectory->FindObject( h_vz_initial->GetName()))->Fill( sevt->vtx[0].z);
  ((TH1I*)gDirectory->FindObject(h_ntrack->GetName()))->Fill(ntrack);
  ((TH1I*)gDirectory->FindObject(h_nclust->GetName()))->Fill(nclust);
  
  ((TH1I*)gDirectory->FindObject(h_evttypebytrks->GetName()))->Fill(evttypebytrks);
  ((TH1I*)gDirectory->FindObject(h_nelectrack->GetName()))->Fill(nelectrack);
  ((TH1I*)gDirectory->FindObject(h_mutrack->GetName()))->Fill(nmutrack);
  ((TH2I*)gDirectory->FindObject(h_neltrk_mutrk->GetName()))->Fill(nelectrack, nmutrack);
  ((TH1I*)gDirectory->FindObject(h_mutrackq->GetName()))->Fill(muchtrkq);
  ((TH1I*)gDirectory->FindObject(h_e1trackq->GetName()))->Fill(eltrkq1);
  ((TH1I*)gDirectory->FindObject(h_e2trackq->GetName()))->Fill(eltrkq2);

  ((TH1F*)gDirectory->FindObject(h_el1trktime->GetName()))        ->Fill( el1trktime );
  ((TH1F*)gDirectory->FindObject(h_el1trkhodtime->GetName()))        ->Fill( el1trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_el1trktimediff->GetName()))    ->Fill( el1trktime - el1trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_ecltime1->GetName()))            ->Fill( el1cltime );
  ((TH1F*)gDirectory->FindObject(h_e1trkcltimediff->GetName()))    ->Fill( el1trktime- el1cltime );

  ((TH1F*)gDirectory->FindObject(h_el2trktime->GetName()))        ->Fill( el2trktime );
  ((TH1F*)gDirectory->FindObject(h_el2trkhodtime->GetName()))        ->Fill( el2trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_el2trktimediff->GetName()))    ->Fill( el2trktime - el2trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_ecltime2->GetName()))            ->Fill( el2cltime );
  ((TH1F*)gDirectory->FindObject(h_e2trkcltimediff->GetName()))    ->Fill( el2trktime-el2cltime );

  ((TH1F*)gDirectory->FindObject(h_e1e2timediff->GetName()))        ->Fill( el1trktime-el2trktime );
  ((TH1F*)gDirectory->FindObject(h_e1e2cltimediff->GetName()))    ->Fill( el1cltime-el2cltime );
  ((TH1F*)gDirectory->FindObject(h_e1e2hodtimediff->GetName()))        ->Fill( el1trkhodtime-el2trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_e1muhodtimediff->GetName()))        ->Fill( el1trkhodtime-muchtrkhodtime );
  ((TH1F*)gDirectory->FindObject(h_e2muhodtimediff->GetName()))        ->Fill( el2trkhodtime-muchtrkhodtime );



  ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_mu_track);
  ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_e1_track);
  ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_e2_track);

  
  ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_mu_track, (*muchlv).P());
  ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_e1_track, (*el1lv).P() );
  ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_e2_track, (*el2lv).P() );
 
  ((TH1F*)gDirectory->FindObject(h_eop_mu->GetName()))->Fill(eop_mu_track);
  ((TH2F*)gDirectory->FindObject(h_eop_mu_vs_p->GetName()))->Fill(eop_mu_track, (*muchlv).P());
  ((TH1F*)gDirectory->FindObject(h_eop_el->GetName()))->Fill(eop_e1_track);
  ((TH2F*)gDirectory->FindObject(h_eop_el_vs_p->GetName()))->Fill(eop_e1_track, (*el1lv).P());
  ((TH1F*)gDirectory->FindObject(h_eop_el->GetName()))->Fill(eop_e2_track);
  ((TH2F*)gDirectory->FindObject(h_eop_el_vs_p->GetName()))->Fill(eop_e2_track, (*el2lv).P());

  ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill( elp1_uncorr[0] );
  ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill( elp1_uncorr[1] );
  ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill( elv1_uncorr[0] );
  ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill( elv1_uncorr[1] );
  ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill( elp1_uncorr[0], elp1_uncorr[1] );
  ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill( elv1_uncorr[0], elv1_uncorr[1] );

  ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( el1lv->P() );
  ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( el2lv->P() );
  ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el1lv->E() );
  ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el2lv->E() );
  ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el1lv->P() );
  ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el2lv->P() );
  ((TH1F*)gDirectory->FindObject(h_e1m->GetName()))->Fill( el1lv->M() );
  ((TH1F*)gDirectory->FindObject(h_e2m->GetName()))->Fill( el2lv->M() );
  ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el1lv->Pt() );
  ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el2lv->Pt() );
  ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
  ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
  ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );
  ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);
  ((TH1F*)gDirectory->FindObject(h_angle_e1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()) );
  ((TH1F*)gDirectory->FindObject(h_angle_e1e2_mu->GetName()))->Fill( (*eelv).Vect().Angle(muchlv->Vect()) );

  ((TH1F*)gDirectory->FindObject(h_muchtrktime->GetName()))        ->Fill(muchtrktime);
  ((TH1F*)gDirectory->FindObject(h_muchtrkhodtime->GetName()))        ->Fill( muchtrkhodtime );
  ((TH1F*)gDirectory->FindObject(h_muchtrktimediff->GetName()))    ->Fill( muchtrktime-muchtrkhodtime );
  //((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( muchtrktime-muchcltime );

  ((TH1F*)gDirectory->FindObject(h_muchlvE->GetName()))->Fill( (*muchlv).E() );
  ((TH1F*)gDirectory->FindObject(h_muchlvP->GetName()))->Fill( (*muchlv).P() );
  ((TH1F*)gDirectory->FindObject(h_muchlvPt->GetName()))->Fill((*muchlv).Pt() );

  ((TH1F*)gDirectory->FindObject(h_KE->GetName()))->Fill((*pK).E());
  ((TH1F*)gDirectory->FindObject(h_KP->GetName()))->Fill((*pK).P());

  ((TH1F*)gDirectory->FindObject(h_nu_P->GetName()))->Fill(  (*neutrino).P() );
  ((TH1F*)gDirectory->FindObject(h_nu_E->GetName()))->Fill(  (*neutrino).E() );
  ((TH1F*)gDirectory->FindObject(h_nu_Pt->GetName()))->Fill( (*neutrino).Pt() );
  ((TH1F*)gDirectory->FindObject(h_nu_m->GetName()))->Fill(  (*neutrino).M() );
  ((TH1F*)gDirectory->FindObject(h_missing_mass->GetName()))->Fill( (*neutrino).M2() );


  ((TH1F*)gDirectory->FindObject(h_mueeE->GetName()))->Fill( (*muvee).E());
  ((TH1F*)gDirectory->FindObject(h_mueeP->GetName()))->Fill( (*muvee).P());
  ((TH1F*)gDirectory->FindObject(h_mueePt->GetName()))->Fill((*muvee).Pt());
  ((TH1F*)gDirectory->FindObject(h_mueeM->GetName()))->Fill ((*muvee).M());
  ((TH1F*)gDirectory->FindObject(h_3trk_invm->GetName()))->Fill ((*muvee).M2());


  ((TH1F*)gDirectory->FindObject(h_bx_much->GetName()))            ->Fill(muchp_uncorr[0]);
  ((TH1F*)gDirectory->FindObject(h_by_much->GetName()))            ->Fill(muchp_uncorr[1]);
  ((TH1F*)gDirectory->FindObject(h_bdxdz_much->GetName()))            ->Fill(muchv_uncorr[0]);
  ((TH1F*)gDirectory->FindObject(h_bdydz_much->GetName()))            ->Fill(muchv_uncorr[1]) ;
  ((TH1F*)gDirectory->FindObject(h_bx_vs_by_much->GetName()))        ->Fill(muchp_uncorr[0], muchp_uncorr[1]);
  ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_much->GetName()))->Fill(muchv_uncorr[0], muchv_uncorr[1]);

  ((TH1F*)gDirectory->FindObject(h_cda_ee->GetName()))->Fill(cda);
  ((TH1F*)gDirectory->FindObject(h_vx->GetName()))->Fill(vertex[0]);
  ((TH1F*)gDirectory->FindObject(h_vy->GetName()))->Fill(vertex[1]);
  ((TH1F*)gDirectory->FindObject(h_vz->GetName()))->Fill(vertex[2]);
  ((TH2F*)gDirectory->FindObject(h_vx_vs_vy->GetName()))->Fill(vertex[0],vertex[1]);

  ((TH1F*)gDirectory->FindObject(h_vx_e1much->GetName()))->Fill(vertex_e1much[0]);
  ((TH1F*)gDirectory->FindObject(h_vy_e1much->GetName()))->Fill(vertex_e1much[1]);
  ((TH1F*)gDirectory->FindObject(h_vz_e1much->GetName()))->Fill(vertex_e1much[2]);
  ((TH1F*)gDirectory->FindObject(h_cda_e1much->GetName()))->Fill(cda_e1much);
  ((TH1F*)gDirectory->FindObject(h_vz_e1much_diff->GetName()))->Fill( vertex[2] - vertex_e1much[2] );

  ((TH1F*)gDirectory->FindObject(h_vx_e2much->GetName()))->Fill(vertex_e2much[0]);
  ((TH1F*)gDirectory->FindObject(h_vy_e2much->GetName()))->Fill(vertex_e2much[1]);
  ((TH1F*)gDirectory->FindObject(h_vz_e2much->GetName()))->Fill(vertex_e2much[2]);
  ((TH1F*)gDirectory->FindObject(h_cda_e2much->GetName()))->Fill(cda_e2much);
  ((TH1F*)gDirectory->FindObject(h_vz_e2much_diff->GetName()))->Fill( vertex[2] - vertex_e2much[2] );
  ((TH1F*)gDirectory->FindObject(h_vz_mu1mu2_diff->GetName()))->Fill( vertex_e1much[2] - vertex_e2much[2] );


  ((TH1F*)gDirectory->FindObject(h_cda_e1much->GetName()))->Fill(cda_e1much);
  ((TH2F*)gDirectory->FindObject(h_vx_vs_vy_e1much->GetName()))->Fill(vertex_e1much[0],vertex_e1much[1]);

  file1->cd();
  h_cutflow->Fill(cutcounter);
  h_cutflow->GetXaxis()->SetBinLabel( cutcounter, cuts[cutcounter-1].c_str() );
  cutcounter++;
  ///THE END OF ************* Momentum cuts ****************** ///

  /// ************* K3pi selection ****************** ///
  dirs[cutcounter-1]->cd();
  if (nelectrack==2 && muchtrkindx!=-1) {
    if (eltrkq1==eltrkq2 &&
	eltrkq1!=muchtrkq   &&
	muchtrkq==1         &&
	eelvmass >= 0.140   )    {
      
      if (IS_MC){
	//if(true_4.P()  < 3. ||true_5.P() < 3.) {cleanup();return 0;}
	((TH1D*)gDirectory->FindObject( h_mc_muvee_M->GetName()))->Fill( Muee_three_track.M());
	((TH1D*)gDirectory->FindObject( h_mc_muvee_P->GetName()))->Fill( Muee_three_track.P());
	((TH1D*)gDirectory->FindObject( h_mc_muvee_Pt->GetName()))->Fill( Muee_three_track.Pt());
	((TH1D*)gDirectory->FindObject( h_mc_elsP->GetName()))->Fill( mc_ee.P());
	((TH1D*)gDirectory->FindObject( h_mc_elsM->GetName()))->Fill( mc_ee.M());
	((TH1D*)gDirectory->FindObject( h_mc_nu_M->GetName()))->Fill( Nu_true.M());
	((TH1D*)gDirectory->FindObject( h_mc_nu_M2->GetName()))->Fill( Nu_true.M2());
	((TH1D*)gDirectory->FindObject( h_mc_nu_P->GetName()))->Fill( Nu_true.P());
	((TH1D*)gDirectory->FindObject( h_mc_nu_Pt->GetName()))->Fill(Nu_true.Pt());
  
	((TH1D*)gDirectory->FindObject(h_true_M_0->GetName()))->Fill( true_0.M());
	((TH1D*)gDirectory->FindObject(h_true_P_0->GetName()))->Fill( true_0.P());
	((TH1D*)gDirectory->FindObject(h_true_Pt_0->GetName()))->Fill(true_0.Pt());
	((TH1D*)gDirectory->FindObject(h_true_E_0->GetName()))->Fill( true_0.E());
	((TH1D*)gDirectory->FindObject(h_true_M_1->GetName()))->Fill( true_1.M());
	((TH1D*)gDirectory->FindObject(h_true_P_1->GetName()))->Fill( true_1.P());
	((TH1D*)gDirectory->FindObject(h_true_Pt_1->GetName()))->Fill(true_1.Pt());
	((TH1D*)gDirectory->FindObject(h_true_E_1->GetName()))->Fill( true_1.E());
	((TH1D*)gDirectory->FindObject(h_true_M_2->GetName()))->Fill( true_2.M());
	((TH1D*)gDirectory->FindObject(h_true_P_2->GetName()))->Fill( true_2.P());
	((TH1D*)gDirectory->FindObject(h_true_Pt_2->GetName()))->Fill(true_2.Pt());
	((TH1D*)gDirectory->FindObject(h_true_E_2->GetName()))->Fill( true_2.E());
	((TH1D*)gDirectory->FindObject(h_true_M_3->GetName()))->Fill( true_3.M());
	((TH1D*)gDirectory->FindObject(h_true_P_3->GetName()))->Fill( true_3.P());
	((TH1D*)gDirectory->FindObject(h_true_Pt_3->GetName()))->Fill(true_3.Pt());
	((TH1D*)gDirectory->FindObject(h_true_E_3->GetName()))->Fill( true_3.E());
	((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_4.M());
	((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_4.P());
	((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_4.Pt());
	((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_4.E());
	((TH1D*)gDirectory->FindObject(h_true_M_5->GetName()))->Fill( true_5.M());
	((TH1D*)gDirectory->FindObject(h_true_P_5->GetName()))->Fill( true_5.P());
	((TH1D*)gDirectory->FindObject(h_true_Pt_5->GetName()))->Fill(true_5.Pt());
	((TH1D*)gDirectory->FindObject(h_true_E_5->GetName()))->Fill( true_5.E());
	((TH1D*)gDirectory->FindObject(h_true_M_6->GetName()))->Fill( true_6.M());
	((TH1D*)gDirectory->FindObject(h_true_P_6->GetName()))->Fill( true_6.P());
	((TH1D*)gDirectory->FindObject(h_true_Pt_6->GetName()))->Fill(true_6.Pt());
	((TH1D*)gDirectory->FindObject(h_true_E_6->GetName()))->Fill( true_6.E());
	((TH1D*)gDirectory->FindObject(h_true_M_7->GetName()))->Fill( true_7.M());
	((TH1D*)gDirectory->FindObject(h_true_P_7->GetName()))->Fill( true_7.P());
	((TH1D*)gDirectory->FindObject(h_true_Pt_7->GetName()))->Fill(true_7.Pt());
	((TH1D*)gDirectory->FindObject(h_true_E_7->GetName()))->Fill( true_7.E());
      
      }

      ((TH1D*)gDirectory->FindObject( h_vz_initial->GetName()))->Fill( sevt->vtx[0].z);
      ((TH1I*)gDirectory->FindObject(h_ntrack->GetName()))->Fill(ntrack);
      ((TH1I*)gDirectory->FindObject(h_nclust->GetName()))->Fill(nclust);
  
      ((TH1I*)gDirectory->FindObject(h_evttypebytrks->GetName()))->Fill(evttypebytrks);
      ((TH1I*)gDirectory->FindObject(h_nelectrack->GetName()))->Fill(nelectrack);
      ((TH1I*)gDirectory->FindObject(h_mutrack->GetName()))->Fill(nmutrack);
      ((TH2I*)gDirectory->FindObject(h_neltrk_mutrk->GetName()))->Fill(nelectrack, nmutrack);
      ((TH1I*)gDirectory->FindObject(h_mutrackq->GetName()))->Fill(muchtrkq);
      ((TH1I*)gDirectory->FindObject(h_e1trackq->GetName()))->Fill(eltrkq1);
      ((TH1I*)gDirectory->FindObject(h_e2trackq->GetName()))->Fill(eltrkq2);

      ((TH1F*)gDirectory->FindObject(h_el1trktime->GetName()))        ->Fill( el1trktime );
      ((TH1F*)gDirectory->FindObject(h_el1trkhodtime->GetName()))        ->Fill( el1trkhodtime );
      ((TH1F*)gDirectory->FindObject(h_el1trktimediff->GetName()))    ->Fill( el1trktime - el1trkhodtime );
      ((TH1F*)gDirectory->FindObject(h_ecltime1->GetName()))            ->Fill( el1cltime );
      ((TH1F*)gDirectory->FindObject(h_e1trkcltimediff->GetName()))    ->Fill( el1trktime- el1cltime );

      ((TH1F*)gDirectory->FindObject(h_el2trktime->GetName()))        ->Fill( el2trktime );
      ((TH1F*)gDirectory->FindObject(h_el2trkhodtime->GetName()))        ->Fill( el2trkhodtime );
      ((TH1F*)gDirectory->FindObject(h_el2trktimediff->GetName()))    ->Fill( el2trktime - el2trkhodtime );
      ((TH1F*)gDirectory->FindObject(h_ecltime2->GetName()))            ->Fill( el2cltime );
      ((TH1F*)gDirectory->FindObject(h_e2trkcltimediff->GetName()))    ->Fill( el2trktime-el2cltime );

      ((TH1F*)gDirectory->FindObject(h_e1e2timediff->GetName()))        ->Fill( el1trktime-el2trktime );
      ((TH1F*)gDirectory->FindObject(h_e1e2cltimediff->GetName()))    ->Fill( el1cltime-el2cltime );
      ((TH1F*)gDirectory->FindObject(h_e1e2hodtimediff->GetName()))        ->Fill( el1trkhodtime-el2trkhodtime );
      ((TH1F*)gDirectory->FindObject(h_e1muhodtimediff->GetName()))        ->Fill( el1trkhodtime-muchtrkhodtime );
      ((TH1F*)gDirectory->FindObject(h_e2muhodtimediff->GetName()))        ->Fill( el2trkhodtime-muchtrkhodtime );

      ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_mu_track);
      ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_e1_track);
      ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_e2_track);

  
      ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_mu_track, (*muchlv).P());
      ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_e1_track, (*el1lv).P() );
      ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_e2_track, (*el2lv).P() );
  
      ((TH1F*)gDirectory->FindObject(h_eop_mu->GetName()))->Fill(eop_mu_track);
      ((TH2F*)gDirectory->FindObject(h_eop_mu_vs_p->GetName()))->Fill(eop_mu_track, (*muchlv).P());
      ((TH1F*)gDirectory->FindObject(h_eop_el->GetName()))->Fill(eop_e1_track);
      ((TH2F*)gDirectory->FindObject(h_eop_el_vs_p->GetName()))->Fill(eop_e1_track, (*el1lv).P());
      ((TH1F*)gDirectory->FindObject(h_eop_el->GetName()))->Fill(eop_e2_track);
      ((TH2F*)gDirectory->FindObject(h_eop_el_vs_p->GetName()))->Fill(eop_e2_track, (*el2lv).P());

      if (eop_e2_track <= eop_e1_track ){
	((TH2F*)gDirectory->FindObject(h_eop_vs_Pt->GetName()))->Fill(eop_e2_track,(*muvee).Pt() );
      } else {
	((TH2F*)gDirectory->FindObject(h_eop_vs_Pt->GetName()))->Fill(eop_e1_track, (*muvee).Pt());
      }

    

      ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill( elp1_uncorr[0] );
      ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill( elp1_uncorr[1] );
      ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill( elv1_uncorr[0] );
      ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill( elv1_uncorr[1] );
      ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill( elp1_uncorr[0], elp1_uncorr[1] );
      ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill( elv1_uncorr[0], elv1_uncorr[1] );

      ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( el1lv->P() );
      ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( el2lv->P() );
      ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el1lv->E() );
      ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el2lv->E() );
      ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el1lv->P() );
      ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el2lv->P() );
      ((TH1F*)gDirectory->FindObject(h_e1m->GetName()))->Fill( el1lv->M() );
      ((TH1F*)gDirectory->FindObject(h_e2m->GetName()))->Fill( el2lv->M() );
      ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el1lv->Pt() );
      ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el2lv->Pt() );
      ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
      ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
      ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );
      ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);
      ((TH1F*)gDirectory->FindObject(h_angle_e1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()) );
      ((TH1F*)gDirectory->FindObject(h_angle_e1e2_mu->GetName()))->Fill( (*eelv).Vect().Angle(muchlv->Vect()) );

      ((TH1F*)gDirectory->FindObject(h_muchtrktime->GetName()))        ->Fill(muchtrktime);
      ((TH1F*)gDirectory->FindObject(h_muchtrkhodtime->GetName()))        ->Fill( muchtrkhodtime );
      ((TH1F*)gDirectory->FindObject(h_muchtrktimediff->GetName()))    ->Fill( muchtrktime-muchtrkhodtime );
      //((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( muchtrktime-muchcltime );

      ((TH1F*)gDirectory->FindObject(h_muchlvE->GetName()))->Fill( (*muchlv).E() );
      ((TH1F*)gDirectory->FindObject(h_muchlvP->GetName()))->Fill( (*muchlv).P() );
      ((TH1F*)gDirectory->FindObject(h_muchlvPt->GetName()))->Fill((*muchlv).Pt() );

      ((TH1F*)gDirectory->FindObject(h_KE->GetName()))->Fill((*pK).E());
      ((TH1F*)gDirectory->FindObject(h_KP->GetName()))->Fill((*pK).P());

      ((TH1F*)gDirectory->FindObject(h_nu_P->GetName()))->Fill(  (*neutrino).P() );
      ((TH1F*)gDirectory->FindObject(h_nu_E->GetName()))->Fill(  (*neutrino).E() );
      ((TH1F*)gDirectory->FindObject(h_nu_Pt->GetName()))->Fill( (*neutrino).Pt() );
      ((TH1F*)gDirectory->FindObject(h_nu_m->GetName()))->Fill(  (*neutrino).M() );
      ((TH1F*)gDirectory->FindObject(h_missing_mass->GetName()))->Fill( (*neutrino).M2() );

      if( muchtrkq+eltrkq1+eltrkq2 == 1){
	((TH1F*)gDirectory->FindObject(h_missing_mass_kp->GetName()))->Fill( (*neutrino).M2() );
      } else {
	((TH1F*)gDirectory->FindObject(h_missing_mass_km->GetName()))->Fill( (*neutrino).M2() );
      }

      ((TH1F*)gDirectory->FindObject(h_mueeE->GetName()))->Fill( (*muvee).E());
      ((TH1F*)gDirectory->FindObject(h_mueeP->GetName()))->Fill( (*muvee).P());
      ((TH1F*)gDirectory->FindObject(h_mueePt->GetName()))->Fill((*muvee).Pt());
      ((TH1F*)gDirectory->FindObject(h_mueeM->GetName()))->Fill ((*muvee).M());
      ((TH1F*)gDirectory->FindObject(h_3trk_invm->GetName()))->Fill ((*muvee).M2());


      ((TH1F*)gDirectory->FindObject(h_bx_much->GetName()))            ->Fill(muchp_uncorr[0]);
      ((TH1F*)gDirectory->FindObject(h_by_much->GetName()))            ->Fill(muchp_uncorr[1]);
      ((TH1F*)gDirectory->FindObject(h_bdxdz_much->GetName()))            ->Fill(muchv_uncorr[0]);
      ((TH1F*)gDirectory->FindObject(h_bdydz_much->GetName()))            ->Fill(muchv_uncorr[1]) ;
      ((TH1F*)gDirectory->FindObject(h_bx_vs_by_much->GetName()))        ->Fill(muchp_uncorr[0], muchp_uncorr[1]);
      ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_much->GetName()))->Fill(muchv_uncorr[0], muchv_uncorr[1]);

      ((TH1F*)gDirectory->FindObject(h_cda_ee->GetName()))->Fill(cda);
      ((TH1F*)gDirectory->FindObject(h_vx->GetName()))->Fill(vertex[0]);
      ((TH1F*)gDirectory->FindObject(h_vy->GetName()))->Fill(vertex[1]);
      ((TH1F*)gDirectory->FindObject(h_vz->GetName()))->Fill(vertex[2]);
      ((TH2F*)gDirectory->FindObject(h_vx_vs_vy->GetName()))->Fill(vertex[0],vertex[1]);

      ((TH1F*)gDirectory->FindObject(h_vx_e1much->GetName()))->Fill(vertex_e1much[0]);
      ((TH1F*)gDirectory->FindObject(h_vy_e1much->GetName()))->Fill(vertex_e1much[1]);
      ((TH1F*)gDirectory->FindObject(h_vz_e1much->GetName()))->Fill(vertex_e1much[2]);
      ((TH1F*)gDirectory->FindObject(h_cda_e1much->GetName()))->Fill(cda_e1much);
      ((TH1F*)gDirectory->FindObject(h_vz_e1much_diff->GetName()))->Fill( vertex[2] - vertex_e1much[2] );

      ((TH1F*)gDirectory->FindObject(h_vx_e2much->GetName()))->Fill(vertex_e2much[0]);
      ((TH1F*)gDirectory->FindObject(h_vy_e2much->GetName()))->Fill(vertex_e2much[1]);
      ((TH1F*)gDirectory->FindObject(h_vz_e2much->GetName()))->Fill(vertex_e2much[2]);
      ((TH1F*)gDirectory->FindObject(h_cda_e2much->GetName()))->Fill(cda_e2much);
      ((TH1F*)gDirectory->FindObject(h_vz_e2much_diff->GetName()))->Fill( vertex[2] - vertex_e2much[2] );
      ((TH1F*)gDirectory->FindObject(h_vz_mu1mu2_diff->GetName()))->Fill( vertex_e1much[2] - vertex_e2much[2] );

      
      ((TH1F*)gDirectory->FindObject(h_cda_e1much->GetName()))->Fill(cda_e1much);
      ((TH2F*)gDirectory->FindObject(h_vx_vs_vy_e1much->GetName()))->Fill(vertex_e1much[0],vertex_e1much[1]);
    }
  }

  file1->cd();
  h_cutflow->Fill(cutcounter);
  h_cutflow->GetXaxis()->SetBinLabel( cutcounter, cuts[cutcounter-1].c_str() );
  cutcounter++;
  ///THE END OF ************* K3pi selection ****************** ///

  /// ************* Muee transverse momentum cut < 0.02 ****************** ///
  //"muee energy > 70GeV test" e-e+ transverse momentum > 0.9 
  //if( (*muvee).E() < 70 )        {cleanup();return 0;}
  //if( (*eelv).Pt() > 0.9 )        {
  //if( sevt->cluster[sevt->track[eltrkindx1].iClus].energy/sevt->track[eltrkindx1].p < 0.95 && sevt->cluster[sevt->track[eltrkindx2].iClus].energy/sevt->track[eltrkindx2].p < 0.95 )        {
  //if (flux_id==1){
  //  if ((*muvee).Pt()*(*muvee).Pt() > 0.0005 ) {cleanup();return 0;}
  //}
  //if (flux_id==0){
  if ((*muvee).Pt()*(*muvee).Pt() < 0.0005 ) {cleanup();return 0;}
  //}
  dirs[cutcounter-1]->cd();

  if (IS_MC){
    //if(true_4.P()  < 3. ||true_5.P() < 3.) {cleanup();return 0;}
    ((TH1D*)gDirectory->FindObject( h_mc_muvee_M->GetName()))->Fill( Muee_three_track.M());
    ((TH1D*)gDirectory->FindObject( h_mc_muvee_P->GetName()))->Fill( Muee_three_track.P());
    ((TH1D*)gDirectory->FindObject( h_mc_muvee_Pt->GetName()))->Fill( Muee_three_track.Pt());
    ((TH1D*)gDirectory->FindObject( h_mc_elsP->GetName()))->Fill( mc_ee.P());
    ((TH1D*)gDirectory->FindObject( h_mc_elsM->GetName()))->Fill( mc_ee.M());
    ((TH1D*)gDirectory->FindObject( h_mc_nu_M->GetName()))->Fill( Nu_true.M());
    ((TH1D*)gDirectory->FindObject( h_mc_nu_M2->GetName()))->Fill( Nu_true.M2());
    ((TH1D*)gDirectory->FindObject( h_mc_nu_P->GetName()))->Fill( Nu_true.P());
    ((TH1D*)gDirectory->FindObject( h_mc_nu_Pt->GetName()))->Fill(Nu_true.Pt());
    
    ((TH1D*)gDirectory->FindObject(h_true_M_0->GetName()))->Fill( true_0.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_0->GetName()))->Fill( true_0.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_0->GetName()))->Fill(true_0.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_0->GetName()))->Fill( true_0.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_1->GetName()))->Fill( true_1.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_1->GetName()))->Fill( true_1.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_1->GetName()))->Fill(true_1.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_1->GetName()))->Fill( true_1.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_2->GetName()))->Fill( true_2.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_2->GetName()))->Fill( true_2.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_2->GetName()))->Fill(true_2.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_2->GetName()))->Fill( true_2.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_3->GetName()))->Fill( true_3.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_3->GetName()))->Fill( true_3.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_3->GetName()))->Fill(true_3.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_3->GetName()))->Fill( true_3.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_4.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_4.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_4.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_4.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_5->GetName()))->Fill( true_5.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_5->GetName()))->Fill( true_5.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_5->GetName()))->Fill(true_5.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_5->GetName()))->Fill( true_5.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_6->GetName()))->Fill( true_6.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_6->GetName()))->Fill( true_6.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_6->GetName()))->Fill(true_6.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_6->GetName()))->Fill( true_6.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_7->GetName()))->Fill( true_7.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_7->GetName()))->Fill( true_7.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_7->GetName()))->Fill(true_7.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_7->GetName()))->Fill( true_7.E());
      
  }

  ((TH1D*)gDirectory->FindObject( h_vz_initial->GetName()))->Fill( sevt->vtx[0].z);
    
  ((TH1I*)gDirectory->FindObject(h_ntrack->GetName()))->Fill(ntrack);
  ((TH1I*)gDirectory->FindObject(h_nclust->GetName()))->Fill(nclust);
  
  ((TH1I*)gDirectory->FindObject(h_evttypebytrks->GetName()))->Fill(evttypebytrks);
  ((TH1I*)gDirectory->FindObject(h_nelectrack->GetName()))->Fill(nelectrack);
  ((TH1I*)gDirectory->FindObject(h_mutrack->GetName()))->Fill(nmutrack);
  ((TH2I*)gDirectory->FindObject(h_neltrk_mutrk->GetName()))->Fill(nelectrack, nmutrack);
  ((TH1I*)gDirectory->FindObject(h_mutrackq->GetName()))->Fill(muchtrkq);
  ((TH1I*)gDirectory->FindObject(h_e1trackq->GetName()))->Fill(eltrkq1);
  ((TH1I*)gDirectory->FindObject(h_e2trackq->GetName()))->Fill(eltrkq2);
    
  ((TH1F*)gDirectory->FindObject(h_el1trktime->GetName()))        ->Fill( el1trktime );
  ((TH1F*)gDirectory->FindObject(h_el1trkhodtime->GetName()))        ->Fill( el1trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_el1trktimediff->GetName()))    ->Fill( el1trktime - el1trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_ecltime1->GetName()))            ->Fill( el1cltime );
  ((TH1F*)gDirectory->FindObject(h_e1trkcltimediff->GetName()))    ->Fill( el1trktime- el1cltime );

  ((TH1F*)gDirectory->FindObject(h_el2trktime->GetName()))        ->Fill( el2trktime );
  ((TH1F*)gDirectory->FindObject(h_el2trkhodtime->GetName()))        ->Fill( el2trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_el2trktimediff->GetName()))    ->Fill( el2trktime - el2trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_ecltime2->GetName()))            ->Fill( el2cltime );
  ((TH1F*)gDirectory->FindObject(h_e2trkcltimediff->GetName()))    ->Fill( el2trktime-el2cltime );

  ((TH1F*)gDirectory->FindObject(h_e1e2timediff->GetName()))        ->Fill( el1trktime-el2trktime );
  ((TH1F*)gDirectory->FindObject(h_e1e2cltimediff->GetName()))    ->Fill( el1cltime-el2cltime );
  ((TH1F*)gDirectory->FindObject(h_e1e2hodtimediff->GetName()))        ->Fill( el1trkhodtime-el2trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_e1muhodtimediff->GetName()))        ->Fill( el1trkhodtime-muchtrkhodtime );
  ((TH1F*)gDirectory->FindObject(h_e2muhodtimediff->GetName()))        ->Fill( el2trkhodtime-muchtrkhodtime );

  ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_mu_track);
  ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_e1_track);
  ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_e2_track);

  
  ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_mu_track, (*muchlv).P());
  ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_e1_track, (*el1lv).P() );
  ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_e2_track, (*el2lv).P() );
  
  ((TH1F*)gDirectory->FindObject(h_eop_mu->GetName()))->Fill(eop_mu_track);
  ((TH2F*)gDirectory->FindObject(h_eop_mu_vs_p->GetName()))->Fill(eop_mu_track, (*muchlv).P());
  ((TH1F*)gDirectory->FindObject(h_eop_el->GetName()))->Fill(eop_e1_track);
  ((TH2F*)gDirectory->FindObject(h_eop_el_vs_p->GetName()))->Fill(eop_e1_track, (*el1lv).P());
  ((TH1F*)gDirectory->FindObject(h_eop_el->GetName()))->Fill(eop_e2_track);
  ((TH2F*)gDirectory->FindObject(h_eop_el_vs_p->GetName()))->Fill(eop_e2_track, (*el2lv).P());

       
  ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill( elp1_uncorr[0] );
  ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill( elp1_uncorr[1] );
  ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill( elv1_uncorr[0] );
  ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill( elv1_uncorr[1] );
  ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill( elp1_uncorr[0], elp1_uncorr[1] );
  ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill( elv1_uncorr[0], elv1_uncorr[1] );

  ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( el1lv->P() );
  ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( el2lv->P() );
  ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el1lv->E() );
  ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el2lv->E() );
  ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el1lv->P() );
  ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el2lv->P() );
  ((TH1F*)gDirectory->FindObject(h_e1m->GetName()))->Fill( el1lv->M() );
  ((TH1F*)gDirectory->FindObject(h_e2m->GetName()))->Fill( el2lv->M() );
  ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el1lv->Pt() );
  ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el2lv->Pt() );
  ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
  ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
  ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );
  ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);
  ((TH1F*)gDirectory->FindObject(h_angle_e1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()) );
  ((TH1F*)gDirectory->FindObject(h_angle_e1e2_mu->GetName()))->Fill( (*eelv).Vect().Angle(muchlv->Vect()) );

  ((TH1F*)gDirectory->FindObject(h_muchtrktime->GetName()))        ->Fill(muchtrktime);
  ((TH1F*)gDirectory->FindObject(h_muchtrkhodtime->GetName()))        ->Fill( muchtrkhodtime );
  ((TH1F*)gDirectory->FindObject(h_muchtrktimediff->GetName()))    ->Fill( muchtrktime-muchtrkhodtime );
  //((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( muchtrktime-muchcltime );

  ((TH1F*)gDirectory->FindObject(h_muchlvE->GetName()))->Fill( (*muchlv).E() );
  ((TH1F*)gDirectory->FindObject(h_muchlvP->GetName()))->Fill( (*muchlv).P() );
  ((TH1F*)gDirectory->FindObject(h_muchlvPt->GetName()))->Fill((*muchlv).Pt() );

  ((TH1F*)gDirectory->FindObject(h_KE->GetName()))->Fill((*pK).E());
  ((TH1F*)gDirectory->FindObject(h_KP->GetName()))->Fill((*pK).P());
    
  ((TH1F*)gDirectory->FindObject(h_nu_P->GetName()))->Fill(  (*neutrino).P() );
  ((TH1F*)gDirectory->FindObject(h_nu_E->GetName()))->Fill(  (*neutrino).E() );
  ((TH1F*)gDirectory->FindObject(h_nu_Pt->GetName()))->Fill( (*neutrino).Pt() );
  ((TH1F*)gDirectory->FindObject(h_nu_m->GetName()))->Fill(  (*neutrino).M() );
  ((TH1F*)gDirectory->FindObject(h_missing_mass->GetName()))->Fill( (*neutrino).M2() );

  if( muchtrkq+eltrkq1+eltrkq2 == 1){
    ((TH1F*)gDirectory->FindObject(h_missing_mass_kp->GetName()))->Fill( (*neutrino).M2() );
  } else {
    ((TH1F*)gDirectory->FindObject(h_missing_mass_km->GetName()))->Fill( (*neutrino).M2() );
  }

  ((TH1F*)gDirectory->FindObject(h_mueeE->GetName()))->Fill( (*muvee).E());
  ((TH1F*)gDirectory->FindObject(h_mueeP->GetName()))->Fill( (*muvee).P());
  ((TH1F*)gDirectory->FindObject(h_mueePt->GetName()))->Fill((*muvee).Pt());
  ((TH1F*)gDirectory->FindObject(h_mueeM->GetName()))->Fill ((*muvee).M());
  ((TH1F*)gDirectory->FindObject(h_3trk_invm->GetName()))->Fill ((*muvee).M2());
    
    
  ((TH1F*)gDirectory->FindObject(h_bx_much->GetName()))            ->Fill(muchp_uncorr[0]);
  ((TH1F*)gDirectory->FindObject(h_by_much->GetName()))            ->Fill(muchp_uncorr[1]);
  ((TH1F*)gDirectory->FindObject(h_bdxdz_much->GetName()))            ->Fill(muchv_uncorr[0]);
  ((TH1F*)gDirectory->FindObject(h_bdydz_much->GetName()))            ->Fill(muchv_uncorr[1]) ;
  ((TH1F*)gDirectory->FindObject(h_bx_vs_by_much->GetName()))        ->Fill(muchp_uncorr[0], muchp_uncorr[1]);
  ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_much->GetName()))->Fill(muchv_uncorr[0], muchv_uncorr[1]);

  ((TH1F*)gDirectory->FindObject(h_cda_ee->GetName()))->Fill(cda);
  ((TH1F*)gDirectory->FindObject(h_vx->GetName()))->Fill(vertex[0]);
  ((TH1F*)gDirectory->FindObject(h_vy->GetName()))->Fill(vertex[1]);
  ((TH1F*)gDirectory->FindObject(h_vz->GetName()))->Fill(vertex[2]);
  ((TH2F*)gDirectory->FindObject(h_vx_vs_vy->GetName()))->Fill(vertex[0],vertex[1]);

  ((TH1F*)gDirectory->FindObject(h_vx_e1much->GetName()))->Fill(vertex_e1much[0]);
  ((TH1F*)gDirectory->FindObject(h_vy_e1much->GetName()))->Fill(vertex_e1much[1]);
  ((TH1F*)gDirectory->FindObject(h_vz_e1much->GetName()))->Fill(vertex_e1much[2]);
  ((TH1F*)gDirectory->FindObject(h_cda_e1much->GetName()))->Fill(cda_e1much);
  ((TH1F*)gDirectory->FindObject(h_vz_e1much_diff->GetName()))->Fill( vertex[2] - vertex_e1much[2] );

  ((TH1F*)gDirectory->FindObject(h_vx_e2much->GetName()))->Fill(vertex_e2much[0]);
  ((TH1F*)gDirectory->FindObject(h_vy_e2much->GetName()))->Fill(vertex_e2much[1]);
  ((TH1F*)gDirectory->FindObject(h_vz_e2much->GetName()))->Fill(vertex_e2much[2]);
  ((TH1F*)gDirectory->FindObject(h_cda_e2much->GetName()))->Fill(cda_e2much);
  ((TH1F*)gDirectory->FindObject(h_vz_e2much_diff->GetName()))->Fill( vertex[2] - vertex_e2much[2] );
  ((TH1F*)gDirectory->FindObject(h_vz_mu1mu2_diff->GetName()))->Fill( vertex_e1much[2] - vertex_e2much[2] );

  

  ((TH1F*)gDirectory->FindObject(h_cda_e1much->GetName()))->Fill(cda_e1much);
  ((TH2F*)gDirectory->FindObject(h_vx_vs_vy_e1much->GetName()))->Fill(vertex_e1much[0],vertex_e1much[1]);

    
  

  file1->cd();
  h_cutflow->Fill(cutcounter);
  h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str()); 
  cutcounter++;
  ///THE END OF *************electron E/p between 0.9 and 0.95 (both)****************** ///
  // e+e- transverse momentum > 0.9 

  /// ************* Final k3pi selection ****************** ///
  dirs[cutcounter-1]->cd();
  if (nelectrack==2 && muchtrkindx!=-1) {
    if (eltrkq1==eltrkq2    &&
	eltrkq1!=muchtrkq   &&
	muchtrkq==1         &&
	eelvmass >= 0.140 
	)    {
      if(eop_e2_track <=0.95 || eop_e1_track <=0.95 ) {cleanup();return 0;}
      if (IS_MC){
	//if(true_4.P()  < 3. ||true_5.P() < 3.) {cleanup();return 0;}
	((TH1D*)gDirectory->FindObject( h_mc_muvee_M->GetName()))->Fill( Muee_three_track.M());
	((TH1D*)gDirectory->FindObject( h_mc_muvee_P->GetName()))->Fill( Muee_three_track.P());
	((TH1D*)gDirectory->FindObject( h_mc_muvee_Pt->GetName()))->Fill( Muee_three_track.Pt());
	((TH1D*)gDirectory->FindObject( h_mc_elsP->GetName()))->Fill( mc_ee.P());
	((TH1D*)gDirectory->FindObject( h_mc_elsM->GetName()))->Fill( mc_ee.M());
	((TH1D*)gDirectory->FindObject( h_mc_nu_M->GetName()))->Fill( Nu_true.M());
	((TH1D*)gDirectory->FindObject( h_mc_nu_M2->GetName()))->Fill( Nu_true.M2());
	((TH1D*)gDirectory->FindObject( h_mc_nu_P->GetName()))->Fill( Nu_true.P());
	((TH1D*)gDirectory->FindObject( h_mc_nu_Pt->GetName()))->Fill(Nu_true.Pt());
  
	((TH1D*)gDirectory->FindObject(h_true_M_0->GetName()))->Fill( true_0.M());
	((TH1D*)gDirectory->FindObject(h_true_P_0->GetName()))->Fill( true_0.P());
	((TH1D*)gDirectory->FindObject(h_true_Pt_0->GetName()))->Fill(true_0.Pt());
	((TH1D*)gDirectory->FindObject(h_true_E_0->GetName()))->Fill( true_0.E());
	((TH1D*)gDirectory->FindObject(h_true_M_1->GetName()))->Fill( true_1.M());
	((TH1D*)gDirectory->FindObject(h_true_P_1->GetName()))->Fill( true_1.P());
	((TH1D*)gDirectory->FindObject(h_true_Pt_1->GetName()))->Fill(true_1.Pt());
	((TH1D*)gDirectory->FindObject(h_true_E_1->GetName()))->Fill( true_1.E());
	((TH1D*)gDirectory->FindObject(h_true_M_2->GetName()))->Fill( true_2.M());
	((TH1D*)gDirectory->FindObject(h_true_P_2->GetName()))->Fill( true_2.P());
	((TH1D*)gDirectory->FindObject(h_true_Pt_2->GetName()))->Fill(true_2.Pt());
	((TH1D*)gDirectory->FindObject(h_true_E_2->GetName()))->Fill( true_2.E());
	((TH1D*)gDirectory->FindObject(h_true_M_3->GetName()))->Fill( true_3.M());
	((TH1D*)gDirectory->FindObject(h_true_P_3->GetName()))->Fill( true_3.P());
	((TH1D*)gDirectory->FindObject(h_true_Pt_3->GetName()))->Fill(true_3.Pt());
	((TH1D*)gDirectory->FindObject(h_true_E_3->GetName()))->Fill( true_3.E());
	((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_4.M());
	((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_4.P());
	((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_4.Pt());
	((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_4.E());
	((TH1D*)gDirectory->FindObject(h_true_M_5->GetName()))->Fill( true_5.M());
	((TH1D*)gDirectory->FindObject(h_true_P_5->GetName()))->Fill( true_5.P());
	((TH1D*)gDirectory->FindObject(h_true_Pt_5->GetName()))->Fill(true_5.Pt());
	((TH1D*)gDirectory->FindObject(h_true_E_5->GetName()))->Fill( true_5.E());
	((TH1D*)gDirectory->FindObject(h_true_M_6->GetName()))->Fill( true_6.M());
	((TH1D*)gDirectory->FindObject(h_true_P_6->GetName()))->Fill( true_6.P());
	((TH1D*)gDirectory->FindObject(h_true_Pt_6->GetName()))->Fill(true_6.Pt());
	((TH1D*)gDirectory->FindObject(h_true_E_6->GetName()))->Fill( true_6.E());
	((TH1D*)gDirectory->FindObject(h_true_M_7->GetName()))->Fill( true_7.M());
	((TH1D*)gDirectory->FindObject(h_true_P_7->GetName()))->Fill( true_7.P());
	((TH1D*)gDirectory->FindObject(h_true_Pt_7->GetName()))->Fill(true_7.Pt());
	((TH1D*)gDirectory->FindObject(h_true_E_7->GetName()))->Fill( true_7.E());
      
      }

      ((TH1D*)gDirectory->FindObject( h_vz_initial->GetName()))->Fill( sevt->vtx[0].z);
      ((TH1I*)gDirectory->FindObject(h_ntrack->GetName()))->Fill(ntrack);
      ((TH1I*)gDirectory->FindObject(h_nclust->GetName()))->Fill(nclust);
  
      ((TH1I*)gDirectory->FindObject(h_evttypebytrks->GetName()))->Fill(evttypebytrks);
      ((TH1I*)gDirectory->FindObject(h_nelectrack->GetName()))->Fill(nelectrack);
      ((TH1I*)gDirectory->FindObject(h_mutrack->GetName()))->Fill(nmutrack);
      ((TH2I*)gDirectory->FindObject(h_neltrk_mutrk->GetName()))->Fill(nelectrack, nmutrack);
      ((TH1I*)gDirectory->FindObject(h_mutrackq->GetName()))->Fill(muchtrkq);
      ((TH1I*)gDirectory->FindObject(h_e1trackq->GetName()))->Fill(eltrkq1);
      ((TH1I*)gDirectory->FindObject(h_e2trackq->GetName()))->Fill(eltrkq2);

      ((TH1F*)gDirectory->FindObject(h_el1trktime->GetName()))        ->Fill( el1trktime );
      ((TH1F*)gDirectory->FindObject(h_el1trkhodtime->GetName()))        ->Fill( el1trkhodtime );
      ((TH1F*)gDirectory->FindObject(h_el1trktimediff->GetName()))    ->Fill( el1trktime - el1trkhodtime );
      ((TH1F*)gDirectory->FindObject(h_ecltime1->GetName()))            ->Fill( el1cltime );
      ((TH1F*)gDirectory->FindObject(h_e1trkcltimediff->GetName()))    ->Fill( el1trktime- el1cltime );

      ((TH1F*)gDirectory->FindObject(h_el2trktime->GetName()))        ->Fill( el2trktime );
      ((TH1F*)gDirectory->FindObject(h_el2trkhodtime->GetName()))        ->Fill( el2trkhodtime );
      ((TH1F*)gDirectory->FindObject(h_el2trktimediff->GetName()))    ->Fill( el2trktime - el2trkhodtime );
      ((TH1F*)gDirectory->FindObject(h_ecltime2->GetName()))            ->Fill( el2cltime );
      ((TH1F*)gDirectory->FindObject(h_e2trkcltimediff->GetName()))    ->Fill( el2trktime-el2cltime );

      ((TH1F*)gDirectory->FindObject(h_e1e2timediff->GetName()))        ->Fill( el1trktime-el2trktime );
      ((TH1F*)gDirectory->FindObject(h_e1e2cltimediff->GetName()))    ->Fill( el1cltime-el2cltime );
      ((TH1F*)gDirectory->FindObject(h_e1e2hodtimediff->GetName()))        ->Fill( el1trkhodtime-el2trkhodtime );
      ((TH1F*)gDirectory->FindObject(h_e1muhodtimediff->GetName()))        ->Fill( el1trkhodtime-muchtrkhodtime );
      ((TH1F*)gDirectory->FindObject(h_e2muhodtimediff->GetName()))        ->Fill( el2trkhodtime-muchtrkhodtime );
      ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_mu_track);
      ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_e1_track);
      ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_e2_track);

  
      ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_mu_track, (*muchlv).P());
      ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_e1_track, (*el1lv).P() );
      ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_e2_track, (*el2lv).P() );
  
      ((TH1F*)gDirectory->FindObject(h_eop_mu->GetName()))->Fill(eop_mu_track);
      ((TH2F*)gDirectory->FindObject(h_eop_mu_vs_p->GetName()))->Fill(eop_mu_track, (*muchlv).P());
      ((TH1F*)gDirectory->FindObject(h_eop_el->GetName()))->Fill(eop_e1_track);
      ((TH2F*)gDirectory->FindObject(h_eop_el_vs_p->GetName()))->Fill(eop_e1_track, (*el1lv).P());
      ((TH1F*)gDirectory->FindObject(h_eop_el->GetName()))->Fill(eop_e2_track);
      ((TH2F*)gDirectory->FindObject(h_eop_el_vs_p->GetName()))->Fill(eop_e2_track, (*el2lv).P());

   
      if (eop_e2_track <= eop_e1_track ){
	((TH2F*)gDirectory->FindObject(h_eop_vs_Pt->GetName()))->Fill(eop_e2_track,(*muvee).Pt() );
      } else {
	((TH2F*)gDirectory->FindObject(h_eop_vs_Pt->GetName()))->Fill(eop_e1_track, (*muvee).Pt());
      }
      

      ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill( elp1_uncorr[0] );
      ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill( elp1_uncorr[1] );
      ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill( elv1_uncorr[0] );
      ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill( elv1_uncorr[1] );
      ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill( elp1_uncorr[0], elp1_uncorr[1] );
      ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill( elv1_uncorr[0], elv1_uncorr[1] );

      ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( el1lv->P() );
      ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( el2lv->P() );
      ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el1lv->E() );
      ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el2lv->E() );
      ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el1lv->P() );
      ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el2lv->P() );
      ((TH1F*)gDirectory->FindObject(h_e1m->GetName()))->Fill( el1lv->M() );
      ((TH1F*)gDirectory->FindObject(h_e2m->GetName()))->Fill( el2lv->M() );
      ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el1lv->Pt() );
      ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el2lv->Pt() );
      ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
      ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
      ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );
      ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);
      ((TH1F*)gDirectory->FindObject(h_angle_e1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()) );
      ((TH1F*)gDirectory->FindObject(h_angle_e1e2_mu->GetName()))->Fill( (*eelv).Vect().Angle(muchlv->Vect()) );

      ((TH1F*)gDirectory->FindObject(h_muchtrktime->GetName()))        ->Fill(muchtrktime);
      ((TH1F*)gDirectory->FindObject(h_muchtrkhodtime->GetName()))        ->Fill( muchtrkhodtime );
      ((TH1F*)gDirectory->FindObject(h_muchtrktimediff->GetName()))    ->Fill( muchtrktime-muchtrkhodtime );
      //((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( muchtrktime-muchcltime );

      ((TH1F*)gDirectory->FindObject(h_muchlvE->GetName()))->Fill( (*muchlv).E() );
      ((TH1F*)gDirectory->FindObject(h_muchlvP->GetName()))->Fill( (*muchlv).P() );
      ((TH1F*)gDirectory->FindObject(h_muchlvPt->GetName()))->Fill((*muchlv).Pt() );

      ((TH1F*)gDirectory->FindObject(h_KE->GetName()))->Fill((*pK).E());
      ((TH1F*)gDirectory->FindObject(h_KP->GetName()))->Fill((*pK).P());

      ((TH1F*)gDirectory->FindObject(h_nu_P->GetName()))->Fill(  (*neutrino).P() );
      ((TH1F*)gDirectory->FindObject(h_nu_E->GetName()))->Fill(  (*neutrino).E() );
      ((TH1F*)gDirectory->FindObject(h_nu_Pt->GetName()))->Fill( (*neutrino).Pt() );
      ((TH1F*)gDirectory->FindObject(h_nu_m->GetName()))->Fill(  (*neutrino).M() );
      ((TH1F*)gDirectory->FindObject(h_missing_mass->GetName()))->Fill( (*neutrino).M2() );

      if( muchtrkq+eltrkq1+eltrkq2 == 1){
	((TH1F*)gDirectory->FindObject(h_missing_mass_kp->GetName()))->Fill( (*neutrino).M2() );
      } else {
	((TH1F*)gDirectory->FindObject(h_missing_mass_km->GetName()))->Fill( (*neutrino).M2() );
      }

      ((TH1F*)gDirectory->FindObject(h_mueeE->GetName()))->Fill( (*muvee).E());
      ((TH1F*)gDirectory->FindObject(h_mueeP->GetName()))->Fill( (*muvee).P());
      ((TH1F*)gDirectory->FindObject(h_mueePt->GetName()))->Fill((*muvee).Pt());
      ((TH1F*)gDirectory->FindObject(h_mueeM->GetName()))->Fill ((*muvee).M());
      ((TH1F*)gDirectory->FindObject(h_3trk_invm->GetName()))->Fill ((*muvee).M2());


      ((TH1F*)gDirectory->FindObject(h_bx_much->GetName()))            ->Fill(muchp_uncorr[0]);
      ((TH1F*)gDirectory->FindObject(h_by_much->GetName()))            ->Fill(muchp_uncorr[1]);
      ((TH1F*)gDirectory->FindObject(h_bdxdz_much->GetName()))            ->Fill(muchv_uncorr[0]);
      ((TH1F*)gDirectory->FindObject(h_bdydz_much->GetName()))            ->Fill(muchv_uncorr[1]) ;
      ((TH1F*)gDirectory->FindObject(h_bx_vs_by_much->GetName()))        ->Fill(muchp_uncorr[0], muchp_uncorr[1]);
      ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_much->GetName()))->Fill(muchv_uncorr[0], muchv_uncorr[1]);

      ((TH1F*)gDirectory->FindObject(h_cda_ee->GetName()))->Fill(cda);
      ((TH1F*)gDirectory->FindObject(h_vx->GetName()))->Fill(vertex[0]);
      ((TH1F*)gDirectory->FindObject(h_vy->GetName()))->Fill(vertex[1]);
      ((TH1F*)gDirectory->FindObject(h_vz->GetName()))->Fill(vertex[2]);
      ((TH2F*)gDirectory->FindObject(h_vx_vs_vy->GetName()))->Fill(vertex[0],vertex[1]);

      ((TH1F*)gDirectory->FindObject(h_vx_e1much->GetName()))->Fill(vertex_e1much[0]);
      ((TH1F*)gDirectory->FindObject(h_vy_e1much->GetName()))->Fill(vertex_e1much[1]);
      ((TH1F*)gDirectory->FindObject(h_vz_e1much->GetName()))->Fill(vertex_e1much[2]);
      ((TH1F*)gDirectory->FindObject(h_cda_e1much->GetName()))->Fill(cda_e1much);
      ((TH1F*)gDirectory->FindObject(h_vz_e1much_diff->GetName()))->Fill( vertex[2] - vertex_e1much[2] );

      ((TH1F*)gDirectory->FindObject(h_vx_e2much->GetName()))->Fill(vertex_e2much[0]);
      ((TH1F*)gDirectory->FindObject(h_vy_e2much->GetName()))->Fill(vertex_e2much[1]);
      ((TH1F*)gDirectory->FindObject(h_vz_e2much->GetName()))->Fill(vertex_e2much[2]);
      ((TH1F*)gDirectory->FindObject(h_cda_e2much->GetName()))->Fill(cda_e2much);
      ((TH1F*)gDirectory->FindObject(h_vz_e2much_diff->GetName()))->Fill( vertex[2] - vertex_e2much[2] );
      ((TH1F*)gDirectory->FindObject(h_vz_mu1mu2_diff->GetName()))->Fill( vertex_e1much[2] - vertex_e2much[2] );

      
      ((TH1F*)gDirectory->FindObject(h_cda_e1much->GetName()))->Fill(cda_e1much);
      ((TH2F*)gDirectory->FindObject(h_vx_vs_vy_e1much->GetName()))->Fill(vertex_e1much[0],vertex_e1much[1]);
    }
  }

  file1->cd();
  h_cutflow->Fill(cutcounter);
  h_cutflow->GetXaxis()->SetBinLabel( cutcounter, cuts[cutcounter-1].c_str() );
  cutcounter++;
  ///THE END OF ************* Final k3pi selection  ****************** ///
  

  /// ************* Particle identification ****************** ///
  if (eltrkq1==eltrkq2)    {cleanup();return 0;}
  dirs[cutcounter-1]->cd();
    
  if (IS_MC){
    //if(true_4.P()  < 3. ||true_5.P() < 3.) {cleanup();return 0;}
    ((TH1D*)gDirectory->FindObject( h_mc_muvee_M->GetName()))->Fill( Muee_three_track.M());
    ((TH1D*)gDirectory->FindObject( h_mc_muvee_P->GetName()))->Fill( Muee_three_track.P());
    ((TH1D*)gDirectory->FindObject( h_mc_muvee_Pt->GetName()))->Fill( Muee_three_track.Pt());
    ((TH1D*)gDirectory->FindObject( h_mc_elsP->GetName()))->Fill( mc_ee.P());
    ((TH1D*)gDirectory->FindObject( h_mc_elsM->GetName()))->Fill( mc_ee.M());
    ((TH1D*)gDirectory->FindObject( h_mc_nu_M->GetName()))->Fill( Nu_true.M());
    ((TH1D*)gDirectory->FindObject( h_mc_nu_M2->GetName()))->Fill( Nu_true.M2());
    ((TH1D*)gDirectory->FindObject( h_mc_nu_P->GetName()))->Fill( Nu_true.P());
    ((TH1D*)gDirectory->FindObject( h_mc_nu_Pt->GetName()))->Fill(Nu_true.Pt());
  
    ((TH1D*)gDirectory->FindObject(h_true_M_0->GetName()))->Fill( true_0.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_0->GetName()))->Fill( true_0.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_0->GetName()))->Fill(true_0.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_0->GetName()))->Fill( true_0.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_1->GetName()))->Fill( true_1.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_1->GetName()))->Fill( true_1.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_1->GetName()))->Fill(true_1.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_1->GetName()))->Fill( true_1.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_2->GetName()))->Fill( true_2.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_2->GetName()))->Fill( true_2.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_2->GetName()))->Fill(true_2.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_2->GetName()))->Fill( true_2.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_3->GetName()))->Fill( true_3.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_3->GetName()))->Fill( true_3.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_3->GetName()))->Fill(true_3.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_3->GetName()))->Fill( true_3.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_4.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_4.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_4.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_4.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_5->GetName()))->Fill( true_5.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_5->GetName()))->Fill( true_5.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_5->GetName()))->Fill(true_5.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_5->GetName()))->Fill( true_5.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_6->GetName()))->Fill( true_6.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_6->GetName()))->Fill( true_6.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_6->GetName()))->Fill(true_6.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_6->GetName()))->Fill( true_6.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_7->GetName()))->Fill( true_7.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_7->GetName()))->Fill( true_7.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_7->GetName()))->Fill(true_7.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_7->GetName()))->Fill( true_7.E());
  
  }


  ((TH1D*)gDirectory->FindObject( h_vz_initial->GetName()))->Fill( sevt->vtx[0].z);
  ((TH1I*)gDirectory->FindObject(h_ntrack->GetName()))->Fill(ntrack);
  ((TH1I*)gDirectory->FindObject(h_nclust->GetName()))->Fill(nclust);

  ((TH1I*)gDirectory->FindObject(h_evttypebytrks->GetName()))->Fill(evttypebytrks);
  ((TH1I*)gDirectory->FindObject(h_nelectrack->GetName()))->Fill(nelectrack);
  ((TH1I*)gDirectory->FindObject(h_mutrack->GetName()))->Fill(nmutrack);
  ((TH2I*)gDirectory->FindObject(h_neltrk_mutrk->GetName()))->Fill(nelectrack, nmutrack);
  ((TH1I*)gDirectory->FindObject(h_mutrackq->GetName()))->Fill(muchtrkq);
  ((TH1I*)gDirectory->FindObject(h_e1trackq->GetName()))->Fill(eltrkq1);
  ((TH1I*)gDirectory->FindObject(h_e2trackq->GetName()))->Fill(eltrkq2);

  ((TH1F*)gDirectory->FindObject(h_el1trktime->GetName()))        ->Fill( el1trktime );
  ((TH1F*)gDirectory->FindObject(h_el1trkhodtime->GetName()))        ->Fill( el1trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_el1trktimediff->GetName()))    ->Fill( el1trktime - el1trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_ecltime1->GetName()))            ->Fill( el1cltime );
  ((TH1F*)gDirectory->FindObject(h_e1trkcltimediff->GetName()))    ->Fill( el1trktime- el1cltime );

  ((TH1F*)gDirectory->FindObject(h_el2trktime->GetName()))        ->Fill( el2trktime );
  ((TH1F*)gDirectory->FindObject(h_el2trkhodtime->GetName()))        ->Fill( el2trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_el2trktimediff->GetName()))    ->Fill( el2trktime - el2trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_ecltime2->GetName()))            ->Fill( el2cltime );
  ((TH1F*)gDirectory->FindObject(h_e2trkcltimediff->GetName()))    ->Fill( el2trktime-el2cltime );

  ((TH1F*)gDirectory->FindObject(h_e1e2timediff->GetName()))        ->Fill( el1trktime-el2trktime );
  ((TH1F*)gDirectory->FindObject(h_e1e2cltimediff->GetName()))    ->Fill( el1cltime-el2cltime );
  ((TH1F*)gDirectory->FindObject(h_e1e2hodtimediff->GetName()))        ->Fill( el1trkhodtime-el2trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_e1muhodtimediff->GetName()))        ->Fill( el1trkhodtime-muchtrkhodtime );
  ((TH1F*)gDirectory->FindObject(h_e2muhodtimediff->GetName()))        ->Fill( el2trkhodtime-muchtrkhodtime );

  ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_mu_track);
  ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_e1_track);
  ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_e2_track);

  
  ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_mu_track, (*muchlv).P());
  ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_e1_track, (*el1lv).P() );
  ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_e2_track, (*el2lv).P() );
  
  ((TH1F*)gDirectory->FindObject(h_eop_mu->GetName()))->Fill(eop_mu_track);
  ((TH2F*)gDirectory->FindObject(h_eop_mu_vs_p->GetName()))->Fill(eop_mu_track, (*muchlv).P());
  ((TH1F*)gDirectory->FindObject(h_eop_el->GetName()))->Fill(eop_e1_track);
  ((TH2F*)gDirectory->FindObject(h_eop_el_vs_p->GetName()))->Fill(eop_e1_track, (*el1lv).P());
  ((TH1F*)gDirectory->FindObject(h_eop_el->GetName()))->Fill(eop_e2_track);
  ((TH2F*)gDirectory->FindObject(h_eop_el_vs_p->GetName()))->Fill(eop_e2_track, (*el2lv).P());

  
  ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill( elp1_uncorr[0] );
  ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill( elp1_uncorr[1] );
  ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill( elv1_uncorr[0] );
  ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill( elv1_uncorr[1] );
  ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill( elp1_uncorr[0], elp1_uncorr[1] );
  ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill( elv1_uncorr[0], elv1_uncorr[1] );

  ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( el1lv->P() );
  ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( el2lv->P() );
  ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el1lv->E() );
  ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el2lv->E() );
  ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el1lv->P() );
  ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el2lv->P() );
  ((TH1F*)gDirectory->FindObject(h_e1m->GetName()))->Fill( el1lv->M() );
  ((TH1F*)gDirectory->FindObject(h_e2m->GetName()))->Fill( el2lv->M() );
  ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el1lv->Pt() );
  ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el2lv->Pt() );
  ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
  ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
  ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );
  ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);
  ((TH1F*)gDirectory->FindObject(h_angle_e1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()) );
  ((TH1F*)gDirectory->FindObject(h_angle_e1e2_mu->GetName()))->Fill( (*eelv).Vect().Angle(muchlv->Vect()) );

  ((TH1F*)gDirectory->FindObject(h_muchtrktime->GetName()))        ->Fill(muchtrktime);
  ((TH1F*)gDirectory->FindObject(h_muchtrkhodtime->GetName()))        ->Fill( muchtrkhodtime );
  ((TH1F*)gDirectory->FindObject(h_muchtrktimediff->GetName()))    ->Fill( muchtrktime-muchtrkhodtime );
  //((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( muchtrktime-muchcltime );

  ((TH1F*)gDirectory->FindObject(h_muchlvE->GetName()))->Fill( (*muchlv).E() );
  ((TH1F*)gDirectory->FindObject(h_muchlvP->GetName()))->Fill( (*muchlv).P() );
  ((TH1F*)gDirectory->FindObject(h_muchlvPt->GetName()))->Fill((*muchlv).Pt() );

  ((TH1F*)gDirectory->FindObject(h_KE->GetName()))->Fill((*pK).E());
  ((TH1F*)gDirectory->FindObject(h_KP->GetName()))->Fill((*pK).P());

  ((TH1F*)gDirectory->FindObject(h_nu_P->GetName()))->Fill(  (*neutrino).P() );
  ((TH1F*)gDirectory->FindObject(h_nu_E->GetName()))->Fill(  (*neutrino).E() );
  ((TH1F*)gDirectory->FindObject(h_nu_Pt->GetName()))->Fill( (*neutrino).Pt() );
  ((TH1F*)gDirectory->FindObject(h_nu_m->GetName()))->Fill(  (*neutrino).M() );
  ((TH1F*)gDirectory->FindObject(h_missing_mass->GetName()))->Fill( (*neutrino).M2() );

  if( muchtrkq+eltrkq1+eltrkq2 == 1){
    ((TH1F*)gDirectory->FindObject(h_missing_mass_kp->GetName()))->Fill( (*neutrino).M2() );
  } else {
    ((TH1F*)gDirectory->FindObject(h_missing_mass_km->GetName()))->Fill( (*neutrino).M2() );
  }

    
  ((TH1F*)gDirectory->FindObject(h_mueeE->GetName()))->Fill( (*muvee).E());
  ((TH1F*)gDirectory->FindObject(h_mueeP->GetName()))->Fill( (*muvee).P());
  ((TH1F*)gDirectory->FindObject(h_mueePt->GetName()))->Fill((*muvee).Pt());
  ((TH1F*)gDirectory->FindObject(h_mueeM->GetName()))->Fill ((*muvee).M());
  ((TH1F*)gDirectory->FindObject(h_3trk_invm->GetName()))->Fill ((*muvee).M2());
    
    
  ((TH1F*)gDirectory->FindObject(h_bx_much->GetName()))            ->Fill(muchp_uncorr[0]);
  ((TH1F*)gDirectory->FindObject(h_by_much->GetName()))            ->Fill(muchp_uncorr[1]);
  ((TH1F*)gDirectory->FindObject(h_bdxdz_much->GetName()))            ->Fill(muchv_uncorr[0]);
  ((TH1F*)gDirectory->FindObject(h_bdydz_much->GetName()))            ->Fill(muchv_uncorr[1]) ;
  ((TH1F*)gDirectory->FindObject(h_bx_vs_by_much->GetName()))        ->Fill(muchp_uncorr[0], muchp_uncorr[1]);
  ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_much->GetName()))->Fill(muchv_uncorr[0], muchv_uncorr[1]);

  ((TH1F*)gDirectory->FindObject(h_cda_ee->GetName()))->Fill(cda);
  ((TH1F*)gDirectory->FindObject(h_vx->GetName()))->Fill(vertex[0]);
  ((TH1F*)gDirectory->FindObject(h_vy->GetName()))->Fill(vertex[1]);
  ((TH1F*)gDirectory->FindObject(h_vz->GetName()))->Fill(vertex[2]);
  ((TH2F*)gDirectory->FindObject(h_vx_vs_vy->GetName()))->Fill(vertex[0],vertex[1]);

  ((TH1F*)gDirectory->FindObject(h_vx_e1much->GetName()))->Fill(vertex_e1much[0]);
  ((TH1F*)gDirectory->FindObject(h_vy_e1much->GetName()))->Fill(vertex_e1much[1]);
  ((TH1F*)gDirectory->FindObject(h_vz_e1much->GetName()))->Fill(vertex_e1much[2]);
  ((TH1F*)gDirectory->FindObject(h_cda_e1much->GetName()))->Fill(cda_e1much);
  ((TH1F*)gDirectory->FindObject(h_vz_e1much_diff->GetName()))->Fill( vertex[2] - vertex_e1much[2] );

  ((TH1F*)gDirectory->FindObject(h_vx_e2much->GetName()))->Fill(vertex_e2much[0]);
  ((TH1F*)gDirectory->FindObject(h_vy_e2much->GetName()))->Fill(vertex_e2much[1]);
  ((TH1F*)gDirectory->FindObject(h_vz_e2much->GetName()))->Fill(vertex_e2much[2]);
  ((TH1F*)gDirectory->FindObject(h_cda_e2much->GetName()))->Fill(cda_e2much);
  ((TH1F*)gDirectory->FindObject(h_vz_e2much_diff->GetName()))->Fill( vertex[2] - vertex_e2much[2] );
  ((TH1F*)gDirectory->FindObject(h_vz_mu1mu2_diff->GetName()))->Fill( vertex_e1much[2] - vertex_e2much[2] );

  
  ((TH1F*)gDirectory->FindObject(h_cda_e1much->GetName()))->Fill(cda_e1much);
  ((TH2F*)gDirectory->FindObject(h_vx_vs_vy_e1much->GetName()))->Fill(vertex_e1much[0],vertex_e1much[1]);
  ((TH1F*)gDirectory->FindObject(h_angle_e1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()) );
  ((TH1F*)gDirectory->FindObject(h_angle_e1e2_mu->GetName()))->Fill( (*eelv).Vect().Angle(muchlv->Vect()) );

  file1->cd();
  h_cutflow->Fill(cutcounter);
  h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());
  cutcounter++;
  ///THE END OF ************* Particle identification ****************** ///

  /// ************* Inv mass -140 ****************** ///
  
  if (eelvmass <= 0.140)  {cleanup();return 0;}
  if(eop_e2_track <=0.95 || eop_e1_track <=0.95 ) {cleanup();return 0;}
  //if (flux_id==1){
  //  if ((*muvee).M() < 0.470 || (*muvee).M() > 0.505 ) {cleanup();return 0;}
  //  if ((*muvee).P() < 54. || (*muvee).P() > 66 ) {cleanup();return 0;}
  //}
  dirs[cutcounter-1]->cd();
  ////COME HERE
  //cout << "---------BEGINNING OF EVENT ----------------"<< endl;
  //cout << "Event number =" << f_cout << "Mu index =" << muchtrkindx << "E1 index =" << eltrkindx1 << "E2 index =" << eltrkindx2 << endl;
  //cout << "Time stamp =" << sevt->timeStamp << endl;
  //cout << "Momentum of:" << endl;
  //cout << "Muon -" << sevt->track[muchtrkindx].p << endl;
  //cout << "Electron1 -" << sevt->track[eltrkindx1].p << endl;
  //cout << "Electron2 -" << sevt->track[eltrkindx2].p << endl;
  //f_cout++;
  //cout << "---------ENDING OF EVENT ----------------"<< endl;


  if (nEvt_counter==0){
    cout << "nEvt=" << sevt->nEvt << endl;
    first_Nevt=sevt->nEvt;
    v_first_e1= new TLorentzVector((*el1lv));
    v_first_e2= new TLorentzVector((*el2lv));
    v_first_mu= new TLorentzVector((*muchlv));
    f_Ntrack=sevt->Ntrack;
    f_NCl=sevt->Ncluster;
  } else {
    second_Nevt = sevt->nEvt;
    v_second_e1= new TLorentzVector((*el1lv));
    v_second_e2= new TLorentzVector((*el2lv));
    v_second_mu= new TLorentzVector((*muchlv));
    s_Ntrack=sevt->Ntrack;
    s_NCl=sevt->Ncluster;
  }
  //cout << "first_Nevt=" << first_Nevt << endl;
  //cout << " second_Nevt=" << second_Nevt << endl;
  int difference= second_Nevt - first_Nevt;
  
  if(nEvt_counter >= 0){

    if (fabs( difference ) <= 1.){
      cout << "Ntrack diff" << s_Ntrack - f_Ntrack << endl;
      cout << "NCluster diff" << s_NCl - f_NCl << endl;
      cout << "e1.P diff" << (*v_first_e1).P() - (*v_second_e1).P() << endl;
      cout << "e2.P diff" << (*v_first_e2).P() - (*v_second_e2).P() << endl;
      cout << "mu.P diff" << (*v_first_mu).P() - (*v_second_mu).P() << endl;
      cout << "e1.P diff" << (*v_first_e1).E() - (*v_second_e1).E() << endl;
      cout << "e2.P diff" << (*v_first_e2).E() - (*v_second_e2).E() << endl;
      cout << "mu.P diff" << (*v_first_mu).E() - (*v_second_mu).E() << endl;
      cout << "TUKA SUM" << endl;
      rep_counter++;
    }
    v_first_e1=v_second_e1;
    v_first_e2=v_second_e2;
    v_first_mu=v_second_mu;
    f_Ntrack  =s_Ntrack;
    f_NCl     =s_NCl;
    first_Nevt=second_Nevt;
  }
  
 

 
 
  
  nEvt_counter++;
  cout << "nEvt_counter=" << nEvt_counter << endl;
  cout << "rep_counter" << rep_counter << endl;

    
  if (IS_MC){
    //if( true_4.P() < 3. && true_5.P() < 3.) {cleanup();return 0;}
    ((TH1D*)gDirectory->FindObject( h_mc_muvee_M->GetName()))->Fill( Muee_three_track.M());
    ((TH1D*)gDirectory->FindObject( h_mc_muvee_P->GetName()))->Fill( Muee_three_track.P());
    ((TH1D*)gDirectory->FindObject( h_mc_muvee_Pt->GetName()))->Fill( Muee_three_track.Pt());
    ((TH1D*)gDirectory->FindObject( h_mc_elsP->GetName()))->Fill( mc_ee.P());
    ((TH1D*)gDirectory->FindObject( h_mc_elsM->GetName()))->Fill( mc_ee.M());
    ((TH1D*)gDirectory->FindObject( h_mc_nu_M->GetName()))->Fill( Nu_true.M());
    ((TH1D*)gDirectory->FindObject( h_mc_nu_M2->GetName()))->Fill( Nu_true.M2());
    ((TH1D*)gDirectory->FindObject( h_mc_nu_P->GetName()))->Fill( Nu_true.P());
    ((TH1D*)gDirectory->FindObject( h_mc_nu_Pt->GetName()))->Fill(Nu_true.Pt());
  
    ((TH1D*)gDirectory->FindObject(h_true_M_0->GetName()))->Fill( true_0.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_0->GetName()))->Fill( true_0.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_0->GetName()))->Fill(true_0.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_0->GetName()))->Fill( true_0.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_1->GetName()))->Fill( true_1.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_1->GetName()))->Fill( true_1.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_1->GetName()))->Fill(true_1.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_1->GetName()))->Fill( true_1.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_2->GetName()))->Fill( true_2.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_2->GetName()))->Fill( true_2.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_2->GetName()))->Fill(true_2.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_2->GetName()))->Fill( true_2.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_3->GetName()))->Fill( true_3.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_3->GetName()))->Fill( true_3.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_3->GetName()))->Fill(true_3.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_3->GetName()))->Fill( true_3.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_4.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_4.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_4.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_4.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_5->GetName()))->Fill( true_5.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_5->GetName()))->Fill( true_5.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_5->GetName()))->Fill(true_5.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_5->GetName()))->Fill( true_5.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_6->GetName()))->Fill( true_6.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_6->GetName()))->Fill( true_6.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_6->GetName()))->Fill(true_6.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_6->GetName()))->Fill( true_6.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_7->GetName()))->Fill( true_7.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_7->GetName()))->Fill( true_7.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_7->GetName()))->Fill(true_7.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_7->GetName()))->Fill( true_7.E());
  
  }

  ((TH1D*)gDirectory->FindObject( h_vz_initial->GetName()))->Fill( sevt->vtx[0].z);
  ((TH1I*)gDirectory->FindObject(h_evttypebytrks->GetName()))->Fill(evttypebytrks);
  ((TH1I*)gDirectory->FindObject(h_nelectrack->GetName()))->Fill(nelectrack);
  ((TH1I*)gDirectory->FindObject(h_mutrack->GetName()))->Fill(nmutrack);
  ((TH2I*)gDirectory->FindObject(h_neltrk_mutrk->GetName()))->Fill(nelectrack, nmutrack);
  ((TH1I*)gDirectory->FindObject(h_mutrackq->GetName()))->Fill(muchtrkq);
  ((TH1I*)gDirectory->FindObject(h_e1trackq->GetName()))->Fill(eltrkq1);
  ((TH1I*)gDirectory->FindObject(h_e2trackq->GetName()))->Fill(eltrkq2);

  ((TH1F*)gDirectory->FindObject(h_el1trktime->GetName()))        ->Fill( el1trktime );
  ((TH1F*)gDirectory->FindObject(h_el1trkhodtime->GetName()))        ->Fill( el1trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_el1trktimediff->GetName()))    ->Fill( el1trktime - el1trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_ecltime1->GetName()))            ->Fill( el1cltime );
  ((TH1F*)gDirectory->FindObject(h_e1trkcltimediff->GetName()))    ->Fill( el1trktime- el1cltime );

  ((TH1F*)gDirectory->FindObject(h_el2trktime->GetName()))        ->Fill( el2trktime );
  ((TH1F*)gDirectory->FindObject(h_el2trkhodtime->GetName()))        ->Fill( el2trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_el2trktimediff->GetName()))    ->Fill( el2trktime - el2trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_ecltime2->GetName()))            ->Fill( el2cltime );
  ((TH1F*)gDirectory->FindObject(h_e2trkcltimediff->GetName()))    ->Fill( el2trktime-el2cltime );

  ((TH1F*)gDirectory->FindObject(h_e1e2timediff->GetName()))        ->Fill( el1trktime-el2trktime );
  ((TH1F*)gDirectory->FindObject(h_e1e2cltimediff->GetName()))    ->Fill( el1cltime-el2cltime );
  ((TH1F*)gDirectory->FindObject(h_e1e2hodtimediff->GetName()))        ->Fill( el1trkhodtime-el2trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_e1muhodtimediff->GetName()))        ->Fill( el1trkhodtime-muchtrkhodtime );
  ((TH1F*)gDirectory->FindObject(h_e2muhodtimediff->GetName()))        ->Fill( el2trkhodtime-muchtrkhodtime );


  ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_mu_track);
  ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_e1_track);
  ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_e2_track);

  
  ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_mu_track, (*muchlv).P());
  ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_e1_track, (*el1lv).P() );
  ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_e2_track, (*el2lv).P() );
  
  ((TH1F*)gDirectory->FindObject(h_eop_mu->GetName()))->Fill(eop_mu_track);
  ((TH2F*)gDirectory->FindObject(h_eop_mu_vs_p->GetName()))->Fill(eop_mu_track, (*muchlv).P());
  ((TH1F*)gDirectory->FindObject(h_eop_el->GetName()))->Fill(eop_e1_track);
  ((TH2F*)gDirectory->FindObject(h_eop_el_vs_p->GetName()))->Fill(eop_e1_track, (*el1lv).P());
  ((TH1F*)gDirectory->FindObject(h_eop_el->GetName()))->Fill(eop_e2_track);
  ((TH2F*)gDirectory->FindObject(h_eop_el_vs_p->GetName()))->Fill(eop_e2_track, (*el2lv).P());


       
  ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill( elp1_uncorr[0] );
  ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill( elp1_uncorr[1] );
  ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill( elv1_uncorr[0] );
  ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill( elv1_uncorr[1] );
  ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill( elp1_uncorr[0], elp1_uncorr[1] );
  ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill( elv1_uncorr[0], elv1_uncorr[1] );

  ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( el1lv->P() );
  ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( el2lv->P() );
  ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el1lv->E() );
  ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el2lv->E() );
  ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el1lv->P() );
  ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el2lv->P() );
  ((TH1F*)gDirectory->FindObject(h_e1m->GetName()))->Fill( el1lv->M() );
  ((TH1F*)gDirectory->FindObject(h_e2m->GetName()))->Fill( el2lv->M() );
  ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el1lv->Pt() );
  ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el2lv->Pt() );
  ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
  ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
  ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );
  ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);
  ((TH1F*)gDirectory->FindObject(h_angle_e1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()) );
  ((TH1F*)gDirectory->FindObject(h_angle_e1e2_mu->GetName()))->Fill( (*eelv).Vect().Angle(muchlv->Vect()) );
  
  ((TH1F*)gDirectory->FindObject(h_muchtrktime->GetName()))        ->Fill(muchtrktime);
  ((TH1F*)gDirectory->FindObject(h_muchtrkhodtime->GetName()))        ->Fill( muchtrkhodtime );
  ((TH1F*)gDirectory->FindObject(h_muchtrktimediff->GetName()))    ->Fill( muchtrktime-muchtrkhodtime );
  //((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( muchtrktime-muchcltime );

  ((TH1F*)gDirectory->FindObject(h_muchlvE->GetName()))->Fill( (*muchlv).E() );
  ((TH1F*)gDirectory->FindObject(h_muchlvP->GetName()))->Fill( (*muchlv).P() );
  ((TH1F*)gDirectory->FindObject(h_muchlvPt->GetName()))->Fill((*muchlv).Pt() );

  ((TH1F*)gDirectory->FindObject(h_KE->GetName()))->Fill((*pK).E());
  ((TH1F*)gDirectory->FindObject(h_KP->GetName()))->Fill((*pK).P());
    
  ((TH1F*)gDirectory->FindObject(h_nu_P->GetName()))->Fill(  (*neutrino).P() );
  ((TH1F*)gDirectory->FindObject(h_nu_E->GetName()))->Fill(  (*neutrino).E() );
  ((TH1F*)gDirectory->FindObject(h_nu_Pt->GetName()))->Fill( (*neutrino).Pt() );
  ((TH1F*)gDirectory->FindObject(h_nu_m->GetName()))->Fill(  (*neutrino).M() );
  ((TH1F*)gDirectory->FindObject(h_missing_mass->GetName()))->Fill( (*neutrino).M2() );

  if( muchtrkq+eltrkq1+eltrkq2 == 1){
    ((TH1F*)gDirectory->FindObject(h_missing_mass_kp->GetName()))->Fill( (*neutrino).M2() );
  } else {
    ((TH1F*)gDirectory->FindObject(h_missing_mass_km->GetName()))->Fill( (*neutrino).M2() );
  }

    
  ((TH1F*)gDirectory->FindObject(h_mueeE->GetName()))->Fill( (*muvee).E());
  ((TH1F*)gDirectory->FindObject(h_mueeP->GetName()))->Fill( (*muvee).P());
  ((TH1F*)gDirectory->FindObject(h_mueePt->GetName()))->Fill((*muvee).Pt());
  ((TH1F*)gDirectory->FindObject(h_mueeM->GetName()))->Fill ((*muvee).M());
  ((TH1F*)gDirectory->FindObject(h_3trk_invm->GetName()))->Fill ((*muvee).M2());
    
    
  ((TH1F*)gDirectory->FindObject(h_bx_much->GetName()))            ->Fill(muchp_uncorr[0]);
  ((TH1F*)gDirectory->FindObject(h_by_much->GetName()))            ->Fill(muchp_uncorr[1]);
  ((TH1F*)gDirectory->FindObject(h_bdxdz_much->GetName()))            ->Fill(muchv_uncorr[0]);
  ((TH1F*)gDirectory->FindObject(h_bdydz_much->GetName()))            ->Fill(muchv_uncorr[1]) ;
  ((TH1F*)gDirectory->FindObject(h_bx_vs_by_much->GetName()))        ->Fill(muchp_uncorr[0], muchp_uncorr[1]);
  ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_much->GetName()))->Fill(muchv_uncorr[0], muchv_uncorr[1]);

  ((TH1F*)gDirectory->FindObject(h_cda_ee->GetName()))->Fill(cda);
  ((TH1F*)gDirectory->FindObject(h_vx->GetName()))->Fill(vertex[0]);
  ((TH1F*)gDirectory->FindObject(h_vy->GetName()))->Fill(vertex[1]);
  ((TH1F*)gDirectory->FindObject(h_vz->GetName()))->Fill(vertex[2]);
  ((TH2F*)gDirectory->FindObject(h_vx_vs_vy->GetName()))->Fill(vertex[0],vertex[1]);

  ((TH1F*)gDirectory->FindObject(h_vx_e1much->GetName()))->Fill(vertex_e1much[0]);
  ((TH1F*)gDirectory->FindObject(h_vy_e1much->GetName()))->Fill(vertex_e1much[1]);
  ((TH1F*)gDirectory->FindObject(h_vz_e1much->GetName()))->Fill(vertex_e1much[2]);
  ((TH1F*)gDirectory->FindObject(h_cda_e1much->GetName()))->Fill(cda_e1much);
  ((TH1F*)gDirectory->FindObject(h_vz_e1much_diff->GetName()))->Fill( vertex[2] - vertex_e1much[2] );

  ((TH1F*)gDirectory->FindObject(h_vx_e2much->GetName()))->Fill(vertex_e2much[0]);
  ((TH1F*)gDirectory->FindObject(h_vy_e2much->GetName()))->Fill(vertex_e2much[1]);
  ((TH1F*)gDirectory->FindObject(h_vz_e2much->GetName()))->Fill(vertex_e2much[2]);
  ((TH1F*)gDirectory->FindObject(h_cda_e2much->GetName()))->Fill(cda_e2much);
  ((TH1F*)gDirectory->FindObject(h_vz_e2much_diff->GetName()))->Fill( vertex[2] - vertex_e2much[2] );
  ((TH1F*)gDirectory->FindObject(h_vz_mu1mu2_diff->GetName()))->Fill( vertex_e1much[2] - vertex_e2much[2] );

  
  ((TH1F*)gDirectory->FindObject(h_cda_e1much->GetName()))->Fill(cda_e1much);
  ((TH2F*)gDirectory->FindObject(h_vx_vs_vy_e1much->GetName()))->Fill(vertex_e1much[0],vertex_e1much[1]);
  
  
  
  file1->cd();
  h_cutflow->Fill(cutcounter);
  h_cutflow->GetXaxis()->SetBinLabel( cutcounter, cuts[cutcounter-1].c_str() );
  cutcounter++;
  ///THE END OF ************* Inv mass -140 ****************** ///

    
  // ~~~ CUT ~~~ CUT ~~~ CUT ~~~
  //"e1e2_hod_Tmatch"
  //if ( el1cltime-el2cltime < -3 || el1cltime-el2cltime > 3 ) {cleanup();return 0;}
  //if (  el1trkhodtime-el2trkhodtime < -3 || el1trkhodtime-el2trkhodtime > 3 ) {cleanup();return 0;}
  //if(IS_DATA)    {
  //    //after final selection this difference has a maximum at -5.5 and sigma 1.5
  //    //5ns is according to Manuel's thesis
  //    if( fabs(el1trktime - el1trkhodtime)>12. )     return 0;
  //    if( fabs(el2trktime - el2trkhodtime)>12. )     return 0;
  //}
  //else    {
  //    if( el1trkhodtime<136. || el1trkhodtime>146. )     {cleanup();return 0;}
  //    if( el2trkhodtime<136. || el2trkhodtime>146. )     {cleanup();return 0;}
  //    //        if( el1trkhodtime<80. || el1trkhodtime>150. )     {cleanup();return 0;}
  //    //        if( el2trkhodtime<80. || el2trkhodtime>150. )     {cleanup();return 0;}
  //}
  //
  //dirs[cutcounter-1]->cd();
    
  //file1->cd();
  //h_cutflow->Fill(cutcounter);
  //h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"e1e2_hod_Tmatch");
  //cutcounter++;


   
  /// **********  charged muon
  //cout << "Mu charge track index = " << muchtrkindx << endl;
  //printf("abcog_params.cogX1p=%f,abcog_params.cogX1n=%f\n",abcog_params.cogX1p,abcog_params.cogX1n);
  //printf("abcog_params.cogX4p=%f,abcog_params.cogX4n=%f\n",abcog_params.cogX4p,abcog_params.cogX4n);

  // ~~~ CUT ~~~ CUT ~~~ CUT ~~~
  //"mu_track"
  //if (muchtrkindx==-1)    {cleanup();return 0;}
  ////if (nmutrack!=1)    {cleanup();return 0;}
  ////
  //dirs[cutcounter-1]->cd();
  //
  //
  //
  //
  //
  //file1->cd();
  //h_cutflow->Fill(cutcounter);
  //h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"mu_track");
  //cutcounter++;
    
  //SO FAR !!!!


  // ~~~ CUT ~~~ CUT ~~~ CUT ~~~
  //"mu_hod_Tmatch"
  //if(IS_DATA)    {
  //    //after final selection this difference has a maximum at -5.5 and sigma 1.5
  //    if( fabs(muchtrktime-muchtrkhodtime) > 12. )    return 0;  //5 is from Manuel's thesis
  //}
  //else    {
  //    if( muchtrkhodtime<135. || muchtrkhodtime>148. )    return 0;
  //}
  ////
  //dirs[cutcounter-1]->cd();
  //file1->cd();
  //h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"mu_hod_Tmatch");
  //h_cutflow->Fill(cutcounter);
  //cutcounter++;

    


  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------
  //                            electron and positron matching
  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------


  // ~~~ CUT ~~~ CUT ~~~ CUT ~~~
  //"e1e2_cl_Tmatch")
  //    if( fabs((long double)(el1cltime-el2cltime))>15. )        {cleanup();return 0;}
  //if( fabs((long double)(el1cltime-el2cltime))>10. )        {cleanup();return 0;}
  //
  //dirs[cutcounter-1]->cd();
  
  //file1->cd();
  //h_cutflow->Fill(cutcounter);
  //h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"e1e2_cl_Tmatch");
  //cutcounter++;



  // ~~~ CUT ~~~ CUT ~~~ CUT ~~~
  //"e1e2_DCHinnerR")
  //if( pow(elp1_uncorr[0],2.) + pow(elp1_uncorr[1],2.) < 196. )    {cleanup();return 0;} //14*14 /it was 12*12/
  //if( pow(elp2_uncorr[0],2.) + pow(elp2_uncorr[1],2.) < 196. )    {cleanup();return 0;}
  ////
  //dirs[cutcounter-1]->cd();
  //
  //
  //file1->cd();
  //h_cutflow->Fill(cutcounter);
  //h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"e1e2_DCHinnerR");
  //cutcounter++;



    
  // ~~~ CUT ~~~ CUT ~~~ CUT ~~~
  //"e1e2_DCHouterR"
  //    if( pow(elp1_uncorr[0],2.) + pow(elp1_uncorr[1],2.) > 13225. )    {cleanup();return 0;} //115*115
  //    if( pow(elp2_uncorr[0],2.) + pow(elp2_uncorr[1],2.) > 13225. )    {cleanup();return 0;}
  //if( pow(elp1_uncorr[0],2.) + pow(elp1_uncorr[1],2.) > 12100. )    {cleanup();return 0;} //110*110
  //if( pow(elp2_uncorr[0],2.) + pow(elp2_uncorr[1],2.) > 12100. )    {cleanup();return 0;}
  //
  //dirs[cutcounter-1]->cd();
    
  
   
  //file1->cd();
  //h_cutflow->Fill(cutcounter);
  //h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"e1e2_DCHouterR");
  //cutcounter++;


  //    cout<<"3"<<endl;


  // ~~~ CUT ~~~ CUT ~~~ CUT ~~~
  //"e1e2_vtx_found"

  //if(vertex_OK!=1)    return 0;

  //dirs[cutcounter-1]->cd();  

  //file1->cd();
  //h_cutflow->Fill(cutcounter);
  //h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"e1e2_vtx_found");
  //cutcounter++;



    
  // ~~~ CUT ~~~ CUT ~~~ CUT ~~~
   
  //    if( vertex[2]<-2000. || vertex[2]>8000. )    {cleanup();return 0;}
  // if (eelvmass <= 0.140)  {cleanup();return 0;}
 
    
  //dirs[cutcounter-1]->cd();    
  //
  //
  //file1->cd();
  //h_cutflow->Fill(cutcounter);
  //h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"ee_mass_140_and_pcut");
  //cutcounter++;
    
  //"e1e2_vtx_good" 
  //if( vertex[2]<-1500. || vertex[2]>7000. )    {cleanup();return 0;}
  //
  //dirs[cutcounter-1]->cd();

  //file1->cd();
  //h_cutflow->Fill(cutcounter);
  //h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"e1e2_vtx_good");
  //cutcounter++;




  // ~~~ CUT ~~~ CUT ~~~ CUT ~~~
  //"e1e2_dist_DCH1"
  //    if(eldistDCH1plane<2.)        return 0;
  //if(eldistDCH1plane<1.)        return 0;
  ////
  //dirs[cutcounter-1]->cd();
  //file1->cd();
  //h_cutflow->Fill(cutcounter);
  //h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"e1e2_dist_DCH1");
  //cutcounter++;


  //    cout<<"6"<<endl;


  // ~~~ CUT ~~~ CUT ~~~ CUT ~~~
  //"e1e2_dist_LKr"
  //    if( eldistLKrplane<15. )        return 0;
  //if( eldistLKrplane<20. )        {cleanup();return 0;}
  ////
  //dirs[cutcounter-1]->cd();
  // file1->cd();
  //h_cutflow->Fill(cutcounter);
  //h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"e1e2_dist_LKr");
  //cutcounter++;

				  

  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------
  //                        charged muon and electron-positron matching
  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------


  // ~~~ CUT ~~~ CUT ~~~ CUT ~~~
  //"mu_ee_Tmatch"
  //if( fabs((long double)(el1trktime-muchtrktime)) > 10. )    {cleanup();return 0;}
  //if( fabs((long double)(el1cltime-muchcltime))   > 10. )    {cleanup();return 0;}
    
  //
  //dirs[cutcounter-1]->cd();
  //
  //file1->cd();
  //h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"mu_ee_Tmatch")
  //h_cutflow->Fill(cutcounter);
  //cutcounter++;



  // ~~~ CUT ~~~ CUT ~~~ CUT ~~~
  //"mu_vtx_good"
  //if( !muchee_vertex_OK)    {cleanup();return 0;}
  //    if( fabs((long double)(vertex[2] - vertex_e1pch[2]))  > 1500. )    {cleanup();return 0;} //old value 300
  //if( vertex_e1much[2]<-1500.  || vertex_e1much[2]>7000. )    {cleanup();return 0;}
  //
  // dirs[cutcounter-1]->cd();

  //
  //
  //file1->cd();
  //h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"mu_vtx_good")
  //h_cutflow->Fill(cutcounter);
  //cutcounter++;



  // ~~~ CUT ~~~ CUT ~~~ CUT ~~~
  //"mu_ee_dist_LKr  & hodtime cut"
  //if(el1trkhodtime-muchtrkhodtime < -2 || el1trkhodtime-muchtrkhodtime > 2 )    {cleanup();return 0;}
  //if(el2trkhodtime-muchtrkhodtime < -2 || el2trkhodtime-muchtrkhodtime > 2 )    {cleanup();return 0;}
  //if(muche1distLKrplane<30.)    {cleanup();return 0;}
  //if(muche2distLKrplane<30.)    {cleanup();return 0;}
  ////
  //dirs[cutcounter-1]->cd();
  //
  //
  //((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill(elp1_uncorr[0]);
  //((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill(elp1_uncorr[1]);
  //((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill(elv1_uncorr[0]);
  //((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill(elv1_uncorr[1]);
  //((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill(elp1_uncorr[0], elp1_uncorr[1] );
  //((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill(elv1_uncorr[0], elv1_uncorr[1] );
  //
  //((TH1F*)gDirectory->FindObject(h_bx_much->GetName()))            ->Fill(muchp_uncorr[0]);
  //((TH1F*)gDirectory->FindObject(h_by_much->GetName()))            ->Fill(muchp_uncorr[1]);
  //((TH1F*)gDirectory->FindObject(h_bdxdz_much->GetName()))            ->Fill(muchv_uncorr[0]);
  //((TH1F*)gDirectory->FindObject(h_bdydz_much->GetName()))            ->Fill(muchv_uncorr[1]) ;
  //((TH1F*)gDirectory->FindObject(h_bx_vs_by_much->GetName()))        ->Fill(muchp_uncorr[0], muchp_uncorr[1]);
  //((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_much->GetName()))->Fill(  muchv_uncorr[0], muchv_uncorr[1]);
  //
  //((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( (*el1lv).P() );
  //((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( (*el2lv).P() );
  //((TH1F*)gDirectory->FindObject(h_e1m->GetName()))->Fill( el1lv->M() );
  //((TH1F*)gDirectory->FindObject(h_e2m->GetName()))->Fill( el2lv->M() );
  //
  //((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill(  (*el1lv).E() );
  //((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill(  (*el2lv).E() );
  //((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill(  (*el1lv).P() );
  //((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill(  (*el2lv).P() );
  //((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( (*el1lv).Pt() );
  //((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( (*el2lv).Pt() );
  //
  //  
  //((TH1F*)gDirectory->FindObject(h_muchlvE->GetName()))->Fill( (*muchlv).E() );
  //((TH1F*)gDirectory->FindObject(h_muchlvP->GetName()))->Fill( (*muchlv).P() );
  //((TH1F*)gDirectory->FindObject(h_muchlvPt->GetName()))->Fill( (*muchlv).Pt() );
  //
  //((TH1I*)gDirectory->FindObject(h_evttypebytrks->GetName()))->Fill(evttypebytrks);
  //((TH1I*)gDirectory->FindObject(h_nelectrack->GetName()))->Fill(nelectrack);
  //((TH1I*)gDirectory->FindObject(h_mutrack->GetName()))->Fill(nmutrack);
  //((TH1I*)gDirectory->FindObject(h_mutrackq->GetName()))->Fill(muchtrkq);
  //((TH2I*)gDirectory->FindObject(h_neltrk_mutrk->GetName()))->Fill(nelectrack, nmutrack);
  //
  //((TH1F*)gDirectory->FindObject(h_eop_mu->GetName()))->Fill(sevt->cluster[sevt->track[muchtrkindx].iClus].energy/sevt->track[muchtrkindx].p);
  //((TH2F*)gDirectory->FindObject(h_eop_mu_vs_p->GetName()))->Fill(sevt->cluster[sevt->track[muchtrkindx].iClus].energy/sevt->track[muchtrkindx].p, sevt->track[muchtrkindx].p);
  //
  //((TH1F*)gDirectory->FindObject(h_eop_el->GetName()))->Fill(sevt->cluster[sevt->track[eltrkindx1].iClus].energy/sevt->track[eltrkindx1].p);
  //((TH2F*)gDirectory->FindObject(h_eop_el_vs_p->GetName()))->Fill(sevt->cluster[sevt->track[eltrkindx1].iClus].energy/sevt->track[eltrkindx1].p, sevt->track[eltrkindx1].p);
  //
  //((TH1F*)gDirectory->FindObject(h_eop_el->GetName()))->Fill(sevt->cluster[sevt->track[eltrkindx2].iClus].energy/sevt->track[eltrkindx2].p);
  //((TH2F*)gDirectory->FindObject(h_eop_el_vs_p->GetName()))->Fill(sevt->cluster[sevt->track[eltrkindx2].iClus].energy/sevt->track[eltrkindx2].p, sevt->track[eltrkindx2].p);
  //	
  //
  //
  //((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
  //((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
  //((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );
  //
  //((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);
  //((TH1F*)gDirectory->FindObject(h_nu_P->GetName()))->Fill(  (*neutrino).P() );
  //((TH1F*)gDirectory->FindObject(h_nu_E->GetName()))->Fill(  (*neutrino).E() );
  //((TH1F*)gDirectory->FindObject(h_nu_Pt->GetName()))->Fill( (*neutrino).Pt() );
  //((TH1F*)gDirectory->FindObject(h_nu_m->GetName()))->Fill(  (*neutrino).M() );
  //
  //((TH1F*)gDirectory->FindObject(h_missing_mass->GetName()))->Fill( (*neutrino).M2() );
  //((TH1F*)gDirectory->FindObject(h_mueeE->GetName()))->Fill((*muvee).E());
  //((TH1F*)gDirectory->FindObject(h_mueeP->GetName()))->Fill((*muvee).P());
  //((TH1F*)gDirectory->FindObject(h_mueePt->GetName()))->Fill((*muvee).Pt());
  //((TH1F*)gDirectory->FindObject(h_mueeM->GetName()))->Fill ((*muvee).M());
  //((TH1F*)gDirectory->FindObject(h_3trk_invm->GetName()))->Fill( (*muvee).M2());
  //
  //((TH1F*)gDirectory->FindObject(h_el1trktime->GetName()))        ->Fill( el1trktime );
  //((TH1F*)gDirectory->FindObject(h_el1trkhodtime->GetName()))        ->Fill( el1trkhodtime );
  //((TH1F*)gDirectory->FindObject(h_el1trktimediff->GetName()))    ->Fill( el1trktime - el1trkhodtime );
  //((TH1F*)gDirectory->FindObject(h_ecltime1->GetName()))            ->Fill( el1cltime );
  //((TH1F*)gDirectory->FindObject(h_e1trkcltimediff->GetName()))    ->Fill( el1trktime- el1cltime );
  //
  //((TH1F*)gDirectory->FindObject(h_el2trktime->GetName()))        ->Fill( el2trktime );
  //((TH1F*)gDirectory->FindObject(h_el2trkhodtime->GetName()))        ->Fill( el2trkhodtime );
  //((TH1F*)gDirectory->FindObject(h_el2trktimediff->GetName()))    ->Fill( el2trktime - el2trkhodtime );
  //((TH1F*)gDirectory->FindObject(h_ecltime2->GetName()))            ->Fill( el2cltime );
  //((TH1F*)gDirectory->FindObject(h_e2trkcltimediff->GetName()))    ->Fill( el2trktime-el2cltime );
  //
  //((TH1F*)gDirectory->FindObject(h_e1e2timediff->GetName()))        ->Fill( el1trktime-el2trktime );
  //((TH1F*)gDirectory->FindObject(h_e1e2cltimediff->GetName()))    ->Fill( el1cltime-el2cltime );
  //
  //((TH1F*)gDirectory->FindObject(h_muchtrktime->GetName()))        ->Fill(muchtrktime);
  //((TH1F*)gDirectory->FindObject(h_muchtrkhodtime->GetName()))        ->Fill( muchtrkhodtime );
  //((TH1F*)gDirectory->FindObject(h_muchcltime->GetName()))            ->Fill(muchcltime);
  //((TH1F*)gDirectory->FindObject(h_muchtrktimediff->GetName()))    ->Fill( muchtrktime-muchtrkhodtime );
  //((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( muchtrktime-muchcltime );
  //((TH1F*)gDirectory->FindObject(h_e1e2hodtimediff->GetName()))        ->Fill( el1trkhodtime-el2trkhodtime );
  //((TH1F*)gDirectory->FindObject(h_e1muhodtimediff->GetName()))        ->Fill( el1trkhodtime-muchtrkhodtime );
  //((TH1F*)gDirectory->FindObject(h_e2muhodtimediff->GetName()))        ->Fill( el2trkhodtime-muchtrkhodtime );
  //
  //
  //((TH1F*)gDirectory->FindObject(h_cda_ee->GetName()))->Fill(cda);
  //((TH1F*)gDirectory->FindObject(h_vx->GetName()))->Fill(vertex[0]);
  //((TH1F*)gDirectory->FindObject(h_vy->GetName()))->Fill(vertex[1]);
  //((TH1F*)gDirectory->FindObject(h_vz->GetName()))->Fill(vertex[2]);
  //((TH2F*)gDirectory->FindObject(h_vx_vs_vy->GetName()))->Fill(vertex[0],vertex[1]);
  //
  //((TH1F*)gDirectory->FindObject(h_vx_e1much->GetName()))->Fill(vertex_e1much[0]);
  //((TH1F*)gDirectory->FindObject(h_vy_e1much->GetName()))->Fill(vertex_e1much[1]);
  //((TH1F*)gDirectory->FindObject(h_vz_e1much->GetName()))->Fill(vertex_e1much[2]);
  //((TH1F*)gDirectory->FindObject(h_cda_e1much->GetName()))->Fill(cda_e1much);
  //((TH1F*)gDirectory->FindObject(h_vz_e1much_diff->GetName()))->Fill( vertex[2] - vertex_e1much[2] );
  //
  //((TH1F*)gDirectory->FindObject(h_cda_e1much->GetName()))->Fill(cda_e1much);
  //((TH2F*)gDirectory->FindObject(h_vx_vs_vy_e1much->GetName()))->Fill(vertex_e1much[0],vertex_e1much[1]);
  //
  //
  //((TH1F*)gDirectory->FindObject(h_muche1distDHC1plane->GetName()))->Fill( muche1distDHC1plane );
  //((TH1F*)gDirectory->FindObject(h_muche2distDHC1plane->GetName()))->Fill( muche2distDHC1plane );
  //
  //
  //((TH1F*)gDirectory->FindObject(h_muche1distLKrplane->GetName()))->Fill(  muche1distLKrplane );
  //((TH1F*)gDirectory->FindObject(h_muche2distLKrplane->GetName()))->Fill(  muche2distLKrplane );
  //
  
  //
  //
  ////    ((TH1F*)gDirectory->FindObject(h_gE->GetName()))->Fill( (*glv).E() );
  ////    ((TH1F*)gDirectory->FindObject(h_gP->GetName()))->Fill( (*glv).P() );
  ////    ((TH1F*)gDirectory->FindObject(h_gPt->GetName()))->Fill( (*glv).Pt() );
  ////    ((TH1F*)gDirectory->FindObject(h_gtrkRadDHC1->GetName()))->Fill( gtrkRadDHC1 );
  ////
  ////    ((TH1F*)gDirectory->FindObject(h_pchlvE->GetName()))->Fill( (*pchlv).E() );
  ////    ((TH1F*)gDirectory->FindObject(h_pchlvP->GetName()))->Fill( (*pchlv).P() );
  ////    ((TH1F*)gDirectory->FindObject(h_pchlvPt->GetName()))->Fill( (*pchlv).Pt() );
  ////
  ////
  ////    ((TH1F*)gDirectory->FindObject(h_angle_e1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()) );
  ////    ((TH1F*)gDirectory->FindObject(h_angle_ph_e1e2->GetName()))->Fill( (*eelv).Vect().Angle(glv->Vect()) );
  ////    ((TH1F*)gDirectory->FindObject(h_angle_pch_pi0->GetName()))->Fill( (*pi0lv).Vect().Angle(pchlv->Vect()) );
  ////    ((TH1F*)gDirectory->FindObject(h_a_e1e2_phe1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()), (*eelv).Vect().Angle(glv->Vect()) );
  ////    ((TH1F*)gDirectory->FindObject(h_a_e1e2_pchpi0->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()), (*pi0lv).Vect().Angle(pchlv->Vect()) );
  ////    ((TH1F*)gDirectory->FindObject(h_a_phe1e2_pchpi0->GetName()))->Fill( (*eelv).Vect().Angle(glv->Vect()), (*pi0lv).Vect().Angle(pchlv->Vect()) );
  ////
  ////    ((TH1F*)gDirectory->FindObject(h_cda_ee->GetName()))->Fill(cda);
  ////    ((TH1F*)gDirectory->FindObject(h_vx->GetName()))->Fill(vertex[0]);
  ////    ((TH1F*)gDirectory->FindObject(h_vy->GetName()))->Fill(vertex[1]);
  ////    ((TH1F*)gDirectory->FindObject(h_vz->GetName()))->Fill(vertex[2]);
  ////    ((TH2F*)gDirectory->FindObject(h_vx_vs_vy->GetName()))->Fill(vertex[0],vertex[1]);
  ////
  ////    ((TH1F*)gDirectory->FindObject(h_cda_e1pch->GetName()))->Fill(cda_e1pch);
  ////    ((TH1F*)gDirectory->FindObject(h_vx_e1pch->GetName()))->Fill(vertex_e1pch[0]);
  ////    ((TH1F*)gDirectory->FindObject(h_vy_e1pch->GetName()))->Fill(vertex_e1pch[1]);
  ////    ((TH1F*)gDirectory->FindObject(h_vz_e1pch->GetName()))->Fill(vertex_e1pch[2]);
  ////    ((TH2F*)gDirectory->FindObject(h_vx_vs_vy_e1pch->GetName()))->Fill(vertex_e1pch[0],vertex_e1pch[1]);
  ////
  ////    ((TH1F*)gDirectory->FindObject(h_vz_e1pch_diff->GetName()))->Fill( vertex[2] - vertex_e1pch[2] );
  //
  ////    ((TH1F*)gDirectory->FindObject(h_pche1distDHC1plane->GetName()))->Fill( pche1distDHC1plane );
  ////    ((TH1F*)gDirectory->FindObject(h_pche2distDHC1plane->GetName()))->Fill( pche2distDHC1plane );
  ////    ((TH1F*)gDirectory->FindObject(h_pche1distLKrplane->GetName()))->Fill( pche1distLKrplane );
  ////    ((TH1F*)gDirectory->FindObject(h_pche2distLKrplane->GetName()))->Fill( pche2distLKrplane );
  ////    ((TH1F*)gDirectory->FindObject(h_gpchdistLKrplane->GetName()))->Fill( gpchdistLKrplane );
  //
  ////
  ////((TH1F*)gDirectory->FindObject(h_pi0E->GetName()))->Fill(  pi0lv->E() );
  ////((TH1F*)gDirectory->FindObject(h_pi0P->GetName()))->Fill(  pi0lv->P() );
  ////((TH1F*)gDirectory->FindObject(h_pi0Pt->GetName()))->Fill( pi0lv->Pt() );
  ////((TH1F*)gDirectory->FindObject(h_mofpi0->GetName()))->Fill( (*pi0lv).M() );
  ////
  ////((TH1F*)gDirectory->FindObject(h_pi0Pcm->GetName()))->Fill( pi0lvcm->P() );
  ////((TH1F*)gDirectory->FindObject(h_mofpi0diff->GetName()))->Fill( (*pi0lvcm).M()-(*pi0lv).M() );
  ////
  ////
  ////((TH1F*)gDirectory->FindObject(h_gmueeE->GetName()))->Fill(  ((*glv)+(*muvee)).E() );
  ////((TH1F*)gDirectory->FindObject(h_gmueeP->GetName()))->Fill(  ((*glv)+(*muvee)).P() );
  ////((TH1F*)gDirectory->FindObject(h_gmueePt->GetName()))->Fill( ((*glv)+(*muvee)).Pt());
  ////((TH1F*)gDirectory->FindObject(h_gmueeM->GetName()))->Fill(  ((*glv)+(*muvee)).M() );
  //
  //file1->cd();
  //h_cutflow->Fill(cutcounter);
  //h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"mu_ee_dist_LKr");
  //cutcounter++;
    
  /// ************* electron E/p between 0.9 and 0.95 (at least one) ****************** ///
  //"muee energy > 70GeV test"  Mu transverse momentum > 0.9
  //if( (*muvee).E() < 70 )        {cleanup();return 0;}
  //if( (*muchlv).Pt() > 0.9 )        {
  //if( sevt->cluster[sevt->track[eltrkindx1].iClus].energy/sevt->track[eltrkindx1].p < 0.95 || sevt->cluster[sevt->track[eltrkindx2].iClus].energy/sevt->track[eltrkindx2].p < 0.95 )        {
  //
  if(mc_ee.M() >= 0.1 )  {cleanup();return 0;}
  dirs[cutcounter-1]->cd();
    
  if (IS_MC){
    //if( true_4.P() < 3. && true_5.P() < 3.) {cleanup();return 0;}
    ((TH1D*)gDirectory->FindObject( h_mc_muvee_M->GetName()))->Fill( Muee_three_track.M());
    ((TH1D*)gDirectory->FindObject( h_mc_muvee_P->GetName()))->Fill( Muee_three_track.P());
    ((TH1D*)gDirectory->FindObject( h_mc_muvee_Pt->GetName()))->Fill( Muee_three_track.Pt());
    ((TH1D*)gDirectory->FindObject( h_mc_elsP->GetName()))->Fill( mc_ee.P());
    ((TH1D*)gDirectory->FindObject( h_mc_elsM->GetName()))->Fill( mc_ee.M());
    ((TH1D*)gDirectory->FindObject( h_mc_nu_M->GetName()))->Fill( Nu_true.M());
    ((TH1D*)gDirectory->FindObject( h_mc_nu_M2->GetName()))->Fill( Nu_true.M2());
    ((TH1D*)gDirectory->FindObject( h_mc_nu_P->GetName()))->Fill( Nu_true.P());
    ((TH1D*)gDirectory->FindObject( h_mc_nu_Pt->GetName()))->Fill(Nu_true.Pt());
  
    ((TH1D*)gDirectory->FindObject(h_true_M_0->GetName()))->Fill( true_0.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_0->GetName()))->Fill( true_0.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_0->GetName()))->Fill(true_0.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_0->GetName()))->Fill( true_0.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_1->GetName()))->Fill( true_1.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_1->GetName()))->Fill( true_1.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_1->GetName()))->Fill(true_1.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_1->GetName()))->Fill( true_1.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_2->GetName()))->Fill( true_2.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_2->GetName()))->Fill( true_2.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_2->GetName()))->Fill(true_2.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_2->GetName()))->Fill( true_2.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_3->GetName()))->Fill( true_3.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_3->GetName()))->Fill( true_3.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_3->GetName()))->Fill(true_3.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_3->GetName()))->Fill( true_3.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_4.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_4.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_4.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_4.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_5->GetName()))->Fill( true_5.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_5->GetName()))->Fill( true_5.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_5->GetName()))->Fill(true_5.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_5->GetName()))->Fill( true_5.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_6->GetName()))->Fill( true_6.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_6->GetName()))->Fill( true_6.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_6->GetName()))->Fill(true_6.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_6->GetName()))->Fill( true_6.E());
    ((TH1D*)gDirectory->FindObject(h_true_M_7->GetName()))->Fill( true_7.M());
    ((TH1D*)gDirectory->FindObject(h_true_P_7->GetName()))->Fill( true_7.P());
    ((TH1D*)gDirectory->FindObject(h_true_Pt_7->GetName()))->Fill(true_7.Pt());
    ((TH1D*)gDirectory->FindObject(h_true_E_7->GetName()))->Fill( true_7.E());
  
  }

  ((TH1I*)gDirectory->FindObject(h_ntrack->GetName()))->Fill(ntrack);
  ((TH1I*)gDirectory->FindObject(h_nclust->GetName()))->Fill(nclust);
  
  ((TH1I*)gDirectory->FindObject(h_evttypebytrks->GetName()))->Fill(evttypebytrks);
  ((TH1I*)gDirectory->FindObject(h_nelectrack->GetName()))->Fill(nelectrack);
  ((TH1I*)gDirectory->FindObject(h_mutrack->GetName()))->Fill(nmutrack);
  ((TH2I*)gDirectory->FindObject(h_neltrk_mutrk->GetName()))->Fill(nelectrack, nmutrack);
  ((TH1I*)gDirectory->FindObject(h_mutrackq->GetName()))->Fill(muchtrkq);
  ((TH1I*)gDirectory->FindObject(h_e1trackq->GetName()))->Fill(eltrkq1);
  ((TH1I*)gDirectory->FindObject(h_e2trackq->GetName()))->Fill(eltrkq2);
    
  ((TH1F*)gDirectory->FindObject(h_el1trktime->GetName()))        ->Fill( el1trktime );
  ((TH1F*)gDirectory->FindObject(h_el1trkhodtime->GetName()))        ->Fill( el1trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_el1trktimediff->GetName()))    ->Fill( el1trktime - el1trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_ecltime1->GetName()))            ->Fill( el1cltime );
  ((TH1F*)gDirectory->FindObject(h_e1trkcltimediff->GetName()))    ->Fill( el1trktime- el1cltime );

  ((TH1F*)gDirectory->FindObject(h_el2trktime->GetName()))        ->Fill( el2trktime );
  ((TH1F*)gDirectory->FindObject(h_el2trkhodtime->GetName()))        ->Fill( el2trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_el2trktimediff->GetName()))    ->Fill( el2trktime - el2trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_ecltime2->GetName()))            ->Fill( el2cltime );
  ((TH1F*)gDirectory->FindObject(h_e2trkcltimediff->GetName()))    ->Fill( el2trktime-el2cltime );

  ((TH1F*)gDirectory->FindObject(h_e1e2timediff->GetName()))        ->Fill( el1trktime-el2trktime );
  ((TH1F*)gDirectory->FindObject(h_e1e2cltimediff->GetName()))    ->Fill( el1cltime-el2cltime );
  ((TH1F*)gDirectory->FindObject(h_e1e2hodtimediff->GetName()))        ->Fill( el1trkhodtime-el2trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_e1muhodtimediff->GetName()))        ->Fill( el1trkhodtime-muchtrkhodtime );
  ((TH1F*)gDirectory->FindObject(h_e2muhodtimediff->GetName()))        ->Fill( el2trkhodtime-muchtrkhodtime );

  
  ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_mu_track);
  ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_e1_track);
  ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_e2_track);

  
  ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_mu_track, (*muchlv).P());
  ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_e1_track, (*el1lv).P() );
  ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_e2_track, (*el2lv).P() );
  
  ((TH1F*)gDirectory->FindObject(h_eop_mu->GetName()))->Fill(eop_mu_track);
  ((TH2F*)gDirectory->FindObject(h_eop_mu_vs_p->GetName()))->Fill(eop_mu_track, (*muchlv).P());
  ((TH1F*)gDirectory->FindObject(h_eop_el->GetName()))->Fill(eop_e1_track);
  ((TH2F*)gDirectory->FindObject(h_eop_el_vs_p->GetName()))->Fill(eop_e1_track, (*el1lv).P());
  ((TH1F*)gDirectory->FindObject(h_eop_el->GetName()))->Fill(eop_e2_track);
  ((TH2F*)gDirectory->FindObject(h_eop_el_vs_p->GetName()))->Fill(eop_e2_track, (*el2lv).P());

  ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill( elp1_uncorr[0] );
  ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill( elp1_uncorr[1] );
  ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill( elv1_uncorr[0] );
  ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill( elv1_uncorr[1] );
  ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill( elp1_uncorr[0], elp1_uncorr[1] );
  ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill( elv1_uncorr[0], elv1_uncorr[1] );

  ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( el1lv->P() );
  ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( el2lv->P() );
  ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el1lv->E() );
  ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el2lv->E() );
  ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el1lv->P() );
  ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el2lv->P() );
  ((TH1F*)gDirectory->FindObject(h_e1m->GetName()))->Fill( el1lv->M() );
  ((TH1F*)gDirectory->FindObject(h_e2m->GetName()))->Fill( el2lv->M() );
  ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el1lv->Pt() );
  ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el2lv->Pt() );
  ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
  ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
  ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );
  ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);
  ((TH1F*)gDirectory->FindObject(h_angle_e1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()) );
  ((TH1F*)gDirectory->FindObject(h_angle_e1e2_mu->GetName()))->Fill( (*eelv).Vect().Angle(muchlv->Vect()) );

  ((TH1F*)gDirectory->FindObject(h_muchtrktime->GetName()))        ->Fill(muchtrktime);
  ((TH1F*)gDirectory->FindObject(h_muchtrkhodtime->GetName()))        ->Fill( muchtrkhodtime );
  ((TH1F*)gDirectory->FindObject(h_muchtrktimediff->GetName()))    ->Fill( muchtrktime-muchtrkhodtime );
  //((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( muchtrktime-muchcltime );

  ((TH1F*)gDirectory->FindObject(h_muchlvE->GetName()))->Fill( (*muchlv).E() );
  ((TH1F*)gDirectory->FindObject(h_muchlvP->GetName()))->Fill( (*muchlv).P() );
  ((TH1F*)gDirectory->FindObject(h_muchlvPt->GetName()))->Fill((*muchlv).Pt() );

  ((TH1F*)gDirectory->FindObject(h_KE->GetName()))->Fill((*pK).E());
  ((TH1F*)gDirectory->FindObject(h_KP->GetName()))->Fill((*pK).P());
    
  ((TH1F*)gDirectory->FindObject(h_nu_P->GetName()))->Fill(  (*neutrino).P() );
  ((TH1F*)gDirectory->FindObject(h_nu_E->GetName()))->Fill(  (*neutrino).E() );
  ((TH1F*)gDirectory->FindObject(h_nu_Pt->GetName()))->Fill( (*neutrino).Pt() );
  ((TH1F*)gDirectory->FindObject(h_nu_m->GetName()))->Fill(  (*neutrino).M() );
  ((TH1F*)gDirectory->FindObject(h_missing_mass->GetName()))->Fill( (*neutrino).M2() );

  if( muchtrkq+eltrkq1+eltrkq2 == 1){
    ((TH1F*)gDirectory->FindObject(h_missing_mass_kp->GetName()))->Fill( (*neutrino).M2() );
  } else {
    ((TH1F*)gDirectory->FindObject(h_missing_mass_km->GetName()))->Fill( (*neutrino).M2() );
  }

  ((TH1F*)gDirectory->FindObject(h_mueeE->GetName()))->Fill( (*muvee).E());
  ((TH1F*)gDirectory->FindObject(h_mueeP->GetName()))->Fill( (*muvee).P());
  ((TH1F*)gDirectory->FindObject(h_mueePt->GetName()))->Fill((*muvee).Pt());
  ((TH1F*)gDirectory->FindObject(h_mueeM->GetName()))->Fill ((*muvee).M());
  ((TH1F*)gDirectory->FindObject(h_3trk_invm->GetName()))->Fill ((*muvee).M2());
    
    
  ((TH1F*)gDirectory->FindObject(h_bx_much->GetName()))            ->Fill(muchp_uncorr[0]);
  ((TH1F*)gDirectory->FindObject(h_by_much->GetName()))            ->Fill(muchp_uncorr[1]);
  ((TH1F*)gDirectory->FindObject(h_bdxdz_much->GetName()))            ->Fill(muchv_uncorr[0]);
  ((TH1F*)gDirectory->FindObject(h_bdydz_much->GetName()))            ->Fill(muchv_uncorr[1]) ;
  ((TH1F*)gDirectory->FindObject(h_bx_vs_by_much->GetName()))        ->Fill(muchp_uncorr[0], muchp_uncorr[1]);
  ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_much->GetName()))->Fill(muchv_uncorr[0], muchv_uncorr[1]);

  ((TH1F*)gDirectory->FindObject(h_cda_ee->GetName()))->Fill(cda);
  ((TH1F*)gDirectory->FindObject(h_vx->GetName()))->Fill(vertex[0]);
  ((TH1F*)gDirectory->FindObject(h_vy->GetName()))->Fill(vertex[1]);
  ((TH1F*)gDirectory->FindObject(h_vz->GetName()))->Fill(vertex[2]);
  ((TH2F*)gDirectory->FindObject(h_vx_vs_vy->GetName()))->Fill(vertex[0],vertex[1]);

  ((TH1F*)gDirectory->FindObject(h_vx_e1much->GetName()))->Fill(vertex_e1much[0]);
  ((TH1F*)gDirectory->FindObject(h_vy_e1much->GetName()))->Fill(vertex_e1much[1]);
  ((TH1F*)gDirectory->FindObject(h_vz_e1much->GetName()))->Fill(vertex_e1much[2]);
  ((TH1F*)gDirectory->FindObject(h_cda_e1much->GetName()))->Fill(cda_e1much);
  ((TH1F*)gDirectory->FindObject(h_vz_e1much_diff->GetName()))->Fill( vertex[2] - vertex_e1much[2] );

  ((TH1F*)gDirectory->FindObject(h_vx_e2much->GetName()))->Fill(vertex_e2much[0]);
  ((TH1F*)gDirectory->FindObject(h_vy_e2much->GetName()))->Fill(vertex_e2much[1]);
  ((TH1F*)gDirectory->FindObject(h_vz_e2much->GetName()))->Fill(vertex_e2much[2]);
  ((TH1F*)gDirectory->FindObject(h_cda_e2much->GetName()))->Fill(cda_e2much);
  ((TH1F*)gDirectory->FindObject(h_vz_e2much_diff->GetName()))->Fill( vertex[2] - vertex_e2much[2] );
  ((TH1F*)gDirectory->FindObject(h_vz_mu1mu2_diff->GetName()))->Fill( vertex_e1much[2] - vertex_e2much[2] );

  
    
  ((TH1F*)gDirectory->FindObject(h_cda_e1much->GetName()))->Fill(cda_e1much);
  ((TH2F*)gDirectory->FindObject(h_vx_vs_vy_e1much->GetName()))->Fill(vertex_e1much[0],vertex_e1much[1]);

    
  // }

  file1->cd();
  h_cutflow->Fill(cutcounter);
  h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str()); 
  cutcounter++;
  ///THE END OF ************* electron E/p between 0.9 and 0.95 (at least one) ****************** ///
  //Mu transverse momentum > 0.9


  

  //------------------------------------------------------------------------------------------------
  //                                        clean up
  cleanup();
  return 0;
}




//Calculates and returns blue tube corrections for slopes of tracks on x-y plane
void applyBlueTubeCorrection(superCmpEvent *sevt, int itrk, double pKaon[3], double vKaon[3], double* trkSlopes, double* trkMidPoints)    {

  double bdxdz_track = sevt->track[itrk].bdxdz;    //Note - here i access same info again, to be changed with a variable!
  double bdydz_track = sevt->track[itrk].bdydz;
  double bdzdz_track = 1.0;
  double bx_track = sevt->track[itrk].bx;         // x coordinate before magnet (genauer: direkt (ca. 8cm) vor DCH1)
  double by_track = sevt->track[itrk].by;         // y coordinate before magnet (genauer: direkt (ca. 8cm) vor DCH1)
  double bz_track =  Geom->DCH.bz;                 // z before magnet

  //cout << "you are stupid !!!!!!!!!!!!!"  << IS_DATA << endl;
  
  double v1[3] = {0.};
  double p1[3] = {0.};
  v1[0] = bdxdz_track;
  v1[1] = bdydz_track;
  v1[2] = bdzdz_track;
  p1[0] = bx_track;
  p1[1] = by_track;
  p1[2] = bz_track;

  double cda=0.;
  bool vertexK_OK=false;
  double vertexK[3] = {0.};
  // Call routine for vertex calculation
  if ( closap_double_(p1, pKaon, v1, vKaon, &cda, vertexK) )
    vertexK_OK = TRUE;
  else
    vertexK_OK = FALSE;

  // Calculate track midpoint coordinates at z = (vertex[2] + DCHbz)/2
  double midpoint[3];
  midpoint[2] = (vertexK[2] + bz_track)/2.;
  midpoint[0] = bx_track - bdxdz_track*(bz_track-midpoint[2]);
  midpoint[1] = by_track - bdydz_track*(bz_track-midpoint[2]);

  //INPUTS:
  // the routine needs a vertex to start with computations,
  //I give here the intersection between track and beam axis (defined by cda point)
  int   nchar   = sevt->track[itrk].q;       // charge of the track
  float tmom    = sevt->track[itrk].p;       // momentum from compact
  float Vxyz[3] = {0.};
  float vpnt[2] = {0.};
  float vdir[2] = {0.};
  Vxyz[0] = vertexK[0];
  Vxyz[1] = vertexK[1];              // floats needed for routine
  Vxyz[2] = vertexK[2];
  vpnt[0] = bx_track;        //sevt->track[i].bx;      // bx from compact
  vpnt[1] = by_track;        //sevt->track[i].by;      // by from compact
  vdir[0] = bdxdz_track;    //sevt->track[i].bdxdz;   // slopes from compact
  vdir[1] = bdydz_track;    //sevt->track[i].bdydz;   // slopes from compact

  blue_tack_(&nchar, &tmom, Vxyz, vpnt, vdir);

  //OUTPUT:
  // as output the routine overwrites vpnt + vdir with bluefield corrected values
  bdxdz_track = vdir[0];
  bdydz_track = vdir[1];

  //Function output
  trkSlopes[0] = bdxdz_track;
  trkSlopes[1] = bdydz_track;
  trkSlopes[2] = bdzdz_track;
  trkMidPoints[0] = midpoint[0];
  trkMidPoints[1] = midpoint[1];
  trkMidPoints[2] = midpoint[2];

}




void cleanup()    {

  if(el1lv)             delete el1lv;     el1lv =0;
  if(el2lv)             delete el2lv;     el2lv =0;
  if(eelv)             delete eelv;      eelv =0;
  if(muchlv)             delete muchlv;     muchlv =0;
  if(muvee)             delete muvee;    muvee =0;
  if(pK)      delete pK; pK=0;
  // if(glv)             delete glv;       glv =0;
  //if(Klv)             delete Klv;       Klv =0;
  //if(el1lvcm)         delete el1lvcm;   el1lvcm =0;
  // if(el2lvcm)         delete el2lvcm;   el2lvcm =0;
  //if(glvcm)             delete glvcm;     glvcm =0;
  // if(pi0lvcm)         delete pi0lvcm;   pi0lvcm =0;
  // if(pi0lv)             delete pi0lv;     pi0lv =0;
  //if(pi0lvKcm)         delete pi0lvKcm;  pi0lvKcm =0;
  //if(pchlvKcm)         delete pchlvKcm;  pchlvKcm =0;
  //if(KlvKcm)             delete KlvKcm;    KlvKcm =0;

  if(trkCorrSlope)     delete trkCorrSlope;     trkCorrSlope=NULL;
  if(trkCorrMidPoint) delete trkCorrMidPoint; trkCorrMidPoint=NULL;
  if(clpos)             delete clpos;             clpos=NULL;

  for (unsigned int i = trkCorrSlopesV.size()-1; i>-1; i--) {
    if (trkCorrSlopesV[i])         delete trkCorrSlopesV[i];
    if (trkCorrMidPointsV[i])     delete trkCorrMidPointsV[i];
  }
  trkCorrSlopesV.clear();
  trkCorrMidPointsV.clear();

  for (unsigned int i = clcorrpos.size()-1; i>-1; i--) {
    if (clcorrpos[i])         delete clcorrpos[i];
  }
  clcorrpos.clear();

  return;
}
