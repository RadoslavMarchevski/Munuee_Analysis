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
TVector3* Kmomentum  =NULL;
TVector3* trkCorrSlope     = NULL;
TVector3* trkCorrMidPoint     = NULL;

vector<TVector3*> clcorrpos;
TVector3* clpos = NULL;

TLorentzVector *el1lv=NULL, *el2lv=NULL, *eelv=NULL;
TLorentzVector *el1lv_pi=NULL, *el2lv_pi=NULL, *eelv_pi=NULL;
TLorentzVector *pi1lv=NULL, *pi2lv=NULL,*pi3lv=NULL;
TLorentzVector *pchlv=NULL, *muchlv=NULL , *muvee=NULL;
TLorentzVector *muvee_pi=NULL ,*muvee_piee=NULL,*muchlv_pi=NULL;
TLorentzVector *pK=NULL, *neutrino=NULL,*k3pi_miss_mass=NULL ;
TLorentzVector *mu_nu=NULL;
TLorentzVector *pieelv=NULL,*k3pic=NULL;
TLorentzVector *glv =NULL;
TLorentzVector *pi0lv=NULL;
TLorentzVector *el1lvcm=NULL, *el2lvcm=NULL,  *glvcm=NULL, *pi0lvcm=NULL;
TLorentzVector *Klv=NULL, *pi0lvKcm=NULL, *pchlvKcm=NULL, *muchlvKcm=NULL, *KlvKcm=NULL,*muveeKcm=NULL;
TLorentzVector *K4Momentum=NULL;
TLorentzVector *v_first_e1=NULL;
TLorentzVector *v_second_e1=NULL;
TLorentzVector *v_first_e2=NULL;
TLorentzVector *v_second_e2=NULL;
TLorentzVector *v_first_mu=NULL;
TLorentzVector *v_second_mu=NULL;
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
  double vKaon[3] = {0.,0.,0.};
  double pKaon[3] = {0.,0.,0.};

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
  int    PU_Q2_Counter =      0;
  int    PU_Q1_Counter =      0;
  int    PU_Q1_NAKL_Counter = 0;
  int    PU_mu1_Counter=      0;

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

  double norm1=0, norm2=0,normmu=0;
  double p1x=0,p1y=0, p1z=0, p2x=0,p2y=0, p2z=0, el1E=0, el2E=0, pK1x=0,pK1y=0, pK1z=0,pKE=0,pmux=0,pmuy=0,pmuz=0;
  double el1E_pi=0,el2E_pi=0;
  double normpi1=0, normpi2=0,normpi3=0;
  double pi1x=0,pi1y=0, pi1z=0, pi2x=0,pi2y=0, pi2z=0, pi1E=0, pi2E=0, pi3x=0,pi3y=0, pi3z=0,pi3E=0;

  int    gclindx=-1;
  double gclenergy=0.;
  double gcltime=0.;
  double gclx=0.;
  double gcly=0.;
  double gclz=0.;
  double gtrkDHC1x, gtrkDHC1y, gtrkRadDHC1 ;


  double pchtrkhodtime, pchE , muchE , muchE_pi, muchtrkhodtime;

  double eldistDCH1plane, el1xLKrplane, el1yLKrplane, el2xLKrplane, el2yLKrplane, eldistLKrplane;
  //  double pidistDCH1plane, pi1xLKrplane, pi1yLKrplane, pi2xLKrplane, pi2yLKrplane, pidistLKrplane;

  double gel1distLKrplane, gel2distLKrplane;
  double pche1distDHC1plane, pche2distDHC1plane, pche1distLKrplane, pche2distLKrplane;
  double muche1distDHC1plane, muche2distDHC1plane, muche1distLKrplane, muche2distLKrplane;

  double pchxLKrplane, pchyLKrplane, gpchdistLKrplane;
  double muchxLKrplane, muchyLKrplane, gmuchdistLKrplane;
  double Dgvtx_x,Dgvtx_y, Dgvtx_z, gnorm,  gpx, gpy, gpz;

  unsigned int vertex_OK ,k3pic12_vertex_OK , k3pic23_vertex_OK, k3pic13_vertex_OK ;
  double vertex[3]={0.},vertex_k3pic12[3]={0.}, vertex_k3pic23[3]={0.},vertex_k3pic13[3]={0.};
  double cda, cda_k3pic12, cda_k3pic23, cda_k3pic13;

  unsigned int muchee_vertex_OK,muchee2_vertex_OK;

  double vertex_e2much[3]={0.}, vertex_e1much[3]={0.};
  double cda_e2much, cda_e1much;

  double elp1_uncorr[3]={0.}, elp1[3]={0.}, elv1_uncorr[3]={0.}, elv1[3]={0.}, elp2_uncorr[3]={0.},elp2[3]={0.}, elv2_uncorr[3]={0.},elv2[3]={0.}, pchp_uncorr[3]={0.},pchp[3]={0.}, pchv_uncorr[3]={0.}, pchv[3]={0.} ,muchp_uncorr[3]={0.},muchp[3]={0.}, muchv_uncorr[3]={0.}, muchv[3]={0.};
  double pilp1_uncorr[3]={0.}, pilp1[3]={0.}, pilv1_uncorr[3]={0.}, pilv1[3]={0.}, pilp2_uncorr[3]={0.},pilp2[3]={0.}, pilv2_uncorr[3]={0.},pilv2[3]={0.},pilp3_uncorr[3]={0.}, pilp3[3]={0.}, pilv3_uncorr[3]={0.}, pilv3[3]={0.};

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

  bool LKrCalCorr=true;
  bool CorrAlphaBeta=true;
  bool Ke2L3trig=false;
  bool Km2L3trig=false;



  int       cutcounter=1;
  int    ntrack = sevt->Ntrack; // number of tracks
  int    nclust = sevt->Ncluster; // number of tracks
  int    nvtx   = sevt->Nvtx;   // number of vtx
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
  double BeamMomentum = 0; //Kaon
  double Beamdxdz, Beamdydz;

  //        if (SQL_Database)    {
  if (Kcharge >= 0)                // K+ beam

    {
      beamCorrX     = abcog_params.pkdxdzp;
      beamOffsetX = abcog_params.pkxoffp;
      beamCorrY     = abcog_params.pkdydzp;
      beamOffsetY = abcog_params.pkyoffp;
      BeamMomentum = abcog_params.pkp * (1 + abcog_params.beta);
    }
  else  // K- beam
    {
      beamCorrX     = abcog_params.pkdxdzm;
      beamOffsetX = abcog_params.pkxoffm;
      beamCorrY     = abcog_params.pkdydzm;
      beamOffsetY = abcog_params.pkyoffm;
      BeamMomentum = abcog_params.pkm * (1 + abcog_params.beta);
      cout << beamCorrX << beamOffsetX << beamCorrY<< beamOffsetY<< endl;
    }
  //        }

  TVector3 Kmomentum (beamCorrX*BeamMomentum,beamCorrY*BeamMomentum,TMath::Sqrt(1. - pow(beamCorrX,2) - pow(beamCorrY,2)));
  K4Momentum = new TLorentzVector(Kmomentum[0],Kmomentum[1],Kmomentum[2],TMath::Sqrt(massKaonC*massKaonC + Kmomentum[0]*Kmomentum[0]+Kmomentum[2]*Kmomentum[2]+Kmomentum[2]*Kmomentum[2]) );

  //  ->SetVect(Kmomentum);
  //K4Momentum->SetE(TMath::Sqrt(massKaonC*massKaonC + Kmomentum*Kmomentum) );
 printf("abcog_params.cogX1p=%f, abcog_params.cogX1n=%f\n",abcog_params.cogX1p,abcog_params.cogX1n);
 //printf("abcog_params.cogX4p=%f,abcog_params.cogX4n=%f\n",abcog_params.cogX4p,abcog_params.cogX4n);
 //printf("abcog_params.cogY1p=%f, abcog_params.cogY1n=%f\n",abcog_params.cogY1p,abcog_params.cogY1n);
 //printf("abcog_params.cogY4p=%f,abcog_params.cogY4n=%f\n",abcog_params.cogY4p,abcog_params.cogY4n);
 ////
 //printf("abcog_params.pkp=%f, abcog_params.pkn=%f\n",abcog_params.pkp,abcog_params.pkm);
 //printf("abcog_params.pkdydzp=%f,abcog_params.pkdydzn=%f\n",abcog_params.pkdydzp,abcog_params.pkdydzm);
 //printf("abcog_params.pkdxdzp=%f, abcog_params.pkdxdzn=%f\n",abcog_params.pkdxdzp,abcog_params.pkdxdzm);
 return 0;
 //printf("abcog_params.cogY4p=%f,abcog_params.cogY4n=%f\n",abcog_params.cogY4p,abcog_params.cogY4n);

  vKaon[0] = beamCorrX;
  vKaon[1] = beamCorrY;
  vKaon[2] = 1;
  pKaon[0] = beamOffsetX + beamCorrX*DCHbz;
  pKaon[1] = beamOffsetY + beamCorrY*DCHbz;
  pKaon[2] = DCHbz;
   // <----------------------------------------------------------------------------------------------

  //void MarioKaon::Calculate4Momentum()
  //{
  //  double BeamMomentum = 0; //Kaon
 //  double Beamdxdz, Beamdydz;
  //
  //  if (evt->track[0].q < 0) // negative track == K- beam
  //    {
  //	Beamdxdz = abcog_params.pkdxdzm;
  //	Beamdydz = abcog_params.pkdydzm;
  //	BeamMomentum = abcog_params.pkm * (1 + abcog_params.beta);
  //    } else // positive track == K+ beam
  //    {
  //	Beamdxdz = abcog_params.pkdxdzp;
  //	Beamdydz = abcog_params.pkdydzp;
  //	BeamMomentum = abcog_params.pkp * (1 + abcog_params.beta);
  //    }
  //  y[3];
  //  y[0] = Beamdxdz; y[1] = Beamdydz; y[2] = TMath::Sqrt(1. - pow(Beamdxdz,2) - pow(Beamdydz,2));
  //
  //  TVector3 KMomentum(y[0]*BeamMomentum,y[1]*BeamMomentum,y[2]*BeamMomentum);
  //  K4Momentum.SetVect(KMomentum);
  //  K4Momentum.SetE(TMath::Sqrt(KaonMass*KaonMass + KMomentum*KMomentum) );
  //}



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
  //for (int i = 0; i < nclust; ++i) {
  //
  //  clenergy = sevt->cluster[i].energy;
  //  clx= sevt->cluster[i].x;
  //  cly= sevt->cluster[i].y;
  //  double x0 = clx;
  //  double y0 = cly;
  //
  //  clz = LKrz + 16.5 + 4.3*log(clenergy);
  //
  //  if(IS_DATA)    {
  //    clx = (x0 + 0.0136 + y0 * 0.87e-3 ) * ( 1 + (clz-LKrz)/10998.);
  //    cly = (y0 + 0.300  + x0 * 0.87e-3 ) * ( 1 + (clz-LKrz)/10998.);
  //  }
  //  else    {
  //    clx = (x0 - 0.013 ) * ( 1 + (clz-LKrz)/10998.);
  //    cly = y0              * ( 1 + (clz-LKrz)/10998.);
  //  }
  //
  //  clpos = new TVector3(clx, cly, clz);
  //  clcorrpos.push_back(clpos);
  //}

  //HistnoBursts();
  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------
  //                                        object selection
  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------


  // ~~~ NO CUT ~~~ NO CUT ~~~ NO CUT ~~~
  // ************* Initial cuts ****************** ///
  dirs[cutcounter-1]->cd();
  // cout << evt->Npart << endl;

  //  ((TH1D*)gDirectory->FindObject( h_vz_initial->GetName()))->Fill( sevt->vtx[0].z);

  if (IS_DATA)
    if( sbur->BadB.Lkr == 1||
	sbur->BadB.Dch == 1||
	//sbur->BadB.Nut == 1||
	sbur->BadB.Mbx == 1||
	sbur->BadB.Muv == 1||
	sbur->BadB.HodC== 1||
	//sbur->BadB.Clk == 1||
	//sbur->BadB.Kab == 1||
	sbur->BadB.Phys== 1){
      sbur->BadB.Skip=1;
    }

  if (IS_MC){
    ((TH1F*)gDirectory->FindObject( h_vz_initial->GetName()))->Fill( k_dvtx[2]);
    if (k_dvtx[2] < -2000 || k_dvtx[2] > 8000) {cleanup();return 0;}
  }
  h_ntrack->Fill(ntrack);
  h_nclust->Fill(nclust);





  if (IS_MC){
    ////if(true_4.P()  < 3. || true_5.P() < 3.) {cleanup();return 0;}
    ((TH1F*)gDirectory->FindObject(h_mc_Npart_tr->GetName()))->Fill(Npart);
    if (Npart >=4 ){

      ((TH1D*)gDirectory->FindObject( h_mc_muvee_M->GetName()))->Fill( Muee_three_track.M());
      ((TH1D*)gDirectory->FindObject( h_mc_muvee_P->GetName()))->Fill( Muee_three_track.P());
      ((TH1D*)gDirectory->FindObject( h_mc_muvee_Pt->GetName()))->Fill( Muee_three_track.Pt());
      ((TH1D*)gDirectory->FindObject( h_mc_elsP->GetName()))->Fill( mc_ee.P());
      ((TH1D*)gDirectory->FindObject( h_mc_elsM->GetName()))->Fill( mc_ee.M());
      ((TH1D*)gDirectory->FindObject( h_mc_nu_M->GetName()))->Fill( Nu_true.M());
      ((TH1D*)gDirectory->FindObject( h_mc_nu_M2->GetName()))->Fill( Nu_true.M2());
      ((TH1D*)gDirectory->FindObject( h_mc_nu_P->GetName()))->Fill( Nu_true.P());
      ((TH1D*)gDirectory->FindObject( h_mc_nu_Pt->GetName()))->Fill(Nu_true.Pt());
      ((TH1D*)gDirectory->FindObject(h_true_M_0->GetName()))->Fill( true_V[0]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_0->GetName()))->Fill( true_V[0]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_0->GetName()))->Fill(true_V[0]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_0->GetName()))->Fill( true_V[0]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_1->GetName()))->Fill( true_V[1]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_1->GetName()))->Fill( true_V[1]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_1->GetName()))->Fill(true_V[1]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_1->GetName()))->Fill( true_V[1]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_2->GetName()))->Fill( true_V[3]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_2->GetName()))->Fill( true_V[3]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_2->GetName()))->Fill(true_V[3]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_2->GetName()))->Fill( true_V[3]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_3->GetName()))->Fill( true_V[3]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_3->GetName()))->Fill( true_V[3]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_3->GetName()))->Fill(true_V[3]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_3->GetName()))->Fill( true_V[3]->E());

      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt0->GetName()))->Fill(pvtx[0]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt1->GetName()))->Fill(pvtx[1]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt2->GetName()))->Fill(pvtx[2]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt3->GetName()))->Fill(pvtx[3]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt0->GetName()))->Fill(dvtx[0]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt1->GetName()))->Fill(dvtx[1]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt2->GetName()))->Fill(dvtx[2]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt3->GetName()))->Fill(dvtx[3]);
       }
    if(Npart==5){
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt4->GetName()))->Fill(pvtx[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt4->GetName()))->Fill(dvtx[4]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt4->GetName()))->Fill(pType[4]);
      ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_V[4]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_V[4]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_V[4]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_V[4]->E());

    }
    if(Npart==6){
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt4->GetName()))->Fill(pvtx[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt4->GetName()))->Fill(dvtx[4]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt4->GetName()))->Fill(pType[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt5->GetName()))->Fill(pvtx[5]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt5->GetName()))->Fill(dvtx[5]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt5->GetName()))->Fill(pType[5]);

      ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_V[4]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_V[4]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_V[4]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_V[4]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_5->GetName()))->Fill( true_V[5]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_5->GetName()))->Fill( true_V[5]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_5->GetName()))->Fill(true_V[5]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_5->GetName()))->Fill( true_V[5]->E());

    }
    if(Npart==7){
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt4->GetName()))->Fill(pvtx[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt4->GetName()))->Fill(pvtx[4]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt4->GetName()))->Fill(pType[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt5->GetName()))->Fill(pvtx[5]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt5->GetName()))->Fill(pvtx[5]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt5->GetName()))->Fill(pType[5]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt6->GetName()))->Fill(pvtx[6]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt6->GetName()))->Fill(pvtx[6]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt6->GetName()))->Fill(pType[6]);

      ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_V[4]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_V[4]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_V[4]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_V[4]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_5->GetName()))->Fill( true_V[5]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_5->GetName()))->Fill( true_V[5]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_5->GetName()))->Fill(true_V[5]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_5->GetName()))->Fill( true_V[5]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_6->GetName()))->Fill( true_V[6]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_6->GetName()))->Fill( true_V[6]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_6->GetName()))->Fill(true_V[6]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_6->GetName()))->Fill( true_V[6]->E());

    }
    if(Npart==8){
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt4->GetName()))->Fill(pvtx[pType[4]]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt4->GetName()))->Fill(dvtx[pType[4]]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt4->GetName()))->Fill(pType[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt5->GetName()))->Fill(pvtx[pType[5]]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt5->GetName()))->Fill(dvtx[pType[5]]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt5->GetName()))->Fill(pType[5]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt6->GetName()))->Fill(pvtx[pType[6]]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt6->GetName()))->Fill( dvtx[pType[6]]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt6->GetName()))->Fill(pType[6]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt7->GetName()))->Fill(pvtx[pType[7]]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt7->GetName()))->Fill( dvtx[pType[7]]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt7->GetName()))->Fill(pType[7]);


      ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_V[4]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_V[4]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_V[4]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_V[4]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_5->GetName()))->Fill( true_V[5]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_5->GetName()))->Fill( true_V[5]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_5->GetName()))->Fill(true_V[5]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_5->GetName()))->Fill( true_V[5]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_6->GetName()))->Fill( true_V[6]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_6->GetName()))->Fill( true_V[6]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_6->GetName()))->Fill(true_V[6]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_6->GetName()))->Fill( true_V[6]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_7->GetName()))->Fill( true_V[7]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_7->GetName()))->Fill( true_V[7]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_7->GetName()))->Fill(true_V[7]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_7->GetName()))->Fill( true_V[7]->E());


     }
    //if(pType[7]!=-1)
    //  ((TH1F*)gDirectory->FindObject(h_mc_vtxdiff_pi_mu2->GetName()))->Fill(pvtx[pType[7]]);
  //}/



  }

  /// ************* tracks  &  clusters  *******************************

  // ~~~ CUT ~~~ CUT ~~~ CUT ~~~
  //if (flux_id ==1 ){
  //Ke2L3trig=sevt->DETstatus[0].LV3ABTrig & 1 ;
  //Km2L3trig=sevt->DETstatus[0].LV3ABTrig & 2 ;
  //if (Km2L3trig)
  //  cout << Km2L3trig << endl;
  if(nvtx != 1){cleanup();return 0;}
  if(IS_DATA){
    ((TH1F*)gDirectory->FindObject( h_vz_initial->GetName()))->Fill( sevt->vtx[0].z);
    if(sevt->vtx[0].z < -2000 || sevt->vtx[0].z > 8000) {cleanup();return 0;}
  }


  if(ntrack!=3 ) {cleanup();return 0;}

  //}
  //"3tracks"

  counter_ofmuons=0;

  for (int i=0; i<ntrack; i++) {
    //if(i==0)
    //  cout << "Loop INITIATING ! :)" << endl;
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

    //cout << "Loop start" << endl;
    //cout << "Event number = " << i << endl;
    //cout << "Imu_track = " << imu_track << endl;
    //if(imu_track > -1)
    //  counter_ofmuons++;
    //cout << "Icl_track = " <<icl_track << endl;
    //if(i==2)
    //  cout << "Number of muons for this event = " << counter_ofmuons << endl;
    //cout << "you are stupid !!!!!!!!!!!!!"  << IS_DATA << endl;


    if(icl_track>-1)
      clst_track     = sevt->cluster[icl_track].status;

    //        if(imu_track>-1)
    // cout<<i<<"    icl_track="<<icl_track<<endl;

    if (quality_track<0.8) { trkCorrSlopesV.push_back(NULL); trkCorrMidPointsV.push_back(NULL); continue; }
    //if (imu_track>-1) { trkCorrSlopesV.push_back(NULL); trkCorrMidPointsV.push_back(NULL); continue; }
    //if ( icl_track<0) { trkCorrSlopesV.push_back(NULL); trkCorrMidPointsV.push_back(NULL); continue; }
    if(IS_DATA)
      if (ddcell_track<2.) { trkCorrSlopesV.push_back(NULL); trkCorrMidPointsV.push_back(NULL); continue; }
    // if(flux_id==1){
    //  if (p_track<5.) { trkCorrSlopesV.push_back(NULL); trkCorrMidPointsV.push_back(NULL); continue; }
    //  if (p_track>50.) { trkCorrSlopesV.push_back(NULL); trkCorrMidPointsV.push_back(NULL); continue; }
    // }
    //if (flux_id==0){

    if (p_track<3.) { trkCorrSlopesV.push_back(NULL); trkCorrMidPointsV.push_back(NULL); continue; }
    if (p_track>50.) { trkCorrSlopesV.push_back(NULL); trkCorrMidPointsV.push_back(NULL); continue; }

    //if(IS_DATA)
      // if (clst_track>3) { trkCorrSlopesV.push_back(NULL); trkCorrMidPointsV.push_back(NULL); continue; }

    if (icl_track>-1)  {
      //calculate distance between the track and its matching cluster in LKr
      gclx=clcorrpos[icl_track]->x();
      gcly=clcorrpos[icl_track]->y();
      gclz=clcorrpos[icl_track]->z();
      trkxcl = sevt->track[i].x + (gclz-DCHz) * sevt->track[i].dxdz ;
      trkycl = sevt->track[i].y + (gclz-DCHz) * sevt->track[i].dydz ;
      dtrkcl = sqrt( pow(trkxcl-gclx, 2.) + pow(trkycl-gcly, 2.) );
    }
    ((TH1F*)gDirectory->FindObject(h_dtrkcl->GetName()))        ->Fill( dtrkcl );


    //if (IS_DATA)
    //  if(dtrkcl>6.)                 { trkCorrSlopesV.push_back(NULL); trkCorrMidPointsV.push_back(NULL); continue; }

    //number of good tracks , cut variable
    ngoodtrack++;


    //calculate Blue tube corrected slopes and track middle points before these corrections
    double trkSlopes[3]     = {0.};
    double trkMidPoints[3]     = {0.};

    applyBlueTubeCorrection(sevt, i, pKaon, vKaon, trkSlopes, trkMidPoints);

    trkCorrSlope    = new TVector3(trkSlopes[0],     trkSlopes[1],     trkSlopes[2]);
    trkCorrMidPoint = new TVector3(trkMidPoints[0], trkMidPoints[1],trkMidPoints[2]);
    trkCorrSlopesV.push_back(trkCorrSlope);
    trkCorrMidPointsV.push_back(trkCorrMidPoint);

    if( p_track > 10.  )    {
      if(pi1indx<0)
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
    if (IS_DATA){
      if(eop_track>c_eop_e  &&  eop_track<c_eop_e_up)    {
	if(eltrkindx1<0)
	  {
	    eltrkindx1=i;
	    eltrkq1 = sevt->track[i].q;
	  }
	else if(eltrkindx2<0 )    {
	  eltrkindx2=i;
	  eltrkq2 = sevt->track[i].q;
	}
	else { ; }
	nelectrack++;
      }
      else if( muchtrkindx < 0 && imu_track !=-1 && eop_track<c_eop_mu_up) {        muchtrkindx=i;
	muchtrkq=sevt->track[i].q;
	nmutrack++;

      }
      else    {
	;
      }
    }

    if (IS_MC){
      if(eop_track>c_eop_e  &&  eop_track<c_eop_e_up)    {
   	//if(eop_track> 0.8)    {
   	if(eltrkindx1<0)
   	  {
   	    eltrkindx1=i;
   	    eltrkq1 = sevt->track[i].q;
   	  }
   	else if(eltrkindx2<0 )    {
   	  eltrkindx2=i;
   	  eltrkq2 = sevt->track[i].q;
   	}
   	else { ; }
   	nelectrack++;
      }
      else if( muchtrkindx < 0) {
   	muchtrkindx=i;
   	muchtrkq=sevt->track[i].q;
   	nmutrack++;
      }
      else    {
   	;
      }
    }
    //Ke4
    //if (IS_MC){
    //  if( sevt->vtx[0].charge * sevt->track[i].q == 1)    {
    //	if(eop_track > 0.95){
    //	  eltrkindx1=i;
    //	  eltrkq1 = sevt->track[i].q;
    //	  nelectrack++;
    //	} else {
    //	  muchtrkindx=i;
    //	  muchtrkq=sevt->track[i].q;
    //	  nmutrack++;
    //	}
    //  }
    //  else if( eltrkindx2 < 0) {
    //	eltrkindx2=i;
    //	eltrkq2=sevt->track[i].q;
    //	nelectrack++;
    //  }
    //  else    {
    //	;
    //  }
    //}

  }




  ((TH1I*)gDirectory->FindObject(h_nrun->GetName()))        ->Fill( sbur->nrun );
  //CHANGE
  //if(IS_DATA){

  if(ngoodtrack!=3) {cleanup();return 0;}




  //cout << "Beginning D---------------------" << endl;
  //cout << Npart << "pT1= " << pType[1] << "pv1= "<< pvtx[1]<< "dv1= "<< dvtx[1]<<  endl;
  //cout << "END D---------------------" << endl;
  //pType.clear();
  //dvtx.clear();
  //pvtx.clear();
  //mu_pvtx.clear();
 // ~~~ CUT ~~~ CUT ~~~ CUT ~~~
  // cout << "WTFFFF" << nmutrack << nelectrack << endl;
  //Munue+e- track assignment selection
  //if(nelectrack+nmutrack==3)    {
  //  if(muchtrkq>0.)    {
  //    if(eltrkq1 != eltrkq2)    {// mu+, e+, e-
  //	evttypebytrks=0;
  //    }
  //    else if(eltrkq1>0)    {    // mu+, e+, e+
  //	evttypebytrks=1;
  //    }
  //    else    {                // mu+, e-, e-
  //	evttypebytrks=2;
  //    }
  //  }
  //  else    {
  //    if(eltrkq1 != eltrkq2)    {// mu-, e+, e-
  //	evttypebytrks=3;
  //    }
  //    else if(eltrkq1>0)    {    // mu-, e+, e+
  //	evttypebytrks=4;
  //    }
  //    else    {                // mu-, e-, e-
  //	evttypebytrks=5;
  //    }
  //  }
  //}
  if(nelectrack+nmutrack==3)    {
    if(nelectrack==2 && nmutrack==1 ){
      if((eltrkq1+eltrkq2+muchtrkq)==3.)    {// mu+, e+, e+
	evttypebytrks=1;
      } else if ((eltrkq1+eltrkq2+muchtrkq)== -3.){ // mu-, e-, e-
	evttypebytrks=5;
      } else if ((eltrkq1+eltrkq2+muchtrkq)== 1. && muchtrkq==1){// mu+, e+, e-
	evttypebytrks=0;
      } else if ((eltrkq1+eltrkq2+muchtrkq)== 1. && muchtrkq==-1.){// mu-, e+, e+
	evttypebytrks=4;
      } else if ((eltrkq1+eltrkq2+muchtrkq)== -1. && muchtrkq==-1.){// mu-, e+, e-
	evttypebytrks=3;
      } else if ((eltrkq1+eltrkq2+muchtrkq)== -1. && muchtrkq==1.){// mu+, e-, e-
	evttypebytrks=2;
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
  //  ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_track);
  ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_track, p_track);
  ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_track);
  //((TH1F*)gDirectory->FindObject(h_dtrkcl->GetName()))        ->Fill( dtrkcl );

  //cout << "LV2 trigger = " <<  (sevt->trigWord & 0 ) << (sevt->trigWord & 1 ) <<  (sevt->trigWord & 2 )<< (sevt->trigWord & 3 ) <<  (sevt->trigWord & 4 ) <<  (sevt->trigWord & 5 ) << endl;
  //cout << "Badb = " << sbur->BadB.Lkr << sbur->BadB.Dch << sbur->BadB.Skip<< sbur->BadB.Mbx<< sbur->BadB.Muv<< sbur->BadB.HodC<< sbur->BadB.Phys << endl;

  //Pion initial histograms

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


  pK1x = 0;
  pK1y = 0;
  pK1z = 60.;
  pKE =  sqrt( pK1z*pK1z + massKaonC*massKaonC);


  pK = new TLorentzVector(pK1x , pK1y , pK1z , pKE);

  //
  //
  //


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

    k3pic = new TLorentzVector( (*pi1lv) + (*pi2lv) + (*pi3lv) );



    k3pi_miss_mass = new TLorentzVector((*K4Momentum) - (*pi1lv) - (*pi2lv) - (*pi3lv) );
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

    //k3pi refferent
    //DCH INNER OUTER RADIUS CUTS FOR
    if (pow(pilp1_uncorr[0],2.) + pow(pilp1_uncorr[1],2.) < 12100. &&
	pow(pilp2_uncorr[0],2.) + pow(pilp2_uncorr[1],2.) < 12100. &&
	pow(pilp3_uncorr[0],2.) + pow(pilp3_uncorr[1],2.) < 12100. &&
	pow(pilp1_uncorr[0],2.) + pow(pilp1_uncorr[1],2.) > 196.   &&
	pow(pilp2_uncorr[0],2.) + pow(pilp2_uncorr[1],2.) > 196.   &&
	pow(pilp3_uncorr[0],2.) + pow(pilp3_uncorr[1],2.) > 196.   &&
	(*k3pic).M() > 0.4837                                      &&
	(*k3pic).M() < 0.5037
	){
      //  if(IS_DATA){
	if(  pi2trkhodtime-pi3trkhodtime  > -2 && pi2trkhodtime-pi3trkhodtime < 2 &&
	     pi2trkhodtime-pi1trkhodtime  > -2 && pi2trkhodtime - pi1trkhodtime < 2 &&
	     pi3trkhodtime-pi1trkhodtime  > -2 && pi3trkhodtime - pi1trkhodtime < 2 &&
	     fabs(pi2trktime - pi2trkhodtime)<12.          &&
	     fabs(pi3trktime - pi3trkhodtime)<12.          &&
	     fabs((long double)(pi2trktime-pi1trktime)) < 12. &&
	     fabs((long double)(pi3trktime-pi1trktime)) < 12. &&
	     fabs(pi1trktime-pi1trkhodtime) < 12.
	     ){

	  // }



	  if(((sevt->trigWord)>>3)&0x1) {
	    // ref. trigg. ok:
	    ((TH1I*)gDirectory->FindObject(h_k3pi_CPRE->GetName()))->Fill(k3pi_evttypebytrks);
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
		((TH1I*)gDirectory->FindObject(h_mu1->GetName()))->Fill(k3pi_evttypebytrks);
	      if( ((sevt->trigWord)>>1)&0x1  )
		((TH1I*)gDirectory->FindObject(h_Q1->GetName()))->Fill(k3pi_evttypebytrks);

	      if( ((sevt->trigWord))&0x1    || ((sevt->trigWord)>>1)&0x1 || ((sevt->trigWord)>>4)&0x1 )
		((TH1I*)gDirectory->FindObject(h_k3pi_full_trig->GetName()))->Fill(k3pi_evttypebytrks);
	    }
	  }






	  //epsilon a
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
	      (sevt->pu[6].chan[5]>>13)&0x1
	      ) {
	    if(((sevt->trigWord))&0x1   || ((sevt->trigWord)>>1)&0x1 || ((sevt->trigWord)>>4)&0x1 ) {
	      //My trigger also should be OK (constructed from PU)
	      ((TH1I*)gDirectory->FindObject(h_k3pi_2vtx->GetName()))->Fill(k3pi_evttypebytrks);
	      if(((sevt->trigWord)>>3)&0x1) {
		((TH1I*)gDirectory->FindObject(h_k3pi_1vtx_mcut->GetName()))->Fill(k3pi_evttypebytrks);
	      }
	    }
	    //tr_eff = hreal/href;
	  }



	  ((TH1I*)gDirectory->FindObject(h_k3pi_evttypebytrks->GetName()))->Fill(k3pi_evttypebytrks);
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

	  if(sevt->track[pi1indx].q * (sevt->track[pi1indx].q + sevt->track[pi2indx].q + sevt->track[pi3indx].q) == -1
	     && sevt->track[pi1indx].iMuon != -1){
	    ((TH2F*)gDirectory->FindObject(h_eop_vs_piodd->GetName()))->Fill(eop_pi1_track, (*pi1lv).P());
	    ((TH1F*)gDirectory->FindObject(h_piodd_eovp->GetName()))->Fill(eop_pi1_track);
	  } else if (sevt->track[pi2indx].q * (sevt->track[pi1indx].q + sevt->track[pi2indx].q + sevt->track[pi3indx].q) == -1
		     && sevt->track[pi2indx].iMuon != -1){
	    ((TH2F*)gDirectory->FindObject(h_eop_vs_piodd->GetName()))->Fill(eop_pi2_track, (*pi2lv).P());
	    ((TH1F*)gDirectory->FindObject(h_piodd_eovp->GetName()))->Fill( eop_pi2_track);
	  }  else if ( sevt->track[pi3indx].q * (sevt->track[pi1indx].q + sevt->track[pi2indx].q + sevt->track[pi3indx].q)== -1
		       && sevt->track[pi3indx].iMuon != -1){
	    ((TH2F*)gDirectory->FindObject(h_eop_vs_piodd->GetName()))->Fill(eop_pi3_track, (*pi3lv).P());
	    ((TH1F*)gDirectory->FindObject(h_piodd_eovp->GetName()))->Fill( eop_pi3_track);
	  }

	  ((TH1F*)gDirectory->FindObject(h_pi1trktime->GetName()))           ->Fill( pi1trktime );
	  ((TH1F*)gDirectory->FindObject(h_pi1trkhodtime->GetName()))        ->Fill( pi1trkhodtime );
	  ((TH1F*)gDirectory->FindObject(h_pi1trktimediff->GetName()))       ->Fill( pi1trktime - pi1trkhodtime );

	  ((TH1F*)gDirectory->FindObject(h_pi2trktime->GetName()))           ->Fill( pi2trktime );
	  ((TH1F*)gDirectory->FindObject(h_pi2trkhodtime->GetName()))        ->Fill( pi2trkhodtime );
	  ((TH1F*)gDirectory->FindObject(h_pi2trktimediff->GetName()))       ->Fill( pi2trktime - pi2trkhodtime );

	  ((TH1F*)gDirectory->FindObject(h_pi3trktime->GetName()))           ->Fill( pi3trktime );
	  ((TH1F*)gDirectory->FindObject(h_pi3trkhodtime->GetName()))        ->Fill( pi3trkhodtime );
	  ((TH1F*)gDirectory->FindObject(h_pi3trktimediff->GetName()))       ->Fill( pi3trktime - pi3trkhodtime );

	  //((TH1F*)gDirectory->FindObject(h_ecltimepi1->GetName()))           ->Fill( pi1cltime );
	  ((TH1F*)gDirectory->FindObject(h_pi12hodtimediff->GetName()))      ->Fill( pi1trkhodtime - pi2trkhodtime );
	  ((TH1F*)gDirectory->FindObject(h_pi13hodtimediff->GetName()))      ->Fill( pi1trkhodtime - pi3trkhodtime );
	  ((TH1F*)gDirectory->FindObject(h_pi23hodtimediff->GetName()))      ->Fill( pi2trkhodtime - pi3trkhodtime );
	  ((TH1F*)gDirectory->FindObject(h_KE->GetName()))->Fill((*K4Momentum).E());
	  ((TH1F*)gDirectory->FindObject(h_KP->GetName()))->Fill((*K4Momentum).P());

	}
    }//DCH inner outer radius and all other cuts
  }
  //  ------------------------------ END OF K3PIC selection part ------------------------------------

  if (eltrkindx1 == -1 || eltrkindx2 == -1 || muchtrkindx==-1)    {cleanup();return 0;}
  if (nelectrack!=2)        {cleanup();return 0;}
  //Making four vectors and time variables of the tracks
  //0000
  //cout << "sevt->vxt.charge = " << sevt->vtx[0].charge << ", charge from particles =" << (eltrkq1+eltrkq2+muchtrkq) << endl;
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

  el1E_pi = sqrt( p1x*p1x+p1y*p1y+p1z*p1z + massPionC*massPionC);
  el2E_pi = sqrt( p2x*p2x+p2y*p2y+p2z*p2z + massPionC*massPionC);


  el1lv_pi = new TLorentzVector(p1x, p1y, p1z, el1E_pi);
  el2lv_pi = new TLorentzVector(p2x, p2y, p2z, el2E_pi);
  eelv_pi  = new TLorentzVector( (*el1lv_pi)+(*el2lv_pi) );

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
  muchE_pi = sqrt( pmux*pmux+pmuy*pmuy+pmuz*pmuz + massPionC*massPionC);
  muchlv_pi = new TLorentzVector(pmux, pmuy, pmuz, muchE_pi);


  muchtrktime         = sevt->track[muchtrkindx].time;
  muchcltime         = sevt->cluster[ sevt->track[muchtrkindx].iClus ].time;
  muchtrkhodtime     = sevt->track[muchtrkindx].hodTime;




  muvee = new TLorentzVector( (*el1lv) + (*el2lv) + (*muchlv) );
  muvee_pi = new TLorentzVector( (*el1lv_pi) + (*el2lv_pi) + (*muchlv_pi) );
  muvee_piee = new TLorentzVector( (*el1lv) + (*el2lv) + (*muchlv_pi) );
  //pK1x = 0;
  //pK1y = 0;
  //pK1z = 60.;
  //pKE =  sqrt( pK1z*pK1z + massKaonC*massKaonC);
  //
  //
  //pK = new TLorentzVector(pK1x , pK1y , pK1z , pKE);


  neutrino = new TLorentzVector((*pK) - (*muchlv) - (*el1lv) - (*el2lv));
  mu_nu = new TLorentzVector( (*neutrino) + (*muchlv));
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


 //epsilon Q
   if(((sevt->trigWord)>>3)&0x1) {
    // ref. trigg. ok:
    ((TH1I*)gDirectory->FindObject(h_CPRE->GetName()))->Fill(evttypebytrks);
    if( (sevt->pu[3].chan[5]>>4)&0x1 ||
	(sevt->pu[4].chan[5]>>4)&0x1 ||//2vtx
	(sevt->pu[5].chan[5]>>4)&0x1 ||  //5 is the central time slice
	(sevt->pu[6].chan[5]>>4)&0x1 ||
	(sevt->pu[3].chan[5]>>8)&0x1 ||//1vtx
	(sevt->pu[4].chan[5]>>8)&0x1 ||
	(sevt->pu[5].chan[5]>>8)&0x1 ||  //5 is the central time slice
	(sevt->pu[6].chan[5]>>8)&0x1 ||
	(sevt->pu[3].chan[5]>>13)&0x1 ||
	(sevt->pu[4].chan[5]>>13)&0x1 ||
	(sevt->pu[5].chan[5]>>13)&0x1 ||  //5 is the central time slice
	(sevt->pu[6].chan[5]>>13)&0x1) {
      if(((sevt->trigWord))&0x1   || ((sevt->trigWord)>>2)&0x1 || ((sevt->trigWord)>>5)&0x1 ) {
	//My trigger also should be OK (constructed from PU)
	((TH1I*)gDirectory->FindObject(h_full_trig->GetName()))->Fill(evttypebytrks);
      }
      //tr_eff = hreal/href;
    }
   }

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

     //Trk hodtimediff
     ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( el1trkhodtime-el2trkhodtime );
     ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( el1trkhodtime-muchtrkhodtime );
     ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( el2trkhodtime-muchtrkhodtime );

     //epsilon a
     //if( (sevt->pu[3].chan[5]>>4)&0x1 ||
     //    (sevt->pu[4].chan[5]>>4)&0x1 ||
     //    (sevt->pu[5].chan[5]>>4)&0x1 ||  //5 is the central time slice
     //    (sevt->pu[6].chan[5]>>4)&0x1 ||
     //    (sevt->pu[3].chan[5]>>8)&0x1 ||
     //    (sevt->pu[4].chan[5]>>8)&0x1 ||
     //    (sevt->pu[5].chan[5]>>8)&0x1 ||  //5 is the central time slice
     //    (sevt->pu[6].chan[5]>>8)&0x1 ||
     //    (sevt->pu[3].chan[5]>>13)&0x1 ||
     //    (sevt->pu[4].chan[5]>>13)&0x1 ||
     //    (sevt->pu[5].chan[5]>>13)&0x1 ||  //5 is the central time slice
     //    (sevt->pu[6].chan[5]>>13)&0x1
     //    ) {
     //tr_eff = hreal/href;
     // }



  //v_second_e1 = new TLorentzVector;



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

  ((TH1F*)gDirectory->FindObject(h_gmueeE->GetName()))->Fill( (*mu_nu).E());
  ((TH1F*)gDirectory->FindObject(h_gmueeP->GetName()))->Fill( (*mu_nu).P());
  ((TH1F*)gDirectory->FindObject(h_gmueePt->GetName()))->Fill((*mu_nu).Pt());
  ((TH1F*)gDirectory->FindObject(h_gmueeM->GetName()))->Fill ((*mu_nu).M());

  ((TH1F*)gDirectory->FindObject(h_bx_much->GetName()))            ->Fill(muchp_uncorr[0]);
  ((TH1F*)gDirectory->FindObject(h_by_much->GetName()))            ->Fill(muchp_uncorr[1]);
  ((TH1F*)gDirectory->FindObject(h_bdxdz_much->GetName()))            ->Fill(muchv_uncorr[0]);
  ((TH1F*)gDirectory->FindObject(h_bdydz_much->GetName()))            ->Fill(muchv_uncorr[1]) ;
  ((TH1F*)gDirectory->FindObject(h_bx_vs_by_much->GetName()))        ->Fill(muchp_uncorr[0], muchp_uncorr[1]);
  ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_much->GetName()))->Fill(muchv_uncorr[0], muchv_uncorr[1]);

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
  if(nclust< 2  ) {cleanup();return 0;}
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
    ////if(true_4.P()  < 3. || true_5.P() < 3.) {cleanup();return 0;}
    ((TH1F*)gDirectory->FindObject(h_mc_Npart_tr->GetName()))->Fill(Npart);
    if (Npart >=4 ){

      ((TH1D*)gDirectory->FindObject(h_mc_muvee_M->GetName()))->Fill( Muee_three_track.M());
      ((TH1D*)gDirectory->FindObject(h_mc_muvee_P->GetName()))->Fill( Muee_three_track.P());
      ((TH1D*)gDirectory->FindObject(h_mc_muvee_Pt->GetName()))->Fill( Muee_three_track.Pt());
      ((TH1D*)gDirectory->FindObject(h_mc_elsP->GetName()))->Fill( mc_ee.P());
      ((TH1D*)gDirectory->FindObject(h_mc_elsM->GetName()))->Fill( mc_ee.M());
      ((TH1D*)gDirectory->FindObject(h_mc_nu_M->GetName()))->Fill( Nu_true.M());
      ((TH1D*)gDirectory->FindObject(h_mc_nu_M2->GetName()))->Fill( Nu_true.M2());
      ((TH1D*)gDirectory->FindObject(h_mc_nu_P->GetName()))->Fill( Nu_true.P());
      ((TH1D*)gDirectory->FindObject(h_mc_nu_Pt->GetName()))->Fill(Nu_true.Pt());
      ((TH1D*)gDirectory->FindObject(h_true_M_0->GetName()))->Fill( true_V[0]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_0->GetName()))->Fill( true_V[0]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_0->GetName()))->Fill(true_V[0]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_0->GetName()))->Fill( true_V[0]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_1->GetName()))->Fill( true_V[1]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_1->GetName()))->Fill( true_V[1]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_1->GetName()))->Fill(true_V[1]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_1->GetName()))->Fill( true_V[1]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_2->GetName()))->Fill( true_V[3]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_2->GetName()))->Fill( true_V[3]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_2->GetName()))->Fill(true_V[3]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_2->GetName()))->Fill( true_V[3]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_3->GetName()))->Fill( true_V[3]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_3->GetName()))->Fill( true_V[3]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_3->GetName()))->Fill(true_V[3]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_3->GetName()))->Fill( true_V[3]->E());

      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt0->GetName()))->Fill(pvtx[0]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt1->GetName()))->Fill(pvtx[1]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt2->GetName()))->Fill(pvtx[2]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt3->GetName()))->Fill(pvtx[3]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt0->GetName()))->Fill(dvtx[0]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt1->GetName()))->Fill(dvtx[1]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt2->GetName()))->Fill(dvtx[2]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt3->GetName()))->Fill(dvtx[3]);
       }
    if(Npart==5){
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt4->GetName()))->Fill(pvtx[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt4->GetName()))->Fill(dvtx[4]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt4->GetName()))->Fill(pType[4]);
      ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_V[4]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_V[4]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_V[4]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_V[4]->E());

    }
    if(Npart==6){
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt4->GetName()))->Fill(pvtx[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt4->GetName()))->Fill(dvtx[4]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt4->GetName()))->Fill(pType[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt5->GetName()))->Fill(pvtx[5]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt5->GetName()))->Fill(dvtx[5]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt5->GetName()))->Fill(pType[5]);

      ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_V[4]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_V[4]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_V[4]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_V[4]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_5->GetName()))->Fill( true_V[5]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_5->GetName()))->Fill( true_V[5]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_5->GetName()))->Fill(true_V[5]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_5->GetName()))->Fill( true_V[5]->E());

    }
    if(Npart==7){
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt4->GetName()))->Fill(pvtx[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt4->GetName()))->Fill(pvtx[4]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt4->GetName()))->Fill(pType[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt5->GetName()))->Fill(pvtx[5]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt5->GetName()))->Fill(pvtx[5]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt5->GetName()))->Fill(pType[5]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt6->GetName()))->Fill(pvtx[6]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt6->GetName()))->Fill(pvtx[6]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt6->GetName()))->Fill(pType[6]);

      ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_V[4]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_V[4]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_V[4]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_V[4]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_5->GetName()))->Fill( true_V[5]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_5->GetName()))->Fill( true_V[5]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_5->GetName()))->Fill(true_V[5]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_5->GetName()))->Fill( true_V[5]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_6->GetName()))->Fill( true_V[6]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_6->GetName()))->Fill( true_V[6]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_6->GetName()))->Fill(true_V[6]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_6->GetName()))->Fill( true_V[6]->E());

    }
    if(Npart==8){
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt4->GetName()))->Fill(pvtx[pType[4]]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt4->GetName()))->Fill(dvtx[pType[4]]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt4->GetName()))->Fill(pType[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt5->GetName()))->Fill(pvtx[pType[5]]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt5->GetName()))->Fill(dvtx[pType[5]]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt5->GetName()))->Fill(pType[5]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt6->GetName()))->Fill(pvtx[pType[6]]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt6->GetName()))->Fill( dvtx[pType[6]]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt6->GetName()))->Fill(pType[6]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt7->GetName()))->Fill(pvtx[pType[7]]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt7->GetName()))->Fill( dvtx[pType[7]]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt7->GetName()))->Fill(pType[7]);


      ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_V[4]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_V[4]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_V[4]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_V[4]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_5->GetName()))->Fill( true_V[5]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_5->GetName()))->Fill( true_V[5]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_5->GetName()))->Fill(true_V[5]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_5->GetName()))->Fill( true_V[5]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_6->GetName()))->Fill( true_V[6]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_6->GetName()))->Fill( true_V[6]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_6->GetName()))->Fill(true_V[6]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_6->GetName()))->Fill( true_V[6]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_7->GetName()))->Fill( true_V[7]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_7->GetName()))->Fill( true_V[7]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_7->GetName()))->Fill(true_V[7]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_7->GetName()))->Fill( true_V[7]->E());


     }
    //if(pType[7]!=-1)
    //  ((TH1F*)gDirectory->FindObject(h_mc_vtxdiff_pi_mu2->GetName()))->Fill(pvtx[pType[7]]);
    }



//epsilon Q
   if(((sevt->trigWord)>>3)&0x1) {
    // ref. trigg. ok:
    ((TH1I*)gDirectory->FindObject(h_CPRE->GetName()))->Fill(evttypebytrks);
    if( (sevt->pu[3].chan[5]>>4)&0x1 ||
	(sevt->pu[4].chan[5]>>4)&0x1 ||//2vtx
	(sevt->pu[5].chan[5]>>4)&0x1 ||  //5 is the central time slice
	(sevt->pu[6].chan[5]>>4)&0x1 ||
	(sevt->pu[3].chan[5]>>8)&0x1 ||//1vtx
	(sevt->pu[4].chan[5]>>8)&0x1 ||
	(sevt->pu[5].chan[5]>>8)&0x1 ||  //5 is the central time slice
	(sevt->pu[6].chan[5]>>8)&0x1 ||
	(sevt->pu[3].chan[5]>>13)&0x1 ||
	(sevt->pu[4].chan[5]>>13)&0x1 ||
	(sevt->pu[5].chan[5]>>13)&0x1 ||  //5 is the central time slice
	(sevt->pu[6].chan[5]>>13)&0x1) {
      if(((sevt->trigWord))&0x1   || ((sevt->trigWord)>>2)&0x1 || ((sevt->trigWord)>>5)&0x1 ) {
	//My trigger also should be OK (constructed from PU)
	((TH1I*)gDirectory->FindObject(h_full_trig->GetName()))->Fill(evttypebytrks);

      }
      //tr_eff = hreal/href;
    }
  }

  //epsilon a
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
      (sevt->pu[6].chan[5]>>13)&0x1
      ) {
    if(((sevt->trigWord))&0x1   || ((sevt->trigWord)>>2)&0x1 || ((sevt->trigWord)>>5)&0x1 ) {
      //My trigger also should be OK (constructed from PU)
      ((TH1I*)gDirectory->FindObject(h_mu1->GetName()))->Fill(evttypebytrks);
      if(((sevt->trigWord)>>3)&0x1) {
	((TH1I*)gDirectory->FindObject(h_Q1->GetName()))->Fill(evttypebytrks);
      }
    }
    //tr_eff = hreal/href;
  }




  ((TH1F*)gDirectory->FindObject( h_vz_initial->GetName()))->Fill( sevt->vtx[0].z);
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
  //((TH1F*)gDirectory->FindObject(h_e1e2timediff->GetName()))        ->Fill( (*muvee_pi).E() );
  ((TH1F*)gDirectory->FindObject(h_e1e2cltimediff->GetName()))    ->Fill( el1cltime-el2cltime );
  ((TH1F*)gDirectory->FindObject(h_e1e2hodtimediff->GetName()))        ->Fill( el1trkhodtime-el2trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_e1muhodtimediff->GetName()))        ->Fill( el1trkhodtime-muchtrkhodtime );
  ((TH1F*)gDirectory->FindObject(h_e2muhodtimediff->GetName()))        ->Fill( el2trkhodtime-muchtrkhodtime );
  //Trk hodtimediff
  ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( el1trkhodtime-el2trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( el1trkhodtime-muchtrkhodtime );
  ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( el2trkhodtime-muchtrkhodtime );

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
   ((TH1F*)gDirectory->FindObject(h_gmueeE->GetName()))->Fill( (*mu_nu).E());
  ((TH1F*)gDirectory->FindObject(h_gmueeP->GetName()))->Fill( (*mu_nu).P());
  ((TH1F*)gDirectory->FindObject(h_gmueePt->GetName()))->Fill((*mu_nu).Pt());
  ((TH1F*)gDirectory->FindObject(h_gmueeM->GetName()))->Fill ((*mu_nu).M());

  ((TH1F*)gDirectory->FindObject(h_bx_much->GetName()))            ->Fill(muchp_uncorr[0]);
  ((TH1F*)gDirectory->FindObject(h_by_much->GetName()))            ->Fill(muchp_uncorr[1]);
  ((TH1F*)gDirectory->FindObject(h_bdxdz_much->GetName()))            ->Fill(muchv_uncorr[0]);
  ((TH1F*)gDirectory->FindObject(h_bdydz_much->GetName()))            ->Fill(muchv_uncorr[1]) ;
  ((TH1F*)gDirectory->FindObject(h_bx_vs_by_much->GetName()))        ->Fill(muchp_uncorr[0], muchp_uncorr[1]);
  ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_much->GetName()))->Fill(muchv_uncorr[0], muchv_uncorr[1]);


  ((TH1F*)gDirectory->FindObject(h_muche1distDHC1plane->GetName()))->Fill( muche1distDHC1plane );
  ((TH1F*)gDirectory->FindObject(h_muche2distDHC1plane->GetName()))->Fill( muche2distDHC1plane );
  ((TH1F*)gDirectory->FindObject(h_muche1distLKrplane->GetName()))->Fill(  muche1distLKrplane );
  ((TH1F*)gDirectory->FindObject(h_muche2distLKrplane->GetName()))->Fill(  muche2distLKrplane );

  ((TH1F*)gDirectory->FindObject(h_eldistDCH1plane->GetName()))->Fill( eldistDCH1plane );
  ((TH1F*)gDirectory->FindObject(h_eldistLKrplane->GetName()))->Fill( eldistLKrplane );
  //cout << (*muvee_pi).E() << (*muvee).E() << endl;

  ((TH1F*)gDirectory->FindObject(h_k3pi_nu_m->GetName()))->Fill( (*muvee_pi).M() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_nu_P->GetName()))->Fill( (*muvee_pi).P() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_nu_E->GetName()))->Fill( (*muvee_pi).E() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_nu_Pt->GetName()))->Fill( (*muvee_pi).Pt() );
 ((TH1F*)gDirectory->FindObject(h_muee_piee_M->GetName()))->Fill( (*muvee_piee).M() );
  ((TH1F*)gDirectory->FindObject(h_muee_piee_P->GetName()))->Fill( (*muvee_piee).P() );
  ((TH1F*)gDirectory->FindObject(h_muee_piee_E->GetName()))->Fill( (*muvee_piee).E() );
  ((TH1F*)gDirectory->FindObject(h_muee_piee_Pt->GetName()))->Fill( (*muvee_piee).Pt() );

  //((TH1F*)gDirectory->FindObject(h_3trk_invm->GetName()))->Fill ((*k3pic).M2());
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_P->GetName()))->Fill( (*muchlv_pi).P() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_P->GetName()))->Fill( (*el1lv_pi).P() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_P->GetName()))->Fill( (*el2lv_pi).P() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_E->GetName()))->Fill( (*muchlv_pi).E() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_E->GetName()))->Fill( (*el1lv_pi).E() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_E->GetName()))->Fill( (*el2lv_pi).E() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_Pt->GetName()))->Fill((*muchlv_pi).Pt() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_Pt->GetName()))->Fill((*el1lv_pi).Pt() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_Pt->GetName()))->Fill((*el2lv_pi).Pt() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_m->GetName()))->Fill( (*muchlv_pi).M() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_m->GetName()))->Fill( (*el1lv_pi).M() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_m->GetName()))->Fill( (*el2lv_pi).M() );


  file1->cd();
  h_cutflow->Fill(cutcounter);
  h_cutflow->GetXaxis()->SetBinLabel( cutcounter, cuts[cutcounter-1].c_str() );
  cutcounter++;
  //THE END OF ************* Geometry cuts ****************** ///




  /// ************* Momentum cuts ****************** ///
  //if(flux_id==1){
  //  ;
  //}
  //if (flux_id==0){
  if ( (*el1lv).P() < 3. || (*el2lv).P() < 3. || (*muchlv).P() < 10.) {cleanup();return 0;}
  //}
  if ((*muvee).P() > 66. ) {cleanup();return 0;}
  if(sevt->vtx[0].chi2 > 20. ) {cleanup();return 0;}
  dirs[cutcounter-1]->cd();

  if (IS_MC){
    ////if(true_4.P()  < 3. || true_5.P() < 3.) {cleanup();return 0;}
    ((TH1F*)gDirectory->FindObject(h_mc_Npart_tr->GetName()))->Fill(Npart);
    if (Npart >=4 ){

      ((TH1D*)gDirectory->FindObject( h_mc_muvee_M->GetName()))->Fill( Muee_three_track.M());
      ((TH1D*)gDirectory->FindObject( h_mc_muvee_P->GetName()))->Fill( Muee_three_track.P());
      ((TH1D*)gDirectory->FindObject( h_mc_muvee_Pt->GetName()))->Fill( Muee_three_track.Pt());
      ((TH1D*)gDirectory->FindObject( h_mc_elsP->GetName()))->Fill( mc_ee.P());
      ((TH1D*)gDirectory->FindObject( h_mc_elsM->GetName()))->Fill( mc_ee.M());
      ((TH1D*)gDirectory->FindObject( h_mc_nu_M->GetName()))->Fill( Nu_true.M());
      ((TH1D*)gDirectory->FindObject( h_mc_nu_M2->GetName()))->Fill( Nu_true.M2());
      ((TH1D*)gDirectory->FindObject( h_mc_nu_P->GetName()))->Fill( Nu_true.P());
      ((TH1D*)gDirectory->FindObject( h_mc_nu_Pt->GetName()))->Fill(Nu_true.Pt());
      ((TH1D*)gDirectory->FindObject(h_true_M_0->GetName()))->Fill( true_V[0]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_0->GetName()))->Fill( true_V[0]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_0->GetName()))->Fill(true_V[0]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_0->GetName()))->Fill( true_V[0]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_1->GetName()))->Fill( true_V[1]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_1->GetName()))->Fill( true_V[1]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_1->GetName()))->Fill(true_V[1]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_1->GetName()))->Fill( true_V[1]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_2->GetName()))->Fill( true_V[3]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_2->GetName()))->Fill( true_V[3]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_2->GetName()))->Fill(true_V[3]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_2->GetName()))->Fill( true_V[3]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_3->GetName()))->Fill( true_V[3]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_3->GetName()))->Fill( true_V[3]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_3->GetName()))->Fill(true_V[3]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_3->GetName()))->Fill( true_V[3]->E());

      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt0->GetName()))->Fill(pvtx[0]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt1->GetName()))->Fill(pvtx[1]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt2->GetName()))->Fill(pvtx[2]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt3->GetName()))->Fill(pvtx[3]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt0->GetName()))->Fill(dvtx[0]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt1->GetName()))->Fill(dvtx[1]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt2->GetName()))->Fill(dvtx[2]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt3->GetName()))->Fill(dvtx[3]);
    }
    if(Npart==5){
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt4->GetName()))->Fill(pvtx[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt4->GetName()))->Fill(dvtx[4]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt4->GetName()))->Fill(pType[4]);
      ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_V[4]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_V[4]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_V[4]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_V[4]->E());

    }
    if(Npart==6){
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt4->GetName()))->Fill(pvtx[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt4->GetName()))->Fill(dvtx[4]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt4->GetName()))->Fill(pType[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt5->GetName()))->Fill(pvtx[5]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt5->GetName()))->Fill(dvtx[5]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt5->GetName()))->Fill(pType[5]);

      ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_V[4]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_V[4]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_V[4]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_V[4]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_5->GetName()))->Fill( true_V[5]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_5->GetName()))->Fill( true_V[5]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_5->GetName()))->Fill(true_V[5]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_5->GetName()))->Fill( true_V[5]->E());

    }
    if(Npart==7){
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt4->GetName()))->Fill(pvtx[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt4->GetName()))->Fill(pvtx[4]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt4->GetName()))->Fill(pType[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt5->GetName()))->Fill(pvtx[5]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt5->GetName()))->Fill(pvtx[5]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt5->GetName()))->Fill(pType[5]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt6->GetName()))->Fill(pvtx[6]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt6->GetName()))->Fill(pvtx[6]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt6->GetName()))->Fill(pType[6]);

      ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_V[4]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_V[4]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_V[4]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_V[4]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_5->GetName()))->Fill( true_V[5]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_5->GetName()))->Fill( true_V[5]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_5->GetName()))->Fill(true_V[5]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_5->GetName()))->Fill( true_V[5]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_6->GetName()))->Fill( true_V[6]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_6->GetName()))->Fill( true_V[6]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_6->GetName()))->Fill(true_V[6]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_6->GetName()))->Fill( true_V[6]->E());

    }
    if(Npart==8){
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt4->GetName()))->Fill(pvtx[pType[4]]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt4->GetName()))->Fill(dvtx[pType[4]]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt4->GetName()))->Fill(pType[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt5->GetName()))->Fill(pvtx[pType[5]]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt5->GetName()))->Fill(dvtx[pType[5]]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt5->GetName()))->Fill(pType[5]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt6->GetName()))->Fill(pvtx[pType[6]]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt6->GetName()))->Fill( dvtx[pType[6]]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt6->GetName()))->Fill(pType[6]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt7->GetName()))->Fill(pvtx[pType[7]]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt7->GetName()))->Fill( dvtx[pType[7]]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt7->GetName()))->Fill(pType[7]);


      ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_V[4]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_V[4]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_V[4]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_V[4]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_5->GetName()))->Fill( true_V[5]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_5->GetName()))->Fill( true_V[5]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_5->GetName()))->Fill(true_V[5]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_5->GetName()))->Fill( true_V[5]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_6->GetName()))->Fill( true_V[6]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_6->GetName()))->Fill( true_V[6]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_6->GetName()))->Fill(true_V[6]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_6->GetName()))->Fill( true_V[6]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_7->GetName()))->Fill( true_V[7]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_7->GetName()))->Fill( true_V[7]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_7->GetName()))->Fill(true_V[7]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_7->GetName()))->Fill( true_V[7]->E());


    }
    //if(pType[7]!=-1)
    //  ((TH1F*)gDirectory->FindObject(h_mc_vtxdiff_pi_mu2->GetName()))->Fill(pvtx[pType[7]]);
  }

//epsilon Q
   if(((sevt->trigWord)>>3)&0x1) {
    // ref. trigg. ok:
    ((TH1I*)gDirectory->FindObject(h_CPRE->GetName()))->Fill(evttypebytrks);
    if( (sevt->pu[3].chan[5]>>4)&0x1 ||
	(sevt->pu[4].chan[5]>>4)&0x1 ||//2vtx
	(sevt->pu[5].chan[5]>>4)&0x1 ||  //5 is the central time slice
	(sevt->pu[6].chan[5]>>4)&0x1 ||
	(sevt->pu[3].chan[5]>>8)&0x1 ||//1vtx
	(sevt->pu[4].chan[5]>>8)&0x1 ||
	(sevt->pu[5].chan[5]>>8)&0x1 ||  //5 is the central time slice
	(sevt->pu[6].chan[5]>>8)&0x1 ||
	(sevt->pu[3].chan[5]>>13)&0x1 ||
	(sevt->pu[4].chan[5]>>13)&0x1 ||
	(sevt->pu[5].chan[5]>>13)&0x1 ||  //5 is the central time slice
	(sevt->pu[6].chan[5]>>13)&0x1) {
      if(((sevt->trigWord))&0x1   || ((sevt->trigWord)>>2)&0x1 || ((sevt->trigWord)>>5)&0x1 ) {
	//My trigger also should be OK (constructed from PU)
	((TH1I*)gDirectory->FindObject(h_full_trig->GetName()))->Fill(evttypebytrks);

      }
      //tr_eff = hreal/href;
    }
  }

  //epsilon a
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
      (sevt->pu[6].chan[5]>>13)&0x1
      ) {
    if(((sevt->trigWord))&0x1   || ((sevt->trigWord)>>2)&0x1 || ((sevt->trigWord)>>5)&0x1 ) {
      //My trigger also should be OK (constructed from PU)
      ((TH1I*)gDirectory->FindObject(h_mu1->GetName()))->Fill(evttypebytrks);
      if(((sevt->trigWord)>>3)&0x1) {
	((TH1I*)gDirectory->FindObject(h_Q1->GetName()))->Fill(evttypebytrks);
      }
    }
    //tr_eff = hreal/href;
  }



  ((TH1F*)gDirectory->FindObject( h_vz_initial->GetName()))->Fill( sevt->vtx[0].z);
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

   //Trk hodtimediff
  ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( el1trkhodtime-el2trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( el1trkhodtime-muchtrkhodtime );
  ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( el2trkhodtime-muchtrkhodtime );


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

   ((TH1F*)gDirectory->FindObject(h_gmueeE->GetName()))->Fill( (*mu_nu).E());
  ((TH1F*)gDirectory->FindObject(h_gmueeP->GetName()))->Fill( (*mu_nu).P());
  ((TH1F*)gDirectory->FindObject(h_gmueePt->GetName()))->Fill((*mu_nu).Pt());
  ((TH1F*)gDirectory->FindObject(h_gmueeM->GetName()))->Fill ((*mu_nu).M());
  ((TH1F*)gDirectory->FindObject(h_bx_much->GetName()))            ->Fill(muchp_uncorr[0]);
  ((TH1F*)gDirectory->FindObject(h_by_much->GetName()))            ->Fill(muchp_uncorr[1]);
  ((TH1F*)gDirectory->FindObject(h_bdxdz_much->GetName()))            ->Fill(muchv_uncorr[0]);
  ((TH1F*)gDirectory->FindObject(h_bdydz_much->GetName()))            ->Fill(muchv_uncorr[1]) ;
  ((TH1F*)gDirectory->FindObject(h_bx_vs_by_much->GetName()))        ->Fill(muchp_uncorr[0], muchp_uncorr[1]);
  ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_much->GetName()))->Fill(muchv_uncorr[0], muchv_uncorr[1]);


  ((TH1F*)gDirectory->FindObject(h_k3pi_nu_m->GetName()))->Fill( (*muvee_pi.M() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_nu_P->GetName()))->Fill( (*muvee_pi).P() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_nu_E->GetName()))->Fill( (*muvee_pi).E() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_nu_Pt->GetName()))->Fill( (*muvee_pi).Pt() );
   ((TH1F*)gDirectory->FindObject(h_muee_piee_M->GetName()))->Fill( (*muvee_piee).M() );
  ((TH1F*)gDirectory->FindObject(h_muee_piee_P->GetName()))->Fill( (*muvee_piee).P() );
  ((TH1F*)gDirectory->FindObject(h_muee_piee_E->GetName()))->Fill( (*muvee_piee).E() );
  ((TH1F*)gDirectory->FindObject(h_muee_piee_Pt->GetName()))->Fill( (*muvee_piee).Pt() );

  //((TH1F*)gDirectory->FindObject(h_3trk_invm->GetName()))->Fill ((*k3pic).M2());
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_P->GetName()))->Fill( (*muchlv_pi).P() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_P->GetName()))->Fill( (*el1lv_pi).P() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_P->GetName()))->Fill( (*el2lv_pi).P() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_E->GetName()))->Fill( (*muchlv_pi).E() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_E->GetName()))->Fill( (*el1lv_pi).E() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_E->GetName()))->Fill( (*el2lv_pi).E() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_Pt->GetName()))->Fill((*muchlv_pi).Pt() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_Pt->GetName()))->Fill((*el1lv_pi).Pt() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_Pt->GetName()))->Fill((*el2lv_pi).Pt() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_m->GetName()))->Fill( (*muchlv_pi).M() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_m->GetName()))->Fill( (*el1lv_pi).M() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_m->GetName()))->Fill( (*el2lv_pi).M() );




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
	//muchtrkq==1         &&
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

       ((TH1F*)gDirectory->FindObject( h_vz_initial->GetName()))->Fill( sevt->vtx[0].z);
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

      //Trk hodtimediff
      ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( el1trkhodtime-el2trkhodtime );
      ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( el1trkhodtime-muchtrkhodtime );
      ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( el2trkhodtime-muchtrkhodtime );

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

       ((TH1F*)gDirectory->FindObject(h_gmueeE->GetName()))->Fill( (*mu_nu).E());
  ((TH1F*)gDirectory->FindObject(h_gmueeP->GetName()))->Fill( (*mu_nu).P());
  ((TH1F*)gDirectory->FindObject(h_gmueePt->GetName()))->Fill((*mu_nu).Pt());
  ((TH1F*)gDirectory->FindObject(h_gmueeM->GetName()))->Fill ((*mu_nu).M());
      ((TH1F*)gDirectory->FindObject(h_bx_much->GetName()))            ->Fill(muchp_uncorr[0]);
      ((TH1F*)gDirectory->FindObject(h_by_much->GetName()))            ->Fill(muchp_uncorr[1]);
      ((TH1F*)gDirectory->FindObject(h_bdxdz_much->GetName()))            ->Fill(muchv_uncorr[0]);
      ((TH1F*)gDirectory->FindObject(h_bdydz_much->GetName()))            ->Fill(muchv_uncorr[1]) ;
      ((TH1F*)gDirectory->FindObject(h_bx_vs_by_much->GetName()))        ->Fill(muchp_uncorr[0], muchp_uncorr[1]);
      ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_much->GetName()))->Fill(muchv_uncorr[0], muchv_uncorr[1]);


      ((TH1F*)gDirectory->FindObject(h_k3pi_nu_m->GetName()))->Fill( (*muvee_pi).M() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_nu_P->GetName()))->Fill( (*muvee_pi).P() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_nu_E->GetName()))->Fill( (*muvee_pi).E() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_nu_Pt->GetName()))->Fill( (*muvee_pi).Pt() );
      ((TH1F*)gDirectory->FindObject(h_muee_piee_M->GetName()))->Fill( (*muvee_piee).M() );
      ((TH1F*)gDirectory->FindObject(h_muee_piee_P->GetName()))->Fill( (*muvee_piee).P() );
      ((TH1F*)gDirectory->FindObject(h_muee_piee_E->GetName()))->Fill( (*muvee_piee).E() );
      ((TH1F*)gDirectory->FindObject(h_muee_piee_Pt->GetName()))->Fill( (*muvee_piee).Pt() );

      //((TH1F*)gDirectory->FindObject(h_3trk_invm->GetName()))->Fill ((*k3pic).M2());
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_P->GetName()))->Fill( (*muchlv_pi).P() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_P->GetName()))->Fill( (*el1lv_pi).P() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_P->GetName()))->Fill( (*el2lv_pi).P() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_E->GetName()))->Fill( (*muchlv_pi).E() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_E->GetName()))->Fill( (*el1lv_pi).E() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_E->GetName()))->Fill( (*el2lv_pi).E() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_Pt->GetName()))->Fill((*muchlv_pi).Pt() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_Pt->GetName()))->Fill((*el1lv_pi).Pt() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_Pt->GetName()))->Fill((*el2lv_pi).Pt() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_m->GetName()))->Fill( (*muchlv_pi).M() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_m->GetName()))->Fill( (*el1lv_pi).M() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_m->GetName()))->Fill( (*el2lv_pi).M() );



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
    ////if(true_4.P()  < 3. || true_5.P() < 3.) {cleanup();return 0;}
    ((TH1F*)gDirectory->FindObject(h_mc_Npart_tr->GetName()))->Fill(Npart);
    if (Npart >=4 ){

      ((TH1D*)gDirectory->FindObject( h_mc_muvee_M->GetName()))->Fill( Muee_three_track.M());
      ((TH1D*)gDirectory->FindObject( h_mc_muvee_P->GetName()))->Fill( Muee_three_track.P());
      ((TH1D*)gDirectory->FindObject( h_mc_muvee_Pt->GetName()))->Fill( Muee_three_track.Pt());
      ((TH1D*)gDirectory->FindObject( h_mc_elsP->GetName()))->Fill( mc_ee.P());
      ((TH1D*)gDirectory->FindObject( h_mc_elsM->GetName()))->Fill( mc_ee.M());
      ((TH1D*)gDirectory->FindObject( h_mc_nu_M->GetName()))->Fill( Nu_true.M());
      ((TH1D*)gDirectory->FindObject( h_mc_nu_M2->GetName()))->Fill( Nu_true.M2());
      ((TH1D*)gDirectory->FindObject( h_mc_nu_P->GetName()))->Fill( Nu_true.P());
      ((TH1D*)gDirectory->FindObject( h_mc_nu_Pt->GetName()))->Fill(Nu_true.Pt());
      ((TH1D*)gDirectory->FindObject(h_true_M_0->GetName()))->Fill( true_V[0]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_0->GetName()))->Fill( true_V[0]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_0->GetName()))->Fill(true_V[0]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_0->GetName()))->Fill( true_V[0]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_1->GetName()))->Fill( true_V[1]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_1->GetName()))->Fill( true_V[1]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_1->GetName()))->Fill(true_V[1]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_1->GetName()))->Fill( true_V[1]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_2->GetName()))->Fill( true_V[3]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_2->GetName()))->Fill( true_V[3]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_2->GetName()))->Fill(true_V[3]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_2->GetName()))->Fill( true_V[3]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_3->GetName()))->Fill( true_V[3]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_3->GetName()))->Fill( true_V[3]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_3->GetName()))->Fill(true_V[3]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_3->GetName()))->Fill( true_V[3]->E());

      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt0->GetName()))->Fill(pvtx[0]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt1->GetName()))->Fill(pvtx[1]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt2->GetName()))->Fill(pvtx[2]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt3->GetName()))->Fill(pvtx[3]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt0->GetName()))->Fill(dvtx[0]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt1->GetName()))->Fill(dvtx[1]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt2->GetName()))->Fill(dvtx[2]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt3->GetName()))->Fill(dvtx[3]);
    }
    if(Npart==5){
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt4->GetName()))->Fill(pvtx[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt4->GetName()))->Fill(dvtx[4]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt4->GetName()))->Fill(pType[4]);
      ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_V[4]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_V[4]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_V[4]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_V[4]->E());

    }
    if(Npart==6){
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt4->GetName()))->Fill(pvtx[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt4->GetName()))->Fill(dvtx[4]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt4->GetName()))->Fill(pType[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt5->GetName()))->Fill(pvtx[5]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt5->GetName()))->Fill(dvtx[5]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt5->GetName()))->Fill(pType[5]);

      ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_V[4]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_V[4]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_V[4]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_V[4]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_5->GetName()))->Fill( true_V[5]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_5->GetName()))->Fill( true_V[5]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_5->GetName()))->Fill(true_V[5]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_5->GetName()))->Fill( true_V[5]->E());

    }
    if(Npart==7){
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt4->GetName()))->Fill(pvtx[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt4->GetName()))->Fill(pvtx[4]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt4->GetName()))->Fill(pType[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt5->GetName()))->Fill(pvtx[5]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt5->GetName()))->Fill(pvtx[5]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt5->GetName()))->Fill(pType[5]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt6->GetName()))->Fill(pvtx[6]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt6->GetName()))->Fill(pvtx[6]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt6->GetName()))->Fill(pType[6]);

      ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_V[4]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_V[4]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_V[4]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_V[4]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_5->GetName()))->Fill( true_V[5]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_5->GetName()))->Fill( true_V[5]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_5->GetName()))->Fill(true_V[5]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_5->GetName()))->Fill( true_V[5]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_6->GetName()))->Fill( true_V[6]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_6->GetName()))->Fill( true_V[6]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_6->GetName()))->Fill(true_V[6]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_6->GetName()))->Fill( true_V[6]->E());

    }
    if(Npart==8){
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt4->GetName()))->Fill(pvtx[pType[4]]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt4->GetName()))->Fill(dvtx[pType[4]]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt4->GetName()))->Fill(pType[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt5->GetName()))->Fill(pvtx[pType[5]]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt5->GetName()))->Fill(dvtx[pType[5]]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt5->GetName()))->Fill(pType[5]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt6->GetName()))->Fill(pvtx[pType[6]]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt6->GetName()))->Fill( dvtx[pType[6]]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt6->GetName()))->Fill(pType[6]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt7->GetName()))->Fill(pvtx[pType[7]]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt7->GetName()))->Fill( dvtx[pType[7]]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt7->GetName()))->Fill(pType[7]);


      ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_V[4]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_V[4]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_V[4]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_V[4]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_5->GetName()))->Fill( true_V[5]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_5->GetName()))->Fill( true_V[5]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_5->GetName()))->Fill(true_V[5]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_5->GetName()))->Fill( true_V[5]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_6->GetName()))->Fill( true_V[6]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_6->GetName()))->Fill( true_V[6]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_6->GetName()))->Fill(true_V[6]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_6->GetName()))->Fill( true_V[6]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_7->GetName()))->Fill( true_V[7]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_7->GetName()))->Fill( true_V[7]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_7->GetName()))->Fill(true_V[7]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_7->GetName()))->Fill( true_V[7]->E());


    }
    //if(pType[7]!=-1)
    //  ((TH1F*)gDirectory->FindObject(h_mc_vtxdiff_pi_mu2->GetName()))->Fill(pvtx[pType[7]]);
  }

   if(((sevt->trigWord)>>3)&0x1) {
    // ref. trigg. ok:
    ((TH1I*)gDirectory->FindObject(h_CPRE->GetName()))->Fill(evttypebytrks);
    if( (sevt->pu[3].chan[5]>>4)&0x1 ||
	(sevt->pu[4].chan[5]>>4)&0x1 ||
	(sevt->pu[5].chan[5]>>4)&0x1 ||  //5 is the central time slice
	(sevt->pu[6].chan[5]>>4)&0x1 ||
	(sevt->pu[3].chan[5]>>8)&0x1 ||
	(sevt->pu[4].chan[5]>>8)&0x1 ||
	(sevt->pu[5].chan[5]>>8)&0x1 ||  //5 is the central time slice
	(sevt->pu[6].chan[5]>>8)&0x1 ) {
      if(((sevt->trigWord))&0x1   || ((sevt->trigWord)>>2)&0x1  ) {
	//My trigger also should be OK (constructed from PU)
	((TH1I*)gDirectory->FindObject(h_full_trig->GetName()))->Fill(evttypebytrks);
      }
      //tr_eff = hreal/href;
    }
  }
   if( (sevt->pu[3].chan[5]>>4)&0x1 ||
       (sevt->pu[4].chan[5]>>4)&0x1 ||
       (sevt->pu[5].chan[5]>>4)&0x1 ||  //5 is the central time slice
       (sevt->pu[6].chan[5]>>4)&0x1 ||
       (sevt->pu[3].chan[5]>>8)&0x1 ||
       (sevt->pu[4].chan[5]>>8)&0x1 ||
       (sevt->pu[5].chan[5]>>8)&0x1 ||  //5 is the central time slice
       (sevt->pu[6].chan[5]>>8)&0x1 ) {
     if(((sevt->trigWord))&0x1   || ((sevt->trigWord)>>2)&0x1  ) {
       //My trigger also should be OK (constructed from PU)
       ((TH1I*)gDirectory->FindObject(h_mu1->GetName()))->Fill(evttypebytrks);
     }
     //tr_eff = hreal/href;
   }
   if(
      (sevt->pu[3].chan[5]>>8)&0x1 ||
      (sevt->pu[4].chan[5]>>8)&0x1 ||   //1vtx
      (sevt->pu[5].chan[5]>>8)&0x1 ||  //5 is the central time slice
      (sevt->pu[6].chan[5]>>8)&0x1 ) {
     if( ((sevt->trigWord)>>2)&0x1  ) {
       //My trigger also should be OK (constructed from PU)
       ((TH1I*)gDirectory->FindObject(h_1vtx_mcut->GetName()))->Fill(evttypebytrks);
     }
     //tr_eff = hreal/href;
   }
   if( (sevt->pu[3].chan[5]>>4)&0x1 ||
       (sevt->pu[4].chan[5]>>4)&0x1 ||//2vtx
       (sevt->pu[5].chan[5]>>4)&0x1 ||  //5 is the central time slice
       (sevt->pu[6].chan[5]>>4)&0x1 ) {
     if(((sevt->trigWord))&0x1 ) {
       //My trigger also should be OK (constructed from PU)
       ((TH1I*)gDirectory->FindObject(h_2vtx->GetName()))->Fill(evttypebytrks);
     }
     //tr_eff = hreal/href;
   }


  ((TH1F*)gDirectory->FindObject( h_vz_initial->GetName()))->Fill( sevt->vtx[0].z);

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

   //Trk hodtimediff
  ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( el1trkhodtime-el2trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( el1trkhodtime-muchtrkhodtime );
  ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( el2trkhodtime-muchtrkhodtime );

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

   ((TH1F*)gDirectory->FindObject(h_gmueeE->GetName()))->Fill( (*mu_nu).E());
  ((TH1F*)gDirectory->FindObject(h_gmueeP->GetName()))->Fill( (*mu_nu).P());
  ((TH1F*)gDirectory->FindObject(h_gmueePt->GetName()))->Fill((*mu_nu).Pt());
  ((TH1F*)gDirectory->FindObject(h_gmueeM->GetName()))->Fill ((*mu_nu).M());
  ((TH1F*)gDirectory->FindObject(h_bx_much->GetName()))            ->Fill(muchp_uncorr[0]);
  ((TH1F*)gDirectory->FindObject(h_by_much->GetName()))            ->Fill(muchp_uncorr[1]);
  ((TH1F*)gDirectory->FindObject(h_bdxdz_much->GetName()))            ->Fill(muchv_uncorr[0]);
  ((TH1F*)gDirectory->FindObject(h_bdydz_much->GetName()))            ->Fill(muchv_uncorr[1]) ;
  ((TH1F*)gDirectory->FindObject(h_bx_vs_by_much->GetName()))        ->Fill(muchp_uncorr[0], muchp_uncorr[1]);
  ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_much->GetName()))->Fill(muchv_uncorr[0], muchv_uncorr[1]);

 ((TH1F*)gDirectory->FindObject(h_k3pi_nu_m->GetName()))->Fill( (*muvee_pi).M() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_nu_P->GetName()))->Fill( (*muvee_pi).P() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_nu_E->GetName()))->Fill( (*muvee_pi).E() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_nu_Pt->GetName()))->Fill( (*muvee_pi).Pt() );
   ((TH1F*)gDirectory->FindObject(h_muee_piee_M->GetName()))->Fill( (*muvee_piee).M() );
  ((TH1F*)gDirectory->FindObject(h_muee_piee_P->GetName()))->Fill( (*muvee_piee).P() );
  ((TH1F*)gDirectory->FindObject(h_muee_piee_E->GetName()))->Fill( (*muvee_piee).E() );
  ((TH1F*)gDirectory->FindObject(h_muee_piee_Pt->GetName()))->Fill( (*muvee_piee).Pt() );

  //((TH1F*)gDirectory->FindObject(h_3trk_invm->GetName()))->Fill ((*k3pic).M2());
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_P->GetName()))->Fill( (*muchlv_pi).P() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_P->GetName()))->Fill( (*el1lv_pi).P() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_P->GetName()))->Fill( (*el2lv_pi).P() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_E->GetName()))->Fill( (*muchlv_pi).E() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_E->GetName()))->Fill( (*el1lv_pi).E() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_E->GetName()))->Fill( (*el2lv_pi).E() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_Pt->GetName()))->Fill((*muchlv_pi).Pt() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_Pt->GetName()))->Fill((*el1lv_pi).Pt() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_Pt->GetName()))->Fill((*el2lv_pi).Pt() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_m->GetName()))->Fill( (*muchlv_pi).M() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_m->GetName()))->Fill( (*el1lv_pi).M() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_m->GetName()))->Fill( (*el2lv_pi).M() );






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
	//muchtrkq==-1        &&
	eelvmass >= 0.140   &&
	(*muvee_pi).M() >= 0.51 &&
	fabs(el1trktime - el1trkhodtime)<12.   &&
	fabs(el2trktime - el2trkhodtime)<12.   &&
	fabs(muchtrktime-muchtrkhodtime)< 12. &&
	fabs((long double)(el1trktime-muchtrktime)) < 10.   &&
	fabs((long double)(el2trktime-muchtrktime)) < 10.   &&
	el1cltime-el2cltime > -3 && el1cltime-el2cltime < 3 &&
     	el1trkhodtime-el2trkhodtime   > -2 && el1trkhodtime-el2trkhodtime    < 2     &&
	el1trkhodtime-muchtrkhodtime  > -2 && el1trkhodtime - muchtrkhodtime < 2 &&
	el2trkhodtime-muchtrkhodtime  > -2 && el2trkhodtime - muchtrkhodtime < 2
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

      ((TH1F*)gDirectory->FindObject( h_vz_initial->GetName()))->Fill( sevt->vtx[0].z);
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

      //Trk hodtimediff
      ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( el1trkhodtime-el2trkhodtime );
      ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( el1trkhodtime-muchtrkhodtime );
      ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( el2trkhodtime-muchtrkhodtime );

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

       ((TH1F*)gDirectory->FindObject(h_gmueeE->GetName()))->Fill( (*mu_nu).E());
  ((TH1F*)gDirectory->FindObject(h_gmueeP->GetName()))->Fill( (*mu_nu).P());
  ((TH1F*)gDirectory->FindObject(h_gmueePt->GetName()))->Fill((*mu_nu).Pt());
  ((TH1F*)gDirectory->FindObject(h_gmueeM->GetName()))->Fill ((*mu_nu).M());
      ((TH1F*)gDirectory->FindObject(h_bx_much->GetName()))            ->Fill(muchp_uncorr[0]);
      ((TH1F*)gDirectory->FindObject(h_by_much->GetName()))            ->Fill(muchp_uncorr[1]);
      ((TH1F*)gDirectory->FindObject(h_bdxdz_much->GetName()))            ->Fill(muchv_uncorr[0]);
      ((TH1F*)gDirectory->FindObject(h_bdydz_much->GetName()))            ->Fill(muchv_uncorr[1]) ;
      ((TH1F*)gDirectory->FindObject(h_bx_vs_by_much->GetName()))        ->Fill(muchp_uncorr[0], muchp_uncorr[1]);
      ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_much->GetName()))->Fill(muchv_uncorr[0], muchv_uncorr[1]);



      ((TH1F*)gDirectory->FindObject(h_k3pi_nu_m->GetName()))->Fill( (*muvee_pi).M() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_nu_P->GetName()))->Fill( (*muvee_pi).P() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_nu_E->GetName()))->Fill( (*muvee_pi).E() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_nu_Pt->GetName()))->Fill( (*muvee_pi).Pt() );
      ((TH1F*)gDirectory->FindObject(h_muee_piee_M->GetName()))->Fill( (*muvee_piee).M() );
      ((TH1F*)gDirectory->FindObject(h_muee_piee_P->GetName()))->Fill( (*muvee_piee).P() );
      ((TH1F*)gDirectory->FindObject(h_muee_piee_E->GetName()))->Fill( (*muvee_piee).E() );
      ((TH1F*)gDirectory->FindObject(h_muee_piee_Pt->GetName()))->Fill( (*muvee_piee).Pt() );

      //((TH1F*)gDirectory->FindObject(h_3trk_invm->GetName()))->Fill ((*k3pic).M2());
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_P->GetName()))->Fill( (*muchlv_pi).P() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_P->GetName()))->Fill( (*el1lv_pi).P() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_P->GetName()))->Fill( (*el2lv_pi).P() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_E->GetName()))->Fill( (*muchlv_pi).E() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_E->GetName()))->Fill( (*el1lv_pi).E() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_E->GetName()))->Fill( (*el2lv_pi).E() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_Pt->GetName()))->Fill((*muchlv_pi).Pt() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_Pt->GetName()))->Fill((*el1lv_pi).Pt() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_Pt->GetName()))->Fill((*el2lv_pi).Pt() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_m->GetName()))->Fill( (*muchlv_pi).M() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_m->GetName()))->Fill( (*el1lv_pi).M() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_m->GetName()))->Fill( (*el2lv_pi).M() );



    }
  }

  file1->cd();
  h_cutflow->Fill(cutcounter);
  h_cutflow->GetXaxis()->SetBinLabel( cutcounter, cuts[cutcounter-1].c_str() );
  cutcounter++;
  ///THE END OF ************* Final k3pi selection  ****************** ///


  /// ************* Particle identification ****************** ///
  if (eltrkq1==eltrkq2)    {cleanup();return 0;}
  if (muchtrkindx==-1 )    {cleanup();return 0;}

  dirs[cutcounter-1]->cd();

  if (IS_MC){
    ////if(true_4.P()  < 3. || true_5.P() < 3.) {cleanup();return 0;}
    ((TH1F*)gDirectory->FindObject(h_mc_Npart_tr->GetName()))->Fill(Npart);
    if (Npart >=4 ){

      ((TH1D*)gDirectory->FindObject( h_mc_muvee_M->GetName()))->Fill( Muee_three_track.M());
      ((TH1D*)gDirectory->FindObject( h_mc_muvee_P->GetName()))->Fill( Muee_three_track.P());
      ((TH1D*)gDirectory->FindObject( h_mc_muvee_Pt->GetName()))->Fill( Muee_three_track.Pt());
      ((TH1D*)gDirectory->FindObject( h_mc_elsP->GetName()))->Fill( mc_ee.P());
      ((TH1D*)gDirectory->FindObject( h_mc_elsM->GetName()))->Fill( mc_ee.M());
      ((TH1D*)gDirectory->FindObject( h_mc_nu_M->GetName()))->Fill( Nu_true.M());
      ((TH1D*)gDirectory->FindObject( h_mc_nu_M2->GetName()))->Fill( Nu_true.M2());
      ((TH1D*)gDirectory->FindObject( h_mc_nu_P->GetName()))->Fill( Nu_true.P());
      ((TH1D*)gDirectory->FindObject( h_mc_nu_Pt->GetName()))->Fill(Nu_true.Pt());
      ((TH1D*)gDirectory->FindObject(h_true_M_0->GetName()))->Fill( true_V[0]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_0->GetName()))->Fill( true_V[0]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_0->GetName()))->Fill(true_V[0]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_0->GetName()))->Fill( true_V[0]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_1->GetName()))->Fill( true_V[1]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_1->GetName()))->Fill( true_V[1]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_1->GetName()))->Fill(true_V[1]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_1->GetName()))->Fill( true_V[1]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_2->GetName()))->Fill( true_V[3]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_2->GetName()))->Fill( true_V[3]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_2->GetName()))->Fill(true_V[3]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_2->GetName()))->Fill( true_V[3]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_3->GetName()))->Fill( true_V[3]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_3->GetName()))->Fill( true_V[3]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_3->GetName()))->Fill(true_V[3]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_3->GetName()))->Fill( true_V[3]->E());

      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt0->GetName()))->Fill(pvtx[0]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt1->GetName()))->Fill(pvtx[1]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt2->GetName()))->Fill(pvtx[2]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt3->GetName()))->Fill(pvtx[3]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt0->GetName()))->Fill(dvtx[0]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt1->GetName()))->Fill(dvtx[1]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt2->GetName()))->Fill(dvtx[2]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt3->GetName()))->Fill(dvtx[3]);
    }
    if(Npart==5){
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt4->GetName()))->Fill(pvtx[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt4->GetName()))->Fill(dvtx[4]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt4->GetName()))->Fill(pType[4]);
      ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_V[4]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_V[4]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_V[4]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_V[4]->E());

    }
    if(Npart==6){
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt4->GetName()))->Fill(pvtx[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt4->GetName()))->Fill(dvtx[4]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt4->GetName()))->Fill(pType[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt5->GetName()))->Fill(pvtx[5]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt5->GetName()))->Fill(dvtx[5]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt5->GetName()))->Fill(pType[5]);

      ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_V[4]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_V[4]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_V[4]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_V[4]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_5->GetName()))->Fill( true_V[5]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_5->GetName()))->Fill( true_V[5]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_5->GetName()))->Fill(true_V[5]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_5->GetName()))->Fill( true_V[5]->E());

    }
    if(Npart==7){
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt4->GetName()))->Fill(pvtx[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt4->GetName()))->Fill(pvtx[4]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt4->GetName()))->Fill(pType[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt5->GetName()))->Fill(pvtx[5]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt5->GetName()))->Fill(pvtx[5]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt5->GetName()))->Fill(pType[5]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt6->GetName()))->Fill(pvtx[6]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt6->GetName()))->Fill(pvtx[6]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt6->GetName()))->Fill(pType[6]);

      ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_V[4]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_V[4]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_V[4]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_V[4]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_5->GetName()))->Fill( true_V[5]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_5->GetName()))->Fill( true_V[5]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_5->GetName()))->Fill(true_V[5]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_5->GetName()))->Fill( true_V[5]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_6->GetName()))->Fill( true_V[6]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_6->GetName()))->Fill( true_V[6]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_6->GetName()))->Fill(true_V[6]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_6->GetName()))->Fill( true_V[6]->E());

    }
    if(Npart==8){
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt4->GetName()))->Fill(pvtx[pType[4]]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt4->GetName()))->Fill(dvtx[pType[4]]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt4->GetName()))->Fill(pType[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt5->GetName()))->Fill(pvtx[pType[5]]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt5->GetName()))->Fill(dvtx[pType[5]]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt5->GetName()))->Fill(pType[5]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt6->GetName()))->Fill(pvtx[pType[6]]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt6->GetName()))->Fill( dvtx[pType[6]]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt6->GetName()))->Fill(pType[6]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt7->GetName()))->Fill(pvtx[pType[7]]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt7->GetName()))->Fill( dvtx[pType[7]]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt7->GetName()))->Fill(pType[7]);


      ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_V[4]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_V[4]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_V[4]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_V[4]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_5->GetName()))->Fill( true_V[5]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_5->GetName()))->Fill( true_V[5]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_5->GetName()))->Fill(true_V[5]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_5->GetName()))->Fill( true_V[5]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_6->GetName()))->Fill( true_V[6]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_6->GetName()))->Fill( true_V[6]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_6->GetName()))->Fill(true_V[6]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_6->GetName()))->Fill( true_V[6]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_7->GetName()))->Fill( true_V[7]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_7->GetName()))->Fill( true_V[7]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_7->GetName()))->Fill(true_V[7]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_7->GetName()))->Fill( true_V[7]->E());


    }
    //if(pType[7]!=-1)
    //  ((TH1F*)gDirectory->FindObject(h_mc_vtxdiff_pi_mu2->GetName()))->Fill(pvtx[pType[7]]);
  }

  if(((sevt->trigWord)>>3)&0x1) {
    // ref. trigg. ok:
    ((TH1I*)gDirectory->FindObject(h_CPRE->GetName()))->Fill(evttypebytrks);
    if( (sevt->pu[3].chan[5]>>4)&0x1 ||
	(sevt->pu[4].chan[5]>>4)&0x1 ||
	(sevt->pu[5].chan[5]>>4)&0x1 ||  //5 is the central time slice
	(sevt->pu[6].chan[5]>>4)&0x1 ||
	(sevt->pu[3].chan[5]>>8)&0x1 ||
	(sevt->pu[4].chan[5]>>8)&0x1 ||
	(sevt->pu[5].chan[5]>>8)&0x1 ||  //5 is the central time slice
	(sevt->pu[6].chan[5]>>8)&0x1 ) {
      if(((sevt->trigWord))&0x1   || ((sevt->trigWord)>>2)&0x1  ) {
	//My trigger also should be OK (constructed from PU)
	((TH1I*)gDirectory->FindObject(h_full_trig->GetName()))->Fill(evttypebytrks);
      }
      //tr_eff = hreal/href;
    }
  }
   if(
      (sevt->pu[3].chan[5]>>8)&0x1 ||
      (sevt->pu[4].chan[5]>>8)&0x1 ||   //1vtx
      (sevt->pu[5].chan[5]>>8)&0x1 ||  //5 is the central time slice
      (sevt->pu[6].chan[5]>>8)&0x1 ) {
     if( ((sevt->trigWord)>>2)&0x1  ) {
       //My trigger also should be OK (constructed from PU)
       ((TH1I*)gDirectory->FindObject(h_1vtx_mcut->GetName()))->Fill(evttypebytrks);
     }
     //tr_eff = hreal/href;
   }
   if( (sevt->pu[3].chan[5]>>4)&0x1 ||
       (sevt->pu[4].chan[5]>>4)&0x1 ||//2vtx
       (sevt->pu[5].chan[5]>>4)&0x1 ||  //5 is the central time slice
       (sevt->pu[6].chan[5]>>4)&0x1 ) {
     if(((sevt->trigWord))&0x1 ) {
       //My trigger also should be OK (constructed from PU)
       ((TH1I*)gDirectory->FindObject(h_2vtx->GetName()))->Fill(evttypebytrks);
     }
     //tr_eff = hreal/href;
   }
   if( (sevt->pu[3].chan[5]>>4)&0x1 ||
       (sevt->pu[4].chan[5]>>4)&0x1 ||
       (sevt->pu[5].chan[5]>>4)&0x1 ||  //5 is the central time slice
       (sevt->pu[6].chan[5]>>4)&0x1 ||
       (sevt->pu[3].chan[5]>>8)&0x1 ||
       (sevt->pu[4].chan[5]>>8)&0x1 ||
       (sevt->pu[5].chan[5]>>8)&0x1 ||  //5 is the central time slice
       (sevt->pu[6].chan[5]>>8)&0x1 ) {
     if(((sevt->trigWord))&0x1   || ((sevt->trigWord)>>2)&0x1  ) {
       //My trigger also should be OK (constructed from PU)
       ((TH1I*)gDirectory->FindObject(h_mu1->GetName()))->Fill(evttypebytrks);
     }
     //tr_eff = hreal/href;
   }

  ((TH1F*)gDirectory->FindObject( h_vz_initial->GetName()))->Fill( sevt->vtx[0].z);
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
   //Trk hodtimediff
  ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( el1trkhodtime-el2trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( el1trkhodtime-muchtrkhodtime );
  ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( el2trkhodtime-muchtrkhodtime );


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
   ((TH1F*)gDirectory->FindObject(h_gmueeE->GetName()))->Fill( (*mu_nu).E());
  ((TH1F*)gDirectory->FindObject(h_gmueeP->GetName()))->Fill( (*mu_nu).P());
  ((TH1F*)gDirectory->FindObject(h_gmueePt->GetName()))->Fill((*mu_nu).Pt());
  ((TH1F*)gDirectory->FindObject(h_gmueeM->GetName()))->Fill ((*mu_nu).M());

  ((TH1F*)gDirectory->FindObject(h_bx_much->GetName()))            ->Fill(muchp_uncorr[0]);
  ((TH1F*)gDirectory->FindObject(h_by_much->GetName()))            ->Fill(muchp_uncorr[1]);
  ((TH1F*)gDirectory->FindObject(h_bdxdz_much->GetName()))            ->Fill(muchv_uncorr[0]);
  ((TH1F*)gDirectory->FindObject(h_bdydz_much->GetName()))            ->Fill(muchv_uncorr[1]) ;
  ((TH1F*)gDirectory->FindObject(h_bx_vs_by_much->GetName()))        ->Fill(muchp_uncorr[0], muchp_uncorr[1]);
  ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_much->GetName()))->Fill(muchv_uncorr[0], muchv_uncorr[1]);
  ((TH1F*)gDirectory->FindObject(h_angle_e1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()) );
  ((TH1F*)gDirectory->FindObject(h_angle_e1e2_mu->GetName()))->Fill( (*eelv).Vect().Angle(muchlv->Vect()) );

   ((TH1F*)gDirectory->FindObject(h_k3pi_nu_m->GetName()))->Fill( (*muvee_pi).M() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_nu_P->GetName()))->Fill( (*muvee_pi).P() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_nu_E->GetName()))->Fill( (*muvee_pi).E() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_nu_Pt->GetName()))->Fill( (*muvee_pi).Pt() );
 ((TH1F*)gDirectory->FindObject(h_muee_piee_M->GetName()))->Fill( (*muvee_piee).M() );
  ((TH1F*)gDirectory->FindObject(h_muee_piee_P->GetName()))->Fill( (*muvee_piee).P() );
  ((TH1F*)gDirectory->FindObject(h_muee_piee_E->GetName()))->Fill( (*muvee_piee).E() );
  ((TH1F*)gDirectory->FindObject(h_muee_piee_Pt->GetName()))->Fill( (*muvee_piee).Pt() );

  //((TH1F*)gDirectory->FindObject(h_3trk_invm->GetName()))->Fill ((*k3pic).M2());
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_P->GetName()))->Fill( (*muchlv_pi).P() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_P->GetName()))->Fill( (*el1lv_pi).P() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_P->GetName()))->Fill( (*el2lv_pi).P() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_E->GetName()))->Fill( (*muchlv_pi).E() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_E->GetName()))->Fill( (*el1lv_pi).E() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_E->GetName()))->Fill( (*el2lv_pi).E() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_Pt->GetName()))->Fill((*muchlv_pi).Pt() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_Pt->GetName()))->Fill((*el1lv_pi).Pt() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_Pt->GetName()))->Fill((*el2lv_pi).Pt() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_m->GetName()))->Fill( (*muchlv_pi).M() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_m->GetName()))->Fill( (*el1lv_pi).M() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_m->GetName()))->Fill( (*el2lv_pi).M() );




  file1->cd();
  h_cutflow->Fill(cutcounter);
  h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());
  cutcounter++;
  ///THE END OF ************* Particle identification ****************** ///

  /// ************* Inv mass -140 ****************** ///

  if (eelvmass <= 0.140)  {cleanup();return 0;}
  if (IS_DATA)
    if(eop_e2_track <=0.95 || eop_e1_track <=0.95 ) {cleanup();return 0;}

  //if (flux_id==1){
  //  if ((*muvee).M() < 0.470 || (*muvee).M() > 0.505 ) {cleanup();return 0;}
  //  if ((*muvee).P() < 54. || (*muvee).P() > 66 ) {cleanup();return 0;}
  //}
  dirs[cutcounter-1]->cd();

  if (IS_MC){
    ////if(true_4.P()  < 3. || true_5.P() < 3.) {cleanup();return 0;}
    ((TH1F*)gDirectory->FindObject(h_mc_Npart_tr->GetName()))->Fill(Npart);
    if (Npart >=4 ){

      ((TH1D*)gDirectory->FindObject( h_mc_muvee_M->GetName()))->Fill( Muee_three_track.M());
      ((TH1D*)gDirectory->FindObject( h_mc_muvee_P->GetName()))->Fill( Muee_three_track.P());
      ((TH1D*)gDirectory->FindObject( h_mc_muvee_Pt->GetName()))->Fill( Muee_three_track.Pt());
      ((TH1D*)gDirectory->FindObject( h_mc_elsP->GetName()))->Fill( mc_ee.P());
      ((TH1D*)gDirectory->FindObject( h_mc_elsM->GetName()))->Fill( mc_ee.M());
      ((TH1D*)gDirectory->FindObject( h_mc_nu_M->GetName()))->Fill( Nu_true.M());
      ((TH1D*)gDirectory->FindObject( h_mc_nu_M2->GetName()))->Fill( Nu_true.M2());
      ((TH1D*)gDirectory->FindObject( h_mc_nu_P->GetName()))->Fill( Nu_true.P());
      ((TH1D*)gDirectory->FindObject( h_mc_nu_Pt->GetName()))->Fill(Nu_true.Pt());
      ((TH1D*)gDirectory->FindObject(h_true_M_0->GetName()))->Fill( true_V[0]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_0->GetName()))->Fill( true_V[0]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_0->GetName()))->Fill(true_V[0]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_0->GetName()))->Fill( true_V[0]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_1->GetName()))->Fill( true_V[1]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_1->GetName()))->Fill( true_V[1]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_1->GetName()))->Fill(true_V[1]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_1->GetName()))->Fill( true_V[1]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_2->GetName()))->Fill( true_V[3]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_2->GetName()))->Fill( true_V[3]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_2->GetName()))->Fill(true_V[3]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_2->GetName()))->Fill( true_V[3]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_3->GetName()))->Fill( true_V[3]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_3->GetName()))->Fill( true_V[3]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_3->GetName()))->Fill(true_V[3]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_3->GetName()))->Fill( true_V[3]->E());

      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt0->GetName()))->Fill(pvtx[0]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt1->GetName()))->Fill(pvtx[1]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt2->GetName()))->Fill(pvtx[2]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt3->GetName()))->Fill(pvtx[3]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt0->GetName()))->Fill(dvtx[0]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt1->GetName()))->Fill(dvtx[1]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt2->GetName()))->Fill(dvtx[2]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt3->GetName()))->Fill(dvtx[3]);
    }
    if(Npart==5){
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt4->GetName()))->Fill(pvtx[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt4->GetName()))->Fill(dvtx[4]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt4->GetName()))->Fill(pType[4]);
      ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_V[4]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_V[4]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_V[4]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_V[4]->E());

    }
    if(Npart==6){
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt4->GetName()))->Fill(pvtx[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt4->GetName()))->Fill(dvtx[4]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt4->GetName()))->Fill(pType[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt5->GetName()))->Fill(pvtx[5]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt5->GetName()))->Fill(dvtx[5]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt5->GetName()))->Fill(pType[5]);

      ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_V[4]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_V[4]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_V[4]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_V[4]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_5->GetName()))->Fill( true_V[5]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_5->GetName()))->Fill( true_V[5]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_5->GetName()))->Fill(true_V[5]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_5->GetName()))->Fill( true_V[5]->E());

    }
    if(Npart==7){
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt4->GetName()))->Fill(pvtx[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt4->GetName()))->Fill(pvtx[4]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt4->GetName()))->Fill(pType[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt5->GetName()))->Fill(pvtx[5]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt5->GetName()))->Fill(pvtx[5]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt5->GetName()))->Fill(pType[5]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt6->GetName()))->Fill(pvtx[6]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt6->GetName()))->Fill(pvtx[6]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt6->GetName()))->Fill(pType[6]);

      ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_V[4]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_V[4]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_V[4]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_V[4]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_5->GetName()))->Fill( true_V[5]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_5->GetName()))->Fill( true_V[5]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_5->GetName()))->Fill(true_V[5]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_5->GetName()))->Fill( true_V[5]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_6->GetName()))->Fill( true_V[6]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_6->GetName()))->Fill( true_V[6]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_6->GetName()))->Fill(true_V[6]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_6->GetName()))->Fill( true_V[6]->E());

    }
    if(Npart==8){
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt4->GetName()))->Fill(pvtx[pType[4]]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt4->GetName()))->Fill(dvtx[pType[4]]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt4->GetName()))->Fill(pType[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt5->GetName()))->Fill(pvtx[pType[5]]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt5->GetName()))->Fill(dvtx[pType[5]]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt5->GetName()))->Fill(pType[5]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt6->GetName()))->Fill(pvtx[pType[6]]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt6->GetName()))->Fill( dvtx[pType[6]]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt6->GetName()))->Fill(pType[6]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt7->GetName()))->Fill(pvtx[pType[7]]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt7->GetName()))->Fill( dvtx[pType[7]]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt7->GetName()))->Fill(pType[7]);


      ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_V[4]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_V[4]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_V[4]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_V[4]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_5->GetName()))->Fill( true_V[5]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_5->GetName()))->Fill( true_V[5]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_5->GetName()))->Fill(true_V[5]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_5->GetName()))->Fill( true_V[5]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_6->GetName()))->Fill( true_V[6]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_6->GetName()))->Fill( true_V[6]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_6->GetName()))->Fill(true_V[6]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_6->GetName()))->Fill( true_V[6]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_7->GetName()))->Fill( true_V[7]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_7->GetName()))->Fill( true_V[7]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_7->GetName()))->Fill(true_V[7]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_7->GetName()))->Fill( true_V[7]->E());


    }
    //if(pType[7]!=-1)
    //  ((TH1F*)gDirectory->FindObject(h_mc_vtxdiff_pi_mu2->GetName()))->Fill(pvtx[pType[7]]);
  }

//epsilon Q
   if(((sevt->trigWord)>>3)&0x1) {
    // ref. trigg. ok:
    ((TH1I*)gDirectory->FindObject(h_CPRE->GetName()))->Fill(evttypebytrks);
    if( (sevt->pu[3].chan[5]>>4)&0x1 ||
	(sevt->pu[4].chan[5]>>4)&0x1 ||//2vtx
	(sevt->pu[5].chan[5]>>4)&0x1 ||  //5 is the central time slice
	(sevt->pu[6].chan[5]>>4)&0x1 ||
	(sevt->pu[3].chan[5]>>8)&0x1 ||//1vtx
	(sevt->pu[4].chan[5]>>8)&0x1 ||
	(sevt->pu[5].chan[5]>>8)&0x1 ||  //5 is the central time slice
	(sevt->pu[6].chan[5]>>8)&0x1 ||
	(sevt->pu[3].chan[5]>>13)&0x1 ||
	(sevt->pu[4].chan[5]>>13)&0x1 ||
	(sevt->pu[5].chan[5]>>13)&0x1 ||  //5 is the central time slice
	(sevt->pu[6].chan[5]>>13)&0x1) {
      if(((sevt->trigWord))&0x1   || ((sevt->trigWord)>>2)&0x1 || ((sevt->trigWord)>>5)&0x1 ) {
	//My trigger also should be OK (constructed from PU)
	((TH1I*)gDirectory->FindObject(h_full_trig->GetName()))->Fill(evttypebytrks);

      }
      //tr_eff = hreal/href;
    }
  }

  //epsilon a
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
      (sevt->pu[6].chan[5]>>13)&0x1
      ) {
    if(((sevt->trigWord))&0x1   || ((sevt->trigWord)>>2)&0x1 || ((sevt->trigWord)>>5)&0x1 ) {
      //My trigger also should be OK (constructed from PU)
      ((TH1I*)gDirectory->FindObject(h_mu1->GetName()))->Fill(evttypebytrks);
      if(((sevt->trigWord)>>3)&0x1) {
	((TH1I*)gDirectory->FindObject(h_Q1->GetName()))->Fill(evttypebytrks);
      }
    }
    //tr_eff = hreal/href;
  }


  //if (nEvt_counter==0){
  //  cout << "nEvt=" << sevt->nEvt << endl;
  //  first_Nevt=sevt->nEvt;
  //  v_first_e1= new TLorentzVector((*el1lv));
  //  v_first_e2= new TLorentzVector((*el2lv));
  //  v_first_mu= new TLorentzVector((*muchlv));
  //  f_Ntrack=sevt->Ntrack;
  //  f_NCl=sevt->Ncluster;
  //} else {
  //  second_Nevt = sevt->nEvt;
  //  v_second_e1= new TLorentzVector((*el1lv));
  //  v_second_e2= new TLorentzVector((*el2lv));
  //  v_second_mu= new TLorentzVector((*muchlv));
  //  s_Ntrack=sevt->Ntrack;
  //  s_NCl=sevt->Ncluster;
  //}
  //cout << "first_Nevt=" << first_Nevt << endl;
  //cout << " second_Nevt=" << second_Nevt << endl;
  int difference= second_Nevt - first_Nevt;

  //if(nEvt_counter >= 0){
  //
  //  if (fabs( difference ) <= 1.){
  //    cout << "Ntrack diff" << s_Ntrack - f_Ntrack << endl;
  //    cout << "NCluster diff" << s_NCl - f_NCl << endl;
  //    cout << "e1.P diff" << (*v_first_e1).P() - (*v_second_e1).P() << endl;
  //    cout << "e2.P diff" << (*v_first_e2).P() - (*v_second_e2).P() << endl;
  //    cout << "mu.P diff" << (*v_first_mu).P() - (*v_second_mu).P() << endl;
  //    cout << "e1.P diff" << (*v_first_e1).E() - (*v_second_e1).E() << endl;
  //    cout << "e2.P diff" << (*v_first_e2).E() - (*v_second_e2).E() << endl;
  //    cout << "mu.P diff" << (*v_first_mu).E() - (*v_second_mu).E() << endl;
  //    cout << "TUKA SUM" << endl;
  //    rep_counter++;
  //  }
  //  v_first_e1=v_second_e1;
  //  v_first_e2=v_second_e2;
  //  v_first_mu=v_second_mu;
  //  f_Ntrack  =s_Ntrack;
  //  f_NCl     =s_NCl;
  //  first_Nevt=second_Nevt;
  //}
  //
  //
  //
  //
  //
  //
  //nEvt_counter++;
  //cout << "nEvt_counter=" << nEvt_counter << endl;
  //cout << "rep_counter" << rep_counter << endl;




  ((TH1F*)gDirectory->FindObject( h_vz_initial->GetName()))->Fill( sevt->vtx[0].z);
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

   //Trk hodtimediff
  ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( el1trkhodtime-el2trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( el1trkhodtime-muchtrkhodtime );
  ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( el2trkhodtime-muchtrkhodtime );

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
   ((TH1F*)gDirectory->FindObject(h_gmueeE->GetName()))->Fill( (*mu_nu).E());
  ((TH1F*)gDirectory->FindObject(h_gmueeP->GetName()))->Fill( (*mu_nu).P());
  ((TH1F*)gDirectory->FindObject(h_gmueePt->GetName()))->Fill((*mu_nu).Pt());
  ((TH1F*)gDirectory->FindObject(h_gmueeM->GetName()))->Fill ((*mu_nu).M());

  ((TH1F*)gDirectory->FindObject(h_bx_much->GetName()))            ->Fill(muchp_uncorr[0]);
  ((TH1F*)gDirectory->FindObject(h_by_much->GetName()))            ->Fill(muchp_uncorr[1]);
  ((TH1F*)gDirectory->FindObject(h_bdxdz_much->GetName()))            ->Fill(muchv_uncorr[0]);
  ((TH1F*)gDirectory->FindObject(h_bdydz_much->GetName()))            ->Fill(muchv_uncorr[1]) ;
  ((TH1F*)gDirectory->FindObject(h_bx_vs_by_much->GetName()))        ->Fill(muchp_uncorr[0], muchp_uncorr[1]);
  ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_much->GetName()))->Fill(muchv_uncorr[0], muchv_uncorr[1]);

  ((TH1F*)gDirectory->FindObject(h_k3pi_nu_m->GetName()))->Fill( (*muvee_pi).M() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_nu_P->GetName()))->Fill( (*muvee_pi).P() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_nu_E->GetName()))->Fill( (*muvee_pi).E() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_nu_Pt->GetName()))->Fill( (*muvee_pi).Pt() );

  ((TH1F*)gDirectory->FindObject(h_muee_piee_M->GetName()))->Fill( (*muvee_piee).M() );
  ((TH1F*)gDirectory->FindObject(h_muee_piee_P->GetName()))->Fill( (*muvee_piee).P() );
  ((TH1F*)gDirectory->FindObject(h_muee_piee_E->GetName()))->Fill( (*muvee_piee).E() );
  ((TH1F*)gDirectory->FindObject(h_muee_piee_Pt->GetName()))->Fill( (*muvee_piee).Pt() );

  //((TH1F*)gDirectory->FindObject(h_3trk_invm->GetName()))->Fill ((*k3pic).M2());
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_P->GetName()))->Fill( (*muchlv_pi).P() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_P->GetName()))->Fill( (*el1lv_pi).P() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_P->GetName()))->Fill( (*el2lv_pi).P() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_E->GetName()))->Fill( (*muchlv_pi).E() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_E->GetName()))->Fill( (*el1lv_pi).E() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_E->GetName()))->Fill( (*el2lv_pi).E() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_Pt->GetName()))->Fill((*muchlv_pi).Pt() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_Pt->GetName()))->Fill((*el1lv_pi).Pt() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_Pt->GetName()))->Fill((*el2lv_pi).Pt() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_m->GetName()))->Fill( (*muchlv_pi).M() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_m->GetName()))->Fill( (*el1lv_pi).M() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_m->GetName()))->Fill( (*el2lv_pi).M() );





  file1->cd();
  h_cutflow->Fill(cutcounter);
  h_cutflow->GetXaxis()->SetBinLabel( cutcounter, cuts[cutcounter-1].c_str() );
  cutcounter++;
  ///THE END OF ************* Inv mass -140 ****************** ///

  /// ************* k3pi inv mass cut (at least one) ****************** ///
  //"muee energy > 70GeV test"  Mu transverse momentum > 0.9
  //if( (*muvee).E() < 70 )        {cleanup();return 0;}
  //if( (*muchlv).Pt() > 0.9 )        {
  //if( sevt->cluster[sevt->track[eltrkindx1].iClus].energy/sevt->track[eltrkindx1].p < 0.95 || sevt->cluster[sevt->track[eltrkindx2].iClus].energy/sevt->track[eltrkindx2].p < 0.95 )        {
  //
  if(  (*muvee_pi).M() <= 0.51 )  {cleanup();return 0;}
  dirs[cutcounter-1]->cd();

  if (IS_MC){
    ////if(true_4.P()  < 3. || true_5.P() < 3.) {cleanup();return 0;}
    ((TH1F*)gDirectory->FindObject(h_mc_Npart_tr->GetName()))->Fill(Npart);
    if (Npart >=4 ){

      ((TH1D*)gDirectory->FindObject( h_mc_muvee_M->GetName()))->Fill( Muee_three_track.M());
      ((TH1D*)gDirectory->FindObject( h_mc_muvee_P->GetName()))->Fill( Muee_three_track.P());
      ((TH1D*)gDirectory->FindObject( h_mc_muvee_Pt->GetName()))->Fill( Muee_three_track.Pt());
      ((TH1D*)gDirectory->FindObject( h_mc_elsP->GetName()))->Fill( mc_ee.P());
      ((TH1D*)gDirectory->FindObject( h_mc_elsM->GetName()))->Fill( mc_ee.M());
      ((TH1D*)gDirectory->FindObject( h_mc_nu_M->GetName()))->Fill( Nu_true.M());
      ((TH1D*)gDirectory->FindObject( h_mc_nu_M2->GetName()))->Fill( Nu_true.M2());
      ((TH1D*)gDirectory->FindObject( h_mc_nu_P->GetName()))->Fill( Nu_true.P());
      ((TH1D*)gDirectory->FindObject( h_mc_nu_Pt->GetName()))->Fill(Nu_true.Pt());
      ((TH1D*)gDirectory->FindObject(h_true_M_0->GetName()))->Fill( true_V[0]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_0->GetName()))->Fill( true_V[0]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_0->GetName()))->Fill(true_V[0]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_0->GetName()))->Fill( true_V[0]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_1->GetName()))->Fill( true_V[1]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_1->GetName()))->Fill( true_V[1]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_1->GetName()))->Fill(true_V[1]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_1->GetName()))->Fill( true_V[1]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_2->GetName()))->Fill( true_V[3]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_2->GetName()))->Fill( true_V[3]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_2->GetName()))->Fill(true_V[3]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_2->GetName()))->Fill( true_V[3]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_3->GetName()))->Fill( true_V[3]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_3->GetName()))->Fill( true_V[3]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_3->GetName()))->Fill(true_V[3]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_3->GetName()))->Fill( true_V[3]->E());

      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt0->GetName()))->Fill(pvtx[0]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt1->GetName()))->Fill(pvtx[1]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt2->GetName()))->Fill(pvtx[2]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt3->GetName()))->Fill(pvtx[3]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt0->GetName()))->Fill(dvtx[0]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt1->GetName()))->Fill(dvtx[1]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt2->GetName()))->Fill(dvtx[2]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt3->GetName()))->Fill(dvtx[3]);
    }
    if(Npart==5){
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt4->GetName()))->Fill(pvtx[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt4->GetName()))->Fill(dvtx[4]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt4->GetName()))->Fill(pType[4]);
      ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_V[4]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_V[4]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_V[4]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_V[4]->E());

    }
    if(Npart==6){
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt4->GetName()))->Fill(pvtx[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt4->GetName()))->Fill(dvtx[4]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt4->GetName()))->Fill(pType[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt5->GetName()))->Fill(pvtx[5]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt5->GetName()))->Fill(dvtx[5]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt5->GetName()))->Fill(pType[5]);

      ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_V[4]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_V[4]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_V[4]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_V[4]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_5->GetName()))->Fill( true_V[5]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_5->GetName()))->Fill( true_V[5]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_5->GetName()))->Fill(true_V[5]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_5->GetName()))->Fill( true_V[5]->E());

    }
    if(Npart==7){
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt4->GetName()))->Fill(pvtx[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt4->GetName()))->Fill(pvtx[4]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt4->GetName()))->Fill(pType[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt5->GetName()))->Fill(pvtx[5]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt5->GetName()))->Fill(pvtx[5]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt5->GetName()))->Fill(pType[5]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt6->GetName()))->Fill(pvtx[6]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt6->GetName()))->Fill(pvtx[6]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt6->GetName()))->Fill(pType[6]);

      ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_V[4]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_V[4]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_V[4]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_V[4]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_5->GetName()))->Fill( true_V[5]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_5->GetName()))->Fill( true_V[5]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_5->GetName()))->Fill(true_V[5]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_5->GetName()))->Fill( true_V[5]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_6->GetName()))->Fill( true_V[6]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_6->GetName()))->Fill( true_V[6]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_6->GetName()))->Fill(true_V[6]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_6->GetName()))->Fill( true_V[6]->E());

    }
    if(Npart==8){
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt4->GetName()))->Fill(pvtx[pType[4]]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt4->GetName()))->Fill(dvtx[pType[4]]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt4->GetName()))->Fill(pType[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt5->GetName()))->Fill(pvtx[pType[5]]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt5->GetName()))->Fill(dvtx[pType[5]]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt5->GetName()))->Fill(pType[5]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt6->GetName()))->Fill(pvtx[pType[6]]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt6->GetName()))->Fill( dvtx[pType[6]]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt6->GetName()))->Fill(pType[6]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt7->GetName()))->Fill(pvtx[pType[7]]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt7->GetName()))->Fill( dvtx[pType[7]]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt7->GetName()))->Fill(pType[7]);


      ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_V[4]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_V[4]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_V[4]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_V[4]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_5->GetName()))->Fill( true_V[5]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_5->GetName()))->Fill( true_V[5]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_5->GetName()))->Fill(true_V[5]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_5->GetName()))->Fill( true_V[5]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_6->GetName()))->Fill( true_V[6]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_6->GetName()))->Fill( true_V[6]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_6->GetName()))->Fill(true_V[6]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_6->GetName()))->Fill( true_V[6]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_7->GetName()))->Fill( true_V[7]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_7->GetName()))->Fill( true_V[7]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_7->GetName()))->Fill(true_V[7]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_7->GetName()))->Fill( true_V[7]->E());


    }
    //if(pType[7]!=-1)
    //  ((TH1F*)gDirectory->FindObject(h_mc_vtxdiff_pi_mu2->GetName()))->Fill(pvtx[pType[7]]);
  }

  //epsilon Q
   if(((sevt->trigWord)>>3)&0x1) {
    // ref. trigg. ok:
    ((TH1I*)gDirectory->FindObject(h_CPRE->GetName()))->Fill(evttypebytrks);
    if( (sevt->pu[3].chan[5]>>4)&0x1 ||
	(sevt->pu[4].chan[5]>>4)&0x1 ||//2vtx
	(sevt->pu[5].chan[5]>>4)&0x1 ||  //5 is the central time slice
	(sevt->pu[6].chan[5]>>4)&0x1 ||
	(sevt->pu[3].chan[5]>>8)&0x1 ||//1vtx
	(sevt->pu[4].chan[5]>>8)&0x1 ||
	(sevt->pu[5].chan[5]>>8)&0x1 ||  //5 is the central time slice
	(sevt->pu[6].chan[5]>>8)&0x1 ||
	(sevt->pu[3].chan[5]>>13)&0x1 ||
	(sevt->pu[4].chan[5]>>13)&0x1 ||
	(sevt->pu[5].chan[5]>>13)&0x1 ||  //5 is the central time slice
	(sevt->pu[6].chan[5]>>13)&0x1) {
      if(((sevt->trigWord))&0x1   || ((sevt->trigWord)>>2)&0x1 || ((sevt->trigWord)>>5)&0x1 ) {
	//My trigger also should be OK (constructed from PU)
	((TH1I*)gDirectory->FindObject(h_full_trig->GetName()))->Fill(evttypebytrks);

      }
      //tr_eff = hreal/href;
    }
  }

  //epsilon a
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
      (sevt->pu[6].chan[5]>>13)&0x1
      ) {
    if(((sevt->trigWord))&0x1   || ((sevt->trigWord)>>2)&0x1 || ((sevt->trigWord)>>5)&0x1 ) {
      //My trigger also should be OK (constructed from PU)
      ((TH1I*)gDirectory->FindObject(h_mu1->GetName()))->Fill(evttypebytrks);
      if(((sevt->trigWord)>>3)&0x1) {
	((TH1I*)gDirectory->FindObject(h_Q1->GetName()))->Fill(evttypebytrks);
      }
    }
    //tr_eff = hreal/href;
  }

  ((TH1F*)gDirectory->FindObject( h_vz_initial->GetName()))->Fill( sevt->vtx[0].z);
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

  //Trk hodtimediff
  ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( el1trkhodtime-el2trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( el1trkhodtime-muchtrkhodtime );
  ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( el2trkhodtime-muchtrkhodtime );

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
   ((TH1F*)gDirectory->FindObject(h_gmueeE->GetName()))->Fill( (*mu_nu).E());
  ((TH1F*)gDirectory->FindObject(h_gmueeP->GetName()))->Fill( (*mu_nu).P());
  ((TH1F*)gDirectory->FindObject(h_gmueePt->GetName()))->Fill((*mu_nu).Pt());
  ((TH1F*)gDirectory->FindObject(h_gmueeM->GetName()))->Fill ((*mu_nu).M());

  ((TH1F*)gDirectory->FindObject(h_bx_much->GetName()))            ->Fill(muchp_uncorr[0]);
  ((TH1F*)gDirectory->FindObject(h_by_much->GetName()))            ->Fill(muchp_uncorr[1]);
  ((TH1F*)gDirectory->FindObject(h_bdxdz_much->GetName()))            ->Fill(muchv_uncorr[0]);
  ((TH1F*)gDirectory->FindObject(h_bdydz_much->GetName()))            ->Fill(muchv_uncorr[1]) ;
  ((TH1F*)gDirectory->FindObject(h_bx_vs_by_much->GetName()))        ->Fill(muchp_uncorr[0], muchp_uncorr[1]);
  ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_much->GetName()))->Fill(muchv_uncorr[0], muchv_uncorr[1]);

  ((TH1F*)gDirectory->FindObject(h_k3pi_nu_m->GetName()))->Fill( (*muvee_pi).M() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_nu_P->GetName()))->Fill( (*muvee_pi).P() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_nu_E->GetName()))->Fill( (*muvee_pi).E() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_nu_Pt->GetName()))->Fill( (*muvee_pi).Pt() );
  ((TH1F*)gDirectory->FindObject(h_muee_piee_M->GetName()))->Fill( (*muvee_piee).M() );
  ((TH1F*)gDirectory->FindObject(h_muee_piee_P->GetName()))->Fill( (*muvee_piee).P() );
  ((TH1F*)gDirectory->FindObject(h_muee_piee_E->GetName()))->Fill( (*muvee_piee).E() );
  ((TH1F*)gDirectory->FindObject(h_muee_piee_Pt->GetName()))->Fill( (*muvee_piee).Pt() );

  //((TH1F*)gDirectory->FindObject(h_3trk_invm->GetName()))->Fill ((*k3pic).M2());
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_P->GetName()))->Fill( (*muchlv_pi).P() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_P->GetName()))->Fill( (*el1lv_pi).P() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_P->GetName()))->Fill( (*el2lv_pi).P() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_E->GetName()))->Fill( (*muchlv_pi).E() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_E->GetName()))->Fill( (*el1lv_pi).E() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_E->GetName()))->Fill( (*el2lv_pi).E() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_Pt->GetName()))->Fill((*muchlv_pi).Pt() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_Pt->GetName()))->Fill((*el1lv_pi).Pt() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_Pt->GetName()))->Fill((*el2lv_pi).Pt() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_m->GetName()))->Fill( (*muchlv_pi).M() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_m->GetName()))->Fill( (*el1lv_pi).M() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_m->GetName()))->Fill( (*el2lv_pi).M() );


  // }

  file1->cd();
  h_cutflow->Fill(cutcounter);
  h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());
  cutcounter++;
  ///THE END OF ************* electron E/p between 0.9 and 0.95 (at least one) ****************** ///
  //Mu transverse momentum > 0.9
  /// ************* SIDE BAND ****************** ///

  //  //5ns is according to Manuel's thesis
  dirs[cutcounter-1]->cd();
  if(IS_DATA){
    if (
	fabs(el1trktime - el1trkhodtime)>12.                &&
	fabs(el2trktime - el2trkhodtime)>12.                &&
	fabs(muchtrktime-muchtrkhodtime)>12.                &&
	fabs((long double)(el1trktime-muchtrktime)) > 10.   &&
	fabs((long double)(el2trktime-muchtrktime)) > 10.   &&
	(el1cltime-el2cltime < -3 || el1cltime-el2cltime > 3) &&
	(el1cltime-el2cltime > -33 || el1cltime-el2cltime < 33) &&
     	(el1trkhodtime-el2trkhodtime   < -2 || el1trkhodtime-el2trkhodtime > 2 )   &&
	(el1trkhodtime-el2trkhodtime   > -22 || el1trkhodtime-el2trkhodtime < 22 )   &&
	(el1trkhodtime-muchtrkhodtime  < -2 || el1trkhodtime - muchtrkhodtime > 2) &&
	(el1trkhodtime-muchtrkhodtime  > -22 || el1trkhodtime - muchtrkhodtime < 22) &&
	(el2trkhodtime-muchtrkhodtime  < -2 || el2trkhodtime - muchtrkhodtime > 2)  &&
	(el2trkhodtime-muchtrkhodtime  > -22 || el2trkhodtime - muchtrkhodtime < 22)
	){


      ((TH1F*)gDirectory->FindObject( h_vz_initial->GetName()))->Fill( sevt->vtx[0].z);
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

      //Trk hodtimediff
      ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( el1trkhodtime-el2trkhodtime );
      ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( el1trkhodtime-muchtrkhodtime );
      ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( el2trkhodtime-muchtrkhodtime );


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

      ((TH1F*)gDirectory->FindObject(h_gmueeE->GetName()))->Fill( (*mu_nu).E());
      ((TH1F*)gDirectory->FindObject(h_gmueeP->GetName()))->Fill( (*mu_nu).P());
      ((TH1F*)gDirectory->FindObject(h_gmueePt->GetName()))->Fill((*mu_nu).Pt());
      ((TH1F*)gDirectory->FindObject(h_gmueeM->GetName()))->Fill ((*mu_nu).M());
      ((TH1F*)gDirectory->FindObject(h_bx_much->GetName()))            ->Fill(muchp_uncorr[0]);
      ((TH1F*)gDirectory->FindObject(h_by_much->GetName()))            ->Fill(muchp_uncorr[1]);
      ((TH1F*)gDirectory->FindObject(h_bdxdz_much->GetName()))            ->Fill(muchv_uncorr[0]);
      ((TH1F*)gDirectory->FindObject(h_bdydz_much->GetName()))            ->Fill(muchv_uncorr[1]) ;
      ((TH1F*)gDirectory->FindObject(h_bx_vs_by_much->GetName()))        ->Fill(muchp_uncorr[0], muchp_uncorr[1]);
      ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_much->GetName()))->Fill(muchv_uncorr[0], muchv_uncorr[1]);

      ((TH1F*)gDirectory->FindObject(h_cda_ee->GetName()))->Fill(cda);

      ((TH1F*)gDirectory->FindObject(h_k3pi_nu_m->GetName()))->Fill( (*muvee_pi).M() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_nu_P->GetName()))->Fill( (*muvee_pi).P() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_nu_E->GetName()))->Fill( (*muvee_pi).E() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_nu_Pt->GetName()))->Fill( (*muvee_pi).Pt() );
      ((TH1F*)gDirectory->FindObject(h_muee_piee_M->GetName()))->Fill( (*muvee_piee).M() );
      ((TH1F*)gDirectory->FindObject(h_muee_piee_P->GetName()))->Fill( (*muvee_piee).P() );
      ((TH1F*)gDirectory->FindObject(h_muee_piee_E->GetName()))->Fill( (*muvee_piee).E() );
      ((TH1F*)gDirectory->FindObject(h_muee_piee_Pt->GetName()))->Fill( (*muvee_piee).Pt() );

      //((TH1F*)gDirectory->FindObject(h_3trk_invm->GetName()))->Fill ((*k3pic).M2());
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_P->GetName()))->Fill( (*muchlv_pi).P() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_P->GetName()))->Fill( (*el1lv_pi).P() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_P->GetName()))->Fill( (*el2lv_pi).P() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_E->GetName()))->Fill( (*muchlv_pi).E() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_E->GetName()))->Fill( (*el1lv_pi).E() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_E->GetName()))->Fill( (*el2lv_pi).E() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_Pt->GetName()))->Fill((*muchlv_pi).Pt() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_Pt->GetName()))->Fill((*el1lv_pi).Pt() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_Pt->GetName()))->Fill((*el2lv_pi).Pt() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_m->GetName()))->Fill( (*muchlv_pi).M() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_m->GetName()))->Fill( (*el1lv_pi).M() );
      ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_m->GetName()))->Fill( (*el2lv_pi).M() );

    }
  }

  file1->cd();
  h_cutflow->Fill(cutcounter);
  h_cutflow->GetXaxis()->SetBinLabel( cutcounter, cuts[cutcounter-1].c_str() );
  cutcounter++;
  ///THE END OF ************* SIDE BAND ****************** ///


   /// ************* Time matching ****************** ///
  //"e1e2_hod_Tmatch"
  //  //5ns is according to Manuel's thesis
  if(IS_DATA){
    if( fabs(el1trktime - el1trkhodtime)>12. )     return 0;
    if( fabs(el2trktime - el2trkhodtime)>12. )     return 0;
    if( fabs(muchtrktime-muchtrkhodtime) > 12. )    return 0;  //5 is from Manuel's thesis
    if( fabs((long double)(el1trktime-muchtrktime)) > 10. )    {cleanup();return 0;}
    if( fabs((long double)(el2trktime-muchtrktime)) > 10. )    {cleanup();return 0;}
    if ( el1cltime-el2cltime < -3 || el1cltime-el2cltime > 3 ) {cleanup();return 0;}

    if (  el1trkhodtime-el2trkhodtime  < -2 || el1trkhodtime-el2trkhodtime > 2 ) {cleanup();return 0;}
    if ( el1trkhodtime-muchtrkhodtime  < -2 || el1trkhodtime - muchtrkhodtime > 2 ) {cleanup();return 0;}
    if ( el2trkhodtime-muchtrkhodtime  < -2 || el2trkhodtime - muchtrkhodtime > 2 ) {cleanup();return 0;}


  }
   //
  //}
   //if( fabs((long double)(el1cltime-muchcltime))   > 10. )    {cleanup();return 0;}

  dirs[cutcounter-1]->cd();
  if (IS_MC){
    ////if(true_4.P()  < 3. || true_5.P() < 3.) {cleanup();return 0;}
    ((TH1F*)gDirectory->FindObject(h_mc_Npart_tr->GetName()))->Fill(Npart);
    if (Npart >=4 ){

      ((TH1D*)gDirectory->FindObject( h_mc_muvee_M->GetName()))->Fill( Muee_three_track.M());
      ((TH1D*)gDirectory->FindObject( h_mc_muvee_P->GetName()))->Fill( Muee_three_track.P());
      ((TH1D*)gDirectory->FindObject( h_mc_muvee_Pt->GetName()))->Fill( Muee_three_track.Pt());
      ((TH1D*)gDirectory->FindObject( h_mc_elsP->GetName()))->Fill( mc_ee.P());
      ((TH1D*)gDirectory->FindObject( h_mc_elsM->GetName()))->Fill( mc_ee.M());
      ((TH1D*)gDirectory->FindObject( h_mc_nu_M->GetName()))->Fill( Nu_true.M());
      ((TH1D*)gDirectory->FindObject( h_mc_nu_M2->GetName()))->Fill( Nu_true.M2());
      ((TH1D*)gDirectory->FindObject( h_mc_nu_P->GetName()))->Fill( Nu_true.P());
      ((TH1D*)gDirectory->FindObject( h_mc_nu_Pt->GetName()))->Fill(Nu_true.Pt());
      ((TH1D*)gDirectory->FindObject(h_true_M_0->GetName()))->Fill( true_V[0]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_0->GetName()))->Fill( true_V[0]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_0->GetName()))->Fill(true_V[0]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_0->GetName()))->Fill( true_V[0]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_1->GetName()))->Fill( true_V[1]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_1->GetName()))->Fill( true_V[1]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_1->GetName()))->Fill(true_V[1]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_1->GetName()))->Fill( true_V[1]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_2->GetName()))->Fill( true_V[3]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_2->GetName()))->Fill( true_V[3]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_2->GetName()))->Fill(true_V[3]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_2->GetName()))->Fill( true_V[3]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_3->GetName()))->Fill( true_V[3]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_3->GetName()))->Fill( true_V[3]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_3->GetName()))->Fill(true_V[3]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_3->GetName()))->Fill( true_V[3]->E());

      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt0->GetName()))->Fill(pvtx[0]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt1->GetName()))->Fill(pvtx[1]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt2->GetName()))->Fill(pvtx[2]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt3->GetName()))->Fill(pvtx[3]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt0->GetName()))->Fill(dvtx[0]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt1->GetName()))->Fill(dvtx[1]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt2->GetName()))->Fill(dvtx[2]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt3->GetName()))->Fill(dvtx[3]);
    }
    if(Npart==5){
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt4->GetName()))->Fill(pvtx[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt4->GetName()))->Fill(dvtx[4]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt4->GetName()))->Fill(pType[4]);
      ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_V[4]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_V[4]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_V[4]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_V[4]->E());

    }
    if(Npart==6){
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt4->GetName()))->Fill(pvtx[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt4->GetName()))->Fill(dvtx[4]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt4->GetName()))->Fill(pType[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt5->GetName()))->Fill(pvtx[5]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt5->GetName()))->Fill(dvtx[5]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt5->GetName()))->Fill(pType[5]);

      ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_V[4]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_V[4]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_V[4]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_V[4]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_5->GetName()))->Fill( true_V[5]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_5->GetName()))->Fill( true_V[5]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_5->GetName()))->Fill(true_V[5]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_5->GetName()))->Fill( true_V[5]->E());

    }
    if(Npart==7){
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt4->GetName()))->Fill(pvtx[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt4->GetName()))->Fill(pvtx[4]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt4->GetName()))->Fill(pType[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt5->GetName()))->Fill(pvtx[5]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt5->GetName()))->Fill(pvtx[5]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt5->GetName()))->Fill(pType[5]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt6->GetName()))->Fill(pvtx[6]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt6->GetName()))->Fill(pvtx[6]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt6->GetName()))->Fill(pType[6]);

      ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_V[4]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_V[4]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_V[4]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_V[4]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_5->GetName()))->Fill( true_V[5]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_5->GetName()))->Fill( true_V[5]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_5->GetName()))->Fill(true_V[5]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_5->GetName()))->Fill( true_V[5]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_6->GetName()))->Fill( true_V[6]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_6->GetName()))->Fill( true_V[6]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_6->GetName()))->Fill(true_V[6]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_6->GetName()))->Fill( true_V[6]->E());

    }
    if(Npart==8){
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt4->GetName()))->Fill(pvtx[pType[4]]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt4->GetName()))->Fill(dvtx[pType[4]]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt4->GetName()))->Fill(pType[4]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt5->GetName()))->Fill(pvtx[pType[5]]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt5->GetName()))->Fill(dvtx[pType[5]]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt5->GetName()))->Fill(pType[5]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt6->GetName()))->Fill(pvtx[pType[6]]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt6->GetName()))->Fill( dvtx[pType[6]]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt6->GetName()))->Fill(pType[6]);
      ((TH1F*)gDirectory->FindObject(h_mc_pvtx_pt7->GetName()))->Fill(pvtx[pType[7]]);
      ((TH1F*)gDirectory->FindObject(h_mc_dvtx_pt7->GetName()))->Fill( dvtx[pType[7]]);
      ((TH1I*)gDirectory->FindObject(h_mc_pt7->GetName()))->Fill(pType[7]);


      ((TH1D*)gDirectory->FindObject(h_true_M_4->GetName()))->Fill( true_V[4]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_4->GetName()))->Fill( true_V[4]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_4->GetName()))->Fill(true_V[4]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_4->GetName()))->Fill( true_V[4]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_5->GetName()))->Fill( true_V[5]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_5->GetName()))->Fill( true_V[5]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_5->GetName()))->Fill(true_V[5]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_5->GetName()))->Fill( true_V[5]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_6->GetName()))->Fill( true_V[6]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_6->GetName()))->Fill( true_V[6]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_6->GetName()))->Fill(true_V[6]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_6->GetName()))->Fill( true_V[6]->E());
      ((TH1D*)gDirectory->FindObject(h_true_M_7->GetName()))->Fill( true_V[7]->M());
      ((TH1D*)gDirectory->FindObject(h_true_P_7->GetName()))->Fill( true_V[7]->P());
      ((TH1D*)gDirectory->FindObject(h_true_Pt_7->GetName()))->Fill(true_V[7]->Pt());
      ((TH1D*)gDirectory->FindObject(h_true_E_7->GetName()))->Fill( true_V[7]->E());


    }
    //if(pType[7]!=-1)
    //  ((TH1F*)gDirectory->FindObject(h_mc_vtxdiff_pi_mu2->GetName()))->Fill(pvtx[pType[7]]);
  }


  if(((sevt->trigWord)>>3)&0x1) {
    // ref. trigg. ok:
    ((TH1I*)gDirectory->FindObject(h_CPRE->GetName()))->Fill(evttypebytrks);
    if( (sevt->pu[3].chan[5]>>4)&0x1 ||
	(sevt->pu[4].chan[5]>>4)&0x1 ||
	(sevt->pu[5].chan[5]>>4)&0x1 ||  //5 is the central time slice
	(sevt->pu[6].chan[5]>>4)&0x1 ||
	(sevt->pu[3].chan[5]>>8)&0x1 ||
	(sevt->pu[4].chan[5]>>8)&0x1 ||
	(sevt->pu[5].chan[5]>>8)&0x1 ||  //5 is the central time slice
	(sevt->pu[6].chan[5]>>8)&0x1 ) {
      if(((sevt->trigWord))&0x1   || ((sevt->trigWord)>>2)&0x1  ) {
	//My trigger also should be OK (constructed from PU)
	((TH1I*)gDirectory->FindObject(h_full_trig->GetName()))->Fill(evttypebytrks);
      }
      //tr_eff = hreal/href;
    }
  }

  if((sevt->pu[3].chan[5]>>8)&0x1 ||
     (sevt->pu[4].chan[5]>>8)&0x1 ||   //1vtx
     (sevt->pu[5].chan[5]>>8)&0x1 ||  //5 is the central time slice
     (sevt->pu[6].chan[5]>>8)&0x1 ) {
    if( ((sevt->trigWord)>>2)&0x1  ) {
      //My trigger also should be OK (constructed from PU)
      ((TH1I*)gDirectory->FindObject(h_1vtx_mcut->GetName()))->Fill(evttypebytrks);
    }
    //tr_eff = hreal/href;
  }
  if( (sevt->pu[3].chan[5]>>4)&0x1 ||
      (sevt->pu[4].chan[5]>>4)&0x1 ||//2vtx
      (sevt->pu[5].chan[5]>>4)&0x1 ||  //5 is the central time slice
      (sevt->pu[6].chan[5]>>4)&0x1 ) {
    if(((sevt->trigWord))&0x1 ) {
      //My trigger also should be OK (constructed from PU)
      ((TH1I*)gDirectory->FindObject(h_2vtx->GetName()))->Fill(evttypebytrks);
    }
    //tr_eff = hreal/href;
  }
  if( (sevt->pu[3].chan[5]>>4)&0x1 ||
      (sevt->pu[4].chan[5]>>4)&0x1 ||
      (sevt->pu[5].chan[5]>>4)&0x1 ||  //5 is the central time slice
      (sevt->pu[6].chan[5]>>4)&0x1 ||
      (sevt->pu[3].chan[5]>>8)&0x1 ||
      (sevt->pu[4].chan[5]>>8)&0x1 ||
      (sevt->pu[5].chan[5]>>8)&0x1 ||  //5 is the central time slice
      (sevt->pu[6].chan[5]>>8)&0x1 ) {
    if(((sevt->trigWord))&0x1   || ((sevt->trigWord)>>2)&0x1  ) {
      //My trigger also should be OK (constructed from PU)
      ((TH1I*)gDirectory->FindObject(h_mu1->GetName()))->Fill(evttypebytrks);
    }
    //tr_eff = hreal/href;
  }
  if( (sevt->pu[3].chan[5]>>13)&0x1 ||
      (sevt->pu[4].chan[5]>>13)&0x1 ||
      (sevt->pu[5].chan[5]>>13)&0x1 ||  //5 is the central time slice
      (sevt->pu[6].chan[5]>>13)&0x1 ) {
    if( ((sevt->trigWord)>>5)&0x1 ) {
      //My trigger also should be OK (constructed from PU)
      ((TH1I*)gDirectory->FindObject(h_Q1_NAKL->GetName()))->Fill(evttypebytrks);
    }
    //tr_eff = hreal/href;
  }



  ((TH1F*)gDirectory->FindObject( h_vz_initial->GetName()))->Fill( sevt->vtx[0].z);
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

  //Trk hodtimediff
  ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( el1trkhodtime-el2trkhodtime );
  ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( el1trkhodtime-muchtrkhodtime );
  ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( el2trkhodtime-muchtrkhodtime );


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

   ((TH1F*)gDirectory->FindObject(h_gmueeE->GetName()))->Fill( (*mu_nu).E());
  ((TH1F*)gDirectory->FindObject(h_gmueeP->GetName()))->Fill( (*mu_nu).P());
  ((TH1F*)gDirectory->FindObject(h_gmueePt->GetName()))->Fill((*mu_nu).Pt());
  ((TH1F*)gDirectory->FindObject(h_gmueeM->GetName()))->Fill ((*mu_nu).M());
  ((TH1F*)gDirectory->FindObject(h_bx_much->GetName()))            ->Fill(muchp_uncorr[0]);
  ((TH1F*)gDirectory->FindObject(h_by_much->GetName()))            ->Fill(muchp_uncorr[1]);
  ((TH1F*)gDirectory->FindObject(h_bdxdz_much->GetName()))            ->Fill(muchv_uncorr[0]);
  ((TH1F*)gDirectory->FindObject(h_bdydz_much->GetName()))            ->Fill(muchv_uncorr[1]) ;
  ((TH1F*)gDirectory->FindObject(h_bx_vs_by_much->GetName()))        ->Fill(muchp_uncorr[0], muchp_uncorr[1]);
  ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_much->GetName()))->Fill(muchv_uncorr[0], muchv_uncorr[1]);

  ((TH1F*)gDirectory->FindObject(h_cda_ee->GetName()))->Fill(cda);

   ((TH1F*)gDirectory->FindObject(h_k3pi_nu_m->GetName()))->Fill( (*muvee_pi).M() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_nu_P->GetName()))->Fill( (*muvee_pi).P() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_nu_E->GetName()))->Fill( (*muvee_pi).E() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_nu_Pt->GetName()))->Fill( (*muvee_pi).Pt() );
   ((TH1F*)gDirectory->FindObject(h_muee_piee_M->GetName()))->Fill( (*muvee_piee).M() );
  ((TH1F*)gDirectory->FindObject(h_muee_piee_P->GetName()))->Fill( (*muvee_piee).P() );
  ((TH1F*)gDirectory->FindObject(h_muee_piee_E->GetName()))->Fill( (*muvee_piee).E() );
  ((TH1F*)gDirectory->FindObject(h_muee_piee_Pt->GetName()))->Fill( (*muvee_piee).Pt() );

  //((TH1F*)gDirectory->FindObject(h_3trk_invm->GetName()))->Fill ((*k3pic).M2());
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_P->GetName()))->Fill( (*muchlv_pi).P() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_P->GetName()))->Fill( (*el1lv_pi).P() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_P->GetName()))->Fill( (*el2lv_pi).P() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_E->GetName()))->Fill( (*muchlv_pi).E() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_E->GetName()))->Fill( (*el1lv_pi).E() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_E->GetName()))->Fill( (*el2lv_pi).E() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_Pt->GetName()))->Fill((*muchlv_pi).Pt() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_Pt->GetName()))->Fill((*el1lv_pi).Pt() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_Pt->GetName()))->Fill((*el2lv_pi).Pt() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi1_m->GetName()))->Fill( (*muchlv_pi).M() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi2_m->GetName()))->Fill( (*el1lv_pi).M() );
  ((TH1F*)gDirectory->FindObject(h_k3pi_pi3_m->GetName()))->Fill( (*el2lv_pi).M() );



  file1->cd();
  h_cutflow->Fill(cutcounter);
  h_cutflow->GetXaxis()->SetBinLabel( cutcounter, cuts[cutcounter-1].c_str() );
  cutcounter++;
  ///THE END OF ************* Time matching ****************** ///



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

  if(el1lv)         {    delete el1lv;     el1lv =NULL;}
  if(el2lv)         {    delete el2lv;     el2lv =NULL;}
  if(el1lv_pi)      {       delete el1lv_pi;     el1lv_pi =NULL;}
  if(el2lv_pi)      {       delete el2lv_pi;     el2lv_pi =NULL;}
  if(k3pi_miss_mass)  {   delete k3pi_miss_mass; k3pi_miss_mass =NULL;}
  if(pi1lv)      {      delete pi1lv;        pi1lv=NULL;}
  if(pi2lv)      {      delete pi2lv; pi2lv=NULL;	}
  if(pi3lv)      {      delete pi3lv; pi3lv=NULL;	}
  if(k3pic)      {       delete k3pic; k3pic =NULL;	}
  if(pieelv)     {        delete pieelv; pieelv =NULL;	}
  if(eelv)       {      delete eelv;      eelv =NULL;	}
  if(eelv_pi)    { delete eelv_pi; eelv_pi =NULL;	}
  if(muchlv)     { delete muchlv;     muchlv =NULL;	}
  if(pchlv)      { delete pchlv;     pchlv =NULL;	}
  if(muchlv_pi)  { delete muchlv_pi; muchlv_pi =NULL;	}
  if(muvee)      { delete muvee;    muvee =NULL;	}
  if(pK)         { delete pK; pK=NULL;			}
  if(neutrino)   { delete neutrino; neutrino=NULL;	}
  if(muvee_piee) { delete muvee_piee; muvee_piee=NULL;	}
  if(muvee_pi)   {     delete muvee_pi; muvee_pi=NULL;	}
  if(v_second_e1){ delete v_second_e1;v_second_e1 =NULL;}
  if(v_second_e2){ delete v_second_e2;v_second_e2 =NULL;}
  if(v_second_mu){ delete v_second_mu;v_second_mu =NULL;}
  if(v_first_e1) {delete v_first_e1;v_first_e1 =NULL;	}
  if(v_first_e2) {delete v_first_e2;v_first_e2 =NULL;	}
  if(v_first_mu) {delete v_first_mu;v_first_mu =NULL;   }
  if(mu_nu) {delete mu_nu;mu_nu =NULL;   }
  if(K4Momentum){ delete K4Momentum ;K4Momentum =NULL;}
  if(Kmomentum) {delete Kmomentum ;Kmomentum =NULL;   }

  if(glv)      {      delete glv;       glv =NULL;                }
  if(Klv)      {      delete Klv;       Klv =NULL;		  }
  if(el1lvcm)  {      delete el1lvcm;   el1lvcm =NULL;		  }
  if(el2lvcm)  {      delete el2lvcm;   el2lvcm =NULL;		  }
  if(glvcm)    {        delete glvcm;     glvcm =NULL;		  }
  if(pi0lvcm)  {      delete pi0lvcm;   pi0lvcm =NULL;		  }
  if(pi0lv)    {        delete pi0lv;     pi0lv =NULL;		  }
   if(muveeKcm){            delete muveeKcm;     muveeKcm =NULL;  }
  if(pi0lvKcm) {       delete pi0lvKcm;  pi0lvKcm =NULL;	  }
  if(pchlvKcm) {       delete pchlvKcm;  pchlvKcm =NULL;	  }
  if(muchlvKcm){        delete muchlvKcm;  muchlvKcm =NULL;	  }
  if(KlvKcm)   {         delete KlvKcm;    KlvKcm =NULL;          }
  //if(K_true)            delete K_true;     K_true=0;
  //if(Pi0_true)            delete Pi0_true;     Pi0_true=0;
  //if(El_true)            delete El_true;     El_true=0;
  //if(Nu_true)            delete Nu_true;     Nu_true=0;
  //if(Muee_three_track) delete Muee_three_track;Muee_three_track =0;
  //if(mc_ee)            delete mc_ee;     mc_ee=0;
  //if(true_0)            delete true_0 ;true_0=0;
  //if(true_1)            delete true_1; true_1=0;
  //if(true_2)            delete true_2; true_2=0;
  //if(true_3)            delete true_3; true_3=0;
  //if(true_4)            delete true_4; true_4=0;
  //if(true_5)            delete true_5; true_5=0;
  //if(true_6)            delete true_6; true_6=0;
  //if(true_7)            delete true_7; true_7=0;

  for (unsigned int j=0 ; j<=7 ; j++) {
    if(true_V[j])            delete true_V[j];  true_V[j] =0;
  }
  //if(trkCorrSlope)     delete trkCorrSlope;     trkCorrSlope=NULL;
  //if(trkCorrMidPoint) delete trkCorrMidPoint; trkCorrMidPoint=NULL;
  /*
  for (unsigned int i = trkCorrSlopesV.size()-1; i>-1; i--) {
    if (trkCorrSlopesV[i])         delete trkCorrSlopesV[i];
    //  if (trkCorrMidPointsV[i])     delete trkCorrMidPointsV[i];

  }

  for (unsigned int i = trkCorrMidPointsV.size()-1; i>-1; i--) {
    if(trkCorrMidPointsV[i])     delete trkCorrMidPointsV[i];
  }

  trkCorrSlopesV.clear();
  trkCorrMidPointsV.clear();
  */
  /*
  for (unsigned int i = clcorrpos.size()-1; i>-1; i--) {
    if (clcorrpos[i])         delete clcorrpos[i];
  }
  clcorrpos.clear();
  */
  //for (i = clcorrpos.begin(); i != clcorrpos.end(); ++i) {
  //  delete (*i);
  //  *i=0;
  //}
  //clcorrpos.clear();

 // for (i = trkCorrMidPointsV.begin(); i != trkCorrMidPointsV.end(); ++i) {
 //   delete (*i);
 //   *i=0;
 // }
 // trkCorrMidPointsV.clear();
  //printf("CLCORRPOS:  %d , %d, %d \n" ,  clcorrpos.empty(),clcorrpos.size(),clcorrpos.capacity());
  //for(int i = 0; i< clcorrpos.size();i++) {
  //  printf("%d\n" , clcorrpos.at(i) );
  //  if(clcorrpos.at(i)) clcorrpos.at(i)->Print();
  //}
  while(!trkCorrMidPointsV.empty() ){delete trkCorrMidPointsV.back(); trkCorrMidPointsV.pop_back();}
  while(!trkCorrSlopesV.empty() ){delete trkCorrSlopesV.back(); trkCorrSlopesV.pop_back();}
  while(!clcorrpos.empty() ){if(clcorrpos.back())   delete clcorrpos.back(); clcorrpos.pop_back();}
  //  if(clpos)             delete clpos;             clpos=NULL;

  // while(!trkCorrMidPointsV.empty() ){delete trkCorrMidPointsV.back(); trkCorrMidPointsV.pop_back();}

  return;
}
