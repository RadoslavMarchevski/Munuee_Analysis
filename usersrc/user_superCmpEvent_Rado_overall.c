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
TLorentzVector *pchlv=NULL, *muchlv=NULL , *muvee=NULL;
TLorentzVector *pKaon_simple=NULL, *muvee_neutrino=NULL;
TLorentzVector *pieelv=NULL;
TLorentzVector *glv =NULL;
TLorentzVector *pi0lv=NULL;
TLorentzVector *el1lvcm=NULL, *el2lvcm=NULL,  *glvcm=NULL, *pi0lvcm=NULL;
TLorentzVector *Klv=NULL, *pi0lvKcm=NULL, *pchlvKcm=NULL, *muchlvKcm=NULL, *KlvKcm=NULL,*muveeKcm=NULL;


void cleanup();

void applyBlueTubeCorrection(superCmpEvent *sevt, int itrk, double pKaon[3], double vKaon[3], double* trkSlopes, double* trkMidPoints);

int user_superCmpEvent(superBurst *sbur,superCmpEvent *sevt) {
    /* WARNING: do not alter things before this line */
    /*---------- Add user C code here ----------*/
    int nuserevt=0;
    int  a, i, j, k;
    int  l, m, n;

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
    double el1trkhodtime=-1e+3;
    double el2trkhodtime=-1e+3;
    double el1trktime=-1e+3;
    double el2trktime=-1e+3;

    double norm1, norm2,normmu;
    double p1x,p1y, p1z, p2x,p2y, p2z, el1E, el2E, pK1x,pK1y, pK1z,pKE,pmux,pmuy,pmuz;

    int    gclindx=-1;
    double gclenergy=0.;
    double gcltime=0.;
    double gclx=0.;
    double gcly=0.;
    double gclz=0.;
    double gtrkDHC1x, gtrkDHC1y, gtrkRadDHC1 ;


    double pchtrkhodtime, pchE , muchE , muchtrkhodtime;
    
    double eldistDCH1plane, el1xLKrplane, el1yLKrplane, el2xLKrplane, el2yLKrplane, eldistLKrplane;
    double gel1distLKrplane, gel2distLKrplane;
    double pche1distDHC1plane, pche2distDHC1plane, pche1distLKrplane, pche2distLKrplane;
    double muche1distDHC1plane, muche2distDHC1plane, muche1distLKrplane, muche2distLKrplane;
    
    double pchxLKrplane, pchyLKrplane, gpchdistLKrplane;
    double muchxLKrplane, muchyLKrplane, gmuchdistLKrplane;
    double Dgvtx_x,Dgvtx_y, Dgvtx_z, gnorm,  gpx, gpy, gpz;

    unsigned int vertex_OK ,pchee_vertex_OK ;
    double vertex[3], vertex_e2pch[3], vertex_e1pch[3];
    double cda, cda_e2pch, cda_e1pch;

    unsigned int muchee_vertex_OK;
    double vertex_e2much[3], vertex_e1much[3];
    double cda_e2much, cda_e1much;

    double elp1_uncorr[3], elp1[3], elv1_uncorr[3], elv1[3], elp2_uncorr[3],elp2[3], elv2_uncorr[3],elv2[3], pchp_uncorr[3],pchp[3], pchv_uncorr[3], pchv[3] ,muchp_uncorr[3],muchp[3], muchv_uncorr[3], muchv[3];

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
    double eleop1=0;
    double eleop2=0;
    double el1cltime=-1e+3;
    double el2cltime=-1e+3;
    double el1P3;
    double el2P3;

    double eelvmass;

    double c_eop_e=0.90;//0.95;
    double c_eop_e_up=1.05;
    double c_eop_mu_up=0.1;
    double c_eop_pi_min=0.2;
    double c_eop_pi_max=0.85;


    bool   LKrCalCorr=true;
    bool   CorrAlphaBeta=true;

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
    //"Overall cuts "
    dirs[cutcounter-1]->cd();

    h_ntrack->Fill(ntrack);
    h_nclust->Fill(nclust);

    //file1->cd();
    //h_cutflow->Fill(cutcounter);
    //h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"no_cuts");
    //cutcounter++;
    


    /// ************* tracks  &  clusters  *******************************

    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~

    if(ntrack!=3 ) {cleanup();return 0;}
    //"3tracks"
    //111
    //dirs[cutcounter-1]->cd();

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
        if(icl_track>-1)
            clst_track     = sevt->cluster[icl_track].status;

        //        if(imu_track>-1)
        //            cout<<i<<"    imu_track="<<imu_track<<"    "<<sevt->muon[imu_track].chi2<<endl;

        if (quality_track<0.8)        { trkCorrSlopesV.push_back(NULL); trkCorrMidPointsV.push_back(NULL); continue; } //0.8 is in Manuel's thesis
        if (icl_track<0)             { trkCorrSlopesV.push_back(NULL); trkCorrMidPointsV.push_back(NULL); continue; }
        //if (imu_track>-1)             { trkCorrSlopesV.push_back(NULL); trkCorrMidPointsV.push_back(NULL); continue; }
        if (ddcell_track<2.)        { trkCorrSlopesV.push_back(NULL); trkCorrMidPointsV.push_back(NULL); continue; }
        if (p_track<3.)                { trkCorrSlopesV.push_back(NULL); trkCorrMidPointsV.push_back(NULL); continue; }
        if (p_track>75.)            { trkCorrSlopesV.push_back(NULL); trkCorrMidPointsV.push_back(NULL); continue; }
        if (clst_track>3)            { trkCorrSlopesV.push_back(NULL); trkCorrMidPointsV.push_back(NULL); continue; }


        //calculate distance between the track and its matching cluster in LKr
        gclx=clcorrpos[icl_track]->x();
        gcly=clcorrpos[icl_track]->y();
        gclz=clcorrpos[icl_track]->z();
        trkxcl = sevt->track[i].x + (gclz-DCHz) * sevt->track[i].dxdz ;
        trkycl = sevt->track[i].y + (gclz-DCHz) * sevt->track[i].dydz ;
        dtrkcl = sqrt( pow(trkxcl-gclx, 2.) + pow(trkycl-gcly, 2.) );

        ((TH1F*)gDirectory->FindObject(h_dtrkcl->GetName()))        ->Fill( dtrkcl );

        if(dtrkcl>6.)                 { trkCorrSlopesV.push_back(NULL); trkCorrMidPointsV.push_back(NULL); continue; }

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
    }

    //        cout<<"1"<<endl;

    //file1->cd();
    //1112
    //h_cutflow->Fill(cutcounter);
    //h_cutflow->GetXaxis()->SetBinLabel( cutcounter, cuts[cutcounter-1].c_str() );// "3tracks_4clusters");
    //cutcounter++;




    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~

    if(ngoodtrack!=3) {cleanup();return 0;}
    //"3goodtracks"
    //dirs[cutcounter-1]->cd();


    dtrkcl=249.;
    for (int i=0; i<nclust; i++) {
        if(sevt->cluster[i].iTrack <0)     {
                gclindx = i;
                gclenergy = sevt->cluster[i].energy;
                gcltime = sevt->cluster[i].time;
                //                gclx = sevt->cluster[i].x;
                //                gcly = sevt->cluster[i].y;
                gclx=clcorrpos[i]->x();
                gcly=clcorrpos[i]->y();
                gclz=clcorrpos[i]->z();

                float dtrkcl_var=0.;
                for (int itrk = 0; itrk < ntrack; ++itrk) {
                    trkxcl         = sevt->track[itrk].x + (gclz-DCHz) * sevt->track[itrk].dxdz ;
                    trkycl         = sevt->track[itrk].y + (gclz-DCHz) * sevt->track[itrk].dydz ;
                    dtrkcl_var     = sqrt( pow(trkxcl-gclx, 2.) + pow(trkycl-gcly, 2.) );
                    if(dtrkcl>dtrkcl_var)
                        dtrkcl=dtrkcl_var;
                }

                break;
        }
    }


    //    cout<< "    "<<endl;
    for (int i=0; i<ntrack; i++) {

        p_track            = sevt->track[i].p;
        icl_track         = sevt->track[i].iClus;

        //TODO: ask Zhenja about this stuff. Do I need them for data 2003/2004 and where to get from?
        //cout<<"EopCorr[periodFlag][CPDindex][CELLindex] = "<<EopCorr[periodFlag][CPDindex][CELLindex]<<endl;
        //    if (IS_DATA)                                 // ab v65 Korrekturen automatisch separat f�r jede (Sub-) Periode
        //          e_track = sevt->cluster[icl_track].energy / EopCorr[periodFlag][CPDindex][CELLindex];  // Ke3 E/p-Korrektur f�r jede Zelle!!!
        //    else
        e_track  = sevt->cluster[icl_track].energy;
	
        eop_track  = e_track / p_track;
        ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_track);
        ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_track, p_track);
	
        //                cout<< "E/p="<<eop_track<<"    ";
        //if(eop_track>c_eop_e)    {

	if(eop_track>c_eop_e  &&  eop_track<c_eop_e_up)    {
	  
	  //            cout<< "electron ";
	  //WORK ON  making e+ , e-!!!!!!
	  if(eltrkindx1<0)    {
	    eltrkindx1=i;
	    eltrkq1 = sevt->track[i].q;
	  }
	  else if(eltrkindx2<0)    {
	    eltrkindx2=i;
	    eltrkq2 = sevt->track[i].q;
	  }
	  else { ; }
	  
	  nelectrack++;
        }
        // muon track
        else if(eop_track<c_eop_mu_up)    {
	  if(muchtrkindx < 0 && imu_track !=-1)    {
	    muchtrkindx=i;
	    muchtrkq=sevt->track[i].q;
	  }
	    
            nmutrack++;
        }
        //track with too high e/p. skip it.
        else    {
	  ;//            cout<< "?? ";;
        }
        //        cout<< " "<<endl;
    }
	
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
    //K3pi track assignment selection
    /* if(nelectrack+nmutrack==3)    {
        if(muchtrkq>0.)    {
            if(eltrkq1 != eltrkq2)    {// pi+, pi+, pi-
                evttypebytrks=0;
            }
            else if(eltrkq1>0)    {    // pi+, pi+, pi+
                evttypebytrks=1;
            }
            else    {                // pi+, pi-, pi-
                evttypebytrks=2;
            }
        }
        else    {
            if(eltrkq1 != eltrkq2)    {// pi-, pi+, pi-
                evttypebytrks=3;
            }
            else if(eltrkq1>0)    {    // pi-, pi+, pi+
                evttypebytrks=4;
            }
            else    {                // pi-, pi-, pi-
                evttypebytrks=5;
            }
        }
	}
    */
     
    ((TH1I*)gDirectory->FindObject(h_evttypebytrks->GetName()))->Fill(evttypebytrks);
    ((TH1I*)gDirectory->FindObject(h_nelectrack->GetName()))->Fill(nelectrack);
    ((TH1I*)gDirectory->FindObject(h_mutrack->GetName()))->Fill(nmutrack);
    ((TH2I*)gDirectory->FindObject(h_neltrk_mutrk->GetName()))->Fill(nelectrack, nmutrack);

    //((TH1F*)gDirectory->FindObject(h_dtrkcl->GetName()))     ->Fill( dtrkcl );

    //file1->cd();
    //h_cutflow->Fill(cutcounter);
    //h_cutflow->GetXaxis()->SetBinLabel( cutcounter, cuts[cutcounter-1].c_str() );// "3goodtracks");
    //cutcounter++;





    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~

    if(dtrkcl<20.)        {cleanup();return 0;}
    //"cluster_isolation"
    //dirs[cutcounter-1]->cd();
    
    //((TH1I*)gDirectory->FindObject(h_evttypebytrks->GetName()))->Fill(evttypebytrks);
    //((TH1I*)gDirectory->FindObject(h_nelectrack->GetName()))->Fill(nelectrack);
    //((TH1I*)gDirectory->FindObject(h_mutrack->GetName()))->Fill(nmutrack);
    //((TH2I*)gDirectory->FindObject(h_neltrk_mutrk->GetName()))->Fill(nelectrack, nmutrack);

    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel( cutcounter, cuts[cutcounter-1].c_str() );// "Initial cuts" // "cluster_isolation");
    cutcounter++;


    
    
    
    /// ************* electron / positron

    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~

    if (nelectrack!=2)        {cleanup();return 0;}
    if (eltrkq1==eltrkq2)    {cleanup();return 0;}
    //"electron_positron"
    dirs[cutcounter-1]->cd();

    el1trktime      = sevt->track[eltrkindx1].time ;
    el1trkhodtime    = sevt->track[eltrkindx1].hodTime ;
    icl_track       = sevt->track[eltrkindx1].iClus;
    el1cltime        = sevt->cluster[icl_track].time;

    ((TH1I*)gDirectory->FindObject(h_evttypebytrks->GetName()))->Fill(evttypebytrks);
    ((TH1I*)gDirectory->FindObject(h_nelectrack->GetName()))->Fill(nelectrack);
    ((TH1I*)gDirectory->FindObject(h_mutrack->GetName()))->Fill(nmutrack);
    ((TH1I*)gDirectory->FindObject(h_mutrackq->GetName()))->Fill(muchtrkq);
    ((TH2I*)gDirectory->FindObject(h_neltrk_mutrk->GetName()))->Fill(nelectrack, nmutrack);


    ((TH1F*)gDirectory->FindObject(h_el1trktime->GetName()))        ->Fill( el1trktime );
    ((TH1F*)gDirectory->FindObject(h_el1trkhodtime->GetName()))        ->Fill( el1trkhodtime );
    ((TH1F*)gDirectory->FindObject(h_el1trktimediff->GetName()))    ->Fill( el1trktime - el1trkhodtime );
    ((TH1F*)gDirectory->FindObject(h_ecltime1->GetName()))            ->Fill( el1cltime );
    ((TH1F*)gDirectory->FindObject(h_e1trkcltimediff->GetName()))    ->Fill( el1trktime- el1cltime );

    //    cout<< "    first "<<eltrkindx1<<"    icl_track="<<icl_track<<"    Dt="<<el1trktime - el1trkhodtime<<endl;

    el2trktime      = sevt->track[eltrkindx2].time ;
    el2trkhodtime    = sevt->track[eltrkindx2].hodTime ;
    icl_track       = sevt->track[eltrkindx2].iClus;
    el2cltime        = sevt->cluster[icl_track].time;
    ((TH1F*)gDirectory->FindObject(h_el2trktime->GetName()))        ->Fill( el2trktime );
    ((TH1F*)gDirectory->FindObject(h_el2trkhodtime->GetName()))        ->Fill( el2trkhodtime );
    ((TH1F*)gDirectory->FindObject(h_el2trktimediff->GetName()))    ->Fill( el2trktime - el2trkhodtime );
    ((TH1F*)gDirectory->FindObject(h_ecltime2->GetName()))            ->Fill( el2cltime );
    ((TH1F*)gDirectory->FindObject(h_e2trkcltimediff->GetName()))    ->Fill( el2trktime-el2cltime );

    

       ((TH1F*)gDirectory->FindObject(h_eop_mu->GetName()))->Fill(sevt->cluster[sevt->track[muchtrkindx].iClus].energy/sevt->track[muchtrkindx].p);
    ((TH2F*)gDirectory->FindObject(h_eop_mu_vs_p->GetName()))->Fill(sevt->cluster[sevt->track[muchtrkindx].iClus].energy/sevt->track[muchtrkindx].p, sevt->track[muchtrkindx].p);

     ((TH1F*)gDirectory->FindObject(h_eop_el->GetName()))->Fill(sevt->cluster[sevt->track[eltrkindx1].iClus].energy/sevt->track[eltrkindx1].p);
    ((TH2F*)gDirectory->FindObject(h_eop_el_vs_p->GetName()))->Fill(sevt->cluster[sevt->track[eltrkindx1].iClus].energy/sevt->track[eltrkindx1].p, sevt->track[eltrkindx1].p);

     ((TH1F*)gDirectory->FindObject(h_eop_el->GetName()))->Fill(sevt->cluster[sevt->track[eltrkindx2].iClus].energy/sevt->track[eltrkindx2].p);
    ((TH2F*)gDirectory->FindObject(h_eop_el_vs_p->GetName()))->Fill(sevt->cluster[sevt->track[eltrkindx2].iClus].energy/sevt->track[eltrkindx2].p, sevt->track[eltrkindx2].p);

    
    //    cout<< "    second "<<eltrkindx2<<"    icl_track="<<icl_track<<"    Dt="<<el2trktime - el2trkhodtime<<endl;

    //e1e2
    ((TH1F*)gDirectory->FindObject(h_e1e2timediff->GetName()))        ->Fill( el1trktime-el2trktime );
    ((TH1F*)gDirectory->FindObject(h_e1e2cltimediff->GetName()))    ->Fill( el1cltime-el2cltime );
    
    ((TH1F*)gDirectory->FindObject(h_e1e2hodtimediff->GetName()))        ->Fill( el1trkhodtime-el2trkhodtime );
    
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


    //blue field corrected values
    elv1[0] = trkCorrSlopesV[eltrkindx1]->x();
    elv1[1] = trkCorrSlopesV[eltrkindx1]->y();
    elv1[2] = 1.;
    elp1[0] = trkCorrMidPointsV[eltrkindx1]->x();
    elp1[1] = trkCorrMidPointsV[eltrkindx1]->y();
    elp1[2] = trkCorrMidPointsV[eltrkindx1]->z();

    elv2[0] = trkCorrSlopesV[eltrkindx2]->x();
    elv2[1] = trkCorrSlopesV[eltrkindx2]->y();
    elv2[2] = 1.;
    elp2[0] = trkCorrMidPointsV[eltrkindx2]->x();
    elp2[1] = trkCorrMidPointsV[eltrkindx2]->y();
    elp2[2] = trkCorrMidPointsV[eltrkindx2]->z();





    ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill( elp1_uncorr[0] );
    ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill( elp1_uncorr[1] );
    ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill( elv1_uncorr[0] );
    ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill( elv1_uncorr[1] );
    ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill( elp1_uncorr[0], elp1_uncorr[1] );
    ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill( elv1_uncorr[0], elv1_uncorr[1] );

    //is this correct? should not one use uncorrected values? find it out.
    // answer: yes, should be uncorrected ones
    norm1 = 1./sqrt(elv1_uncorr[0]*elv1_uncorr[0] + elv1_uncorr[1]*elv1_uncorr[1] + elv1_uncorr[2]*elv1_uncorr[2] );
    norm2 = 1./sqrt(elv2_uncorr[0]*elv2_uncorr[0] + elv2_uncorr[1]*elv2_uncorr[1] + elv2_uncorr[2]*elv2_uncorr[2] );
    p1x = norm1 * elv1_uncorr[0] * sevt->track[eltrkindx1].p ;
    p1y = norm1 * elv1_uncorr[1] * sevt->track[eltrkindx1].p ;
    p1z = norm1 *                  sevt->track[eltrkindx1].p ;
    p2x = norm2 * elv2_uncorr[0] * sevt->track[eltrkindx2].p ;
    p2y = norm2 * elv2_uncorr[1] * sevt->track[eltrkindx2].p ;
    p2z = norm2 *                  sevt->track[eltrkindx2].p ;
    el1E = sqrt( p1x*p1x+p1y*p1y+p1z*p1z + massPionC*massPionC);
    el2E = sqrt( p2x*p2x+p2y*p2y+p2z*p2z + massPionC*massPionC);


    el1lv = new TLorentzVector(p1x, p1y, p1z, el1E);
    el2lv = new TLorentzVector(p2x, p2y, p2z, el2E);
    eelv  = new TLorentzVector( (*el1lv)+(*el2lv) );

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

    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).E() );
    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).Pt() );


    //
    //variable of interest
    eelvmass = (*eelv).M();

    ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);

    ((TH1F*)gDirectory->FindObject(h_angle_e1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()) );

    //???
    if ( closap_double_(elp1, elp2, elv1, elv2, &cda, vertex) )
        vertex_OK = TRUE;
    else
        vertex_OK = FALSE;


    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"electron_positron");
    cutcounter++;



    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~
    //"e1e2_hod_Tmatch"
    if ( el1cltime-el2cltime < -3 || el1cltime-el2cltime > 3 ) {cleanup();return 0;}
    if (  el1trkhodtime-el2trkhodtime < -3 || el1trkhodtime-el2trkhodtime > 3 ) {cleanup();return 0;}
    if(IS_DATA)    {
        //after final selection this difference has a maximum at -5.5 and sigma 1.5
        //5ns is according to Manuel's thesis
        if( fabs(el1trktime - el1trkhodtime)>12. )     return 0;
        if( fabs(el2trktime - el2trkhodtime)>12. )     return 0;
    }
    else    {
        if( el1trkhodtime<136. || el1trkhodtime>146. )     {cleanup();return 0;}
        if( el2trkhodtime<136. || el2trkhodtime>146. )     {cleanup();return 0;}
        //        if( el1trkhodtime<80. || el1trkhodtime>150. )     {cleanup();return 0;}
        //        if( el2trkhodtime<80. || el2trkhodtime>150. )     {cleanup();return 0;}
    }
    //
    dirs[cutcounter-1]->cd();
    
    ((TH1I*)gDirectory->FindObject(h_evttypebytrks->GetName()))->Fill(evttypebytrks);
    ((TH1I*)gDirectory->FindObject(h_nelectrack->GetName()))->Fill(nelectrack);
    //((TH1I*)gDirectory->FindObject(h_mutrack->GetName()))->Fill(nmutrack);
    ((TH1I*)gDirectory->FindObject(h_mutrack->GetName()))->Fill(muchtrkindx);
    ((TH1I*)gDirectory->FindObject(h_mutrackq->GetName()))->Fill(muchtrkq);
    ((TH2I*)gDirectory->FindObject(h_neltrk_mutrk->GetName()))->Fill(nelectrack, nmutrack);

    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"e1e2_hod_Tmatch");
    cutcounter++;


   
    /// **********  charged muon
    //cout << "Mu charge track index = " << muchtrkindx << endl;
    //printf("abcog_params.cogX1p=%f,abcog_params.cogX1n=%f\n",abcog_params.cogX1p,abcog_params.cogX1n);
    //printf("abcog_params.cogX4p=%f,abcog_params.cogX4n=%f\n",abcog_params.cogX4p,abcog_params.cogX4n);

    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~
    //"mu_track"
    if (muchtrkindx==-1)    {cleanup();return 0;}
    //if (nmutrack!=1)    {cleanup();return 0;}
    //
    dirs[cutcounter-1]->cd();

       
    ((TH1I*)gDirectory->FindObject(h_evttypebytrks->GetName()))->Fill(evttypebytrks);
    ((TH1I*)gDirectory->FindObject(h_nelectrack->GetName()))->Fill(nelectrack);
    //((TH1I*)gDirectory->FindObject(h_mutrack->GetName()))->Fill(nmutrack);
    ((TH1I*)gDirectory->FindObject(h_mutrack->GetName()))->Fill(muchtrkindx);
    ((TH1I*)gDirectory->FindObject(h_mutrackq->GetName()))->Fill(muchtrkq);
    ((TH2I*)gDirectory->FindObject(h_neltrk_mutrk->GetName()))->Fill(nelectrack, nmutrack);
    ((TH1F*)gDirectory->FindObject(h_dtrkcl->GetName()))     ->Fill( dtrkcl );
    

    ((TH1F*)gDirectory->FindObject(h_e1e2timediff->GetName()))        ->Fill( el1trktime-el2trktime );
    ((TH1F*)gDirectory->FindObject(h_e1e2cltimediff->GetName()))    ->Fill( el1cltime-el2cltime );

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
    
    ((TH1F*)gDirectory->FindObject(h_eop_mu->GetName()))->Fill(sevt->cluster[sevt->track[muchtrkindx].iClus].energy/sevt->track[muchtrkindx].p);
    ((TH2F*)gDirectory->FindObject(h_eop_mu_vs_p->GetName()))->Fill(sevt->cluster[sevt->track[muchtrkindx].iClus].energy/sevt->track[muchtrkindx].p, sevt->track[muchtrkindx].p);

     ((TH1F*)gDirectory->FindObject(h_eop_el->GetName()))->Fill(sevt->cluster[sevt->track[eltrkindx1].iClus].energy/sevt->track[eltrkindx1].p);
     ((TH2F*)gDirectory->FindObject(h_eop_el_vs_p->GetName()))->Fill(sevt->cluster[sevt->track[eltrkindx1].iClus].energy/sevt->track[eltrkindx1].p, sevt->track[eltrkindx1].p);

     ((TH1F*)gDirectory->FindObject(h_eop_el->GetName()))->Fill(sevt->cluster[sevt->track[eltrkindx2].iClus].energy/sevt->track[eltrkindx2].p);
    ((TH2F*)gDirectory->FindObject(h_eop_el_vs_p->GetName()))->Fill(sevt->cluster[sevt->track[eltrkindx2].iClus].energy/sevt->track[eltrkindx2].p, sevt->track[eltrkindx2].p);
    
    muchtrktime         = sevt->track[muchtrkindx].time;
    muchcltime         = sevt->cluster[ sevt->track[muchtrkindx].iClus ].time;
    muchtrkhodtime     = sevt->track[muchtrkindx].hodTime;
    icl_track       = sevt->track[muchtrkindx].iClus;

    ((TH1F*)gDirectory->FindObject(h_muchtrktime->GetName()))        ->Fill(muchtrktime);
    ((TH1F*)gDirectory->FindObject(h_muchtrkhodtime->GetName()))        ->Fill( muchtrkhodtime );
    ((TH1F*)gDirectory->FindObject(h_muchcltime->GetName()))            ->Fill(muchcltime);
    ((TH1F*)gDirectory->FindObject(h_muchtrktimediff->GetName()))    ->Fill( muchtrktime-muchtrkhodtime );
    ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( muchtrktime-muchcltime );
    ((TH1F*)gDirectory->FindObject(h_e1e2hodtimediff->GetName()))        ->Fill( el1trkhodtime-el2trkhodtime );
    ((TH1F*)gDirectory->FindObject(h_e1muhodtimediff->GetName()))        ->Fill( el1trkhodtime-muchtrkhodtime );
    ((TH1F*)gDirectory->FindObject(h_e2muhodtimediff->GetName()))        ->Fill( el2trkhodtime-muchtrkhodtime );



    
    muchv_uncorr[0] = sevt->track[muchtrkindx].bdxdz;
    muchv_uncorr[1] = sevt->track[muchtrkindx].bdydz;
    muchv_uncorr[2] = 1.;
    muchp_uncorr[0] = sevt->track[muchtrkindx].bx;
    muchp_uncorr[1] = sevt->track[muchtrkindx].by;
    muchp_uncorr[2] = DCHbz;
   
    normmu = 1./sqrt(muchv_uncorr[0]*muchv_uncorr[0] + muchv_uncorr[1]*muchv_uncorr[1] + muchv_uncorr[2]*muchv_uncorr[2] );
    pmux = normmu * muchv_uncorr[0] * sevt->track[muchtrkindx].p ;
    pmuy = normmu * muchv_uncorr[1] *  sevt->track[muchtrkindx].p ;
    pmuz = normmu *                   sevt->track[muchtrkindx].p ;
    muchE = sqrt( pmux*pmux+pmuy*pmuy+pmuz*pmuz + massPionC*massPionC);

    pK1x = 0;
    pK1y = 0;                                   
    pK1z = 60.;                                     
    pKE =  sqrt( 60.*60. + massKaonC*massKaonC);                             


    pKaon_simple = new TLorentzVector(pK1x , pK1y , pK1z , pKE);
    muchlv = new TLorentzVector(pmux, pmuy, pmuz, muchE);
    muvee = new TLorentzVector( (*el1lv) + (*el2lv) + (*muchlv) );

    muvee_neutrino = new TLorentzVector((*pKaon_simple) - (*muchlv) - (*el1lv) - (*el2lv));
    // find vertex of mu and ee
    muchv[0] = trkCorrSlopesV[muchtrkindx]->x();
    muchv[1] = trkCorrSlopesV[muchtrkindx]->y();
    muchv[2] = 1.;
    muchp[0] = trkCorrMidPointsV[muchtrkindx]->x();
    muchp[1] = trkCorrMidPointsV[muchtrkindx]->y();
    muchp[2] = trkCorrMidPointsV[muchtrkindx]->z();

    if ((*muchlv).P() < 10.)    {cleanup();return 0;}
    

    if ( closap_double_(elp1, muchp, elv1, muchv, &cda_e1much, vertex_e1much) )
        muchee_vertex_OK = TRUE;
    else
        muchee_vertex_OK = FALSE;

    
    ((TH1F*)gDirectory->FindObject(h_nu_P->GetName()))->Fill(  (*muvee_neutrino).P() );
    ((TH1F*)gDirectory->FindObject(h_nu_E->GetName()))->Fill(  (*muvee_neutrino).E() );
    ((TH1F*)gDirectory->FindObject(h_nu_Pt->GetName()))->Fill( (*muvee_neutrino).Pt() );
    ((TH1F*)gDirectory->FindObject(h_nu_m->GetName()))->Fill(  (*muvee_neutrino).M() );
    ((TH1F*)gDirectory->FindObject(h_missing_mass->GetName()))->Fill( (*muvee_neutrino).M2() );
    //
    ((TH1F*)gDirectory->FindObject(h_KE->GetName()))->Fill((*pKaon_simple).E());
    ((TH1F*)gDirectory->FindObject(h_KP->GetName()))->Fill((*pKaon_simple).P());
    ((TH1F*)gDirectory->FindObject(h_muchlvE->GetName()))->Fill( (*muchlv).E() );
    ((TH1F*)gDirectory->FindObject(h_muchlvP->GetName()))->Fill( (*muchlv).P() );
    ((TH1F*)gDirectory->FindObject(h_muchlvPt->GetName()))->Fill((*muchlv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_angle_e1e2->GetName()))->Fill(     (*el1lv).Vect().Angle(el2lv->Vect()) );
    
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

    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).E() );
    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).Pt() );
    ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);

    
    ((TH1F*)gDirectory->FindObject(h_mueeE->GetName()))->Fill( (*muvee).E());
    ((TH1F*)gDirectory->FindObject(h_mueeP->GetName()))->Fill( (*muvee).P());
    ((TH1F*)gDirectory->FindObject(h_mueePt->GetName()))->Fill((*muvee).Pt());
    ((TH1F*)gDirectory->FindObject(h_mueeM->GetName()))->Fill ((*muvee).M());
    ((TH1F*)gDirectory->FindObject(h_3trk_invm->GetName()))->Fill ((*muvee).M2());
    
    ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill( elp1_uncorr[0] );
    ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill( elp1_uncorr[1] );
    ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill( elv1_uncorr[0] );
    ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill( elv1_uncorr[1] );
    ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill( elp1_uncorr[0], elp1_uncorr[1] );
    ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill( elv1_uncorr[0], elv1_uncorr[1] );
    

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

    ((TH1F*)gDirectory->FindObject(h_cda_e1much->GetName()))->Fill(cda_e1much);
    ((TH2F*)gDirectory->FindObject(h_vx_vs_vy_e1much->GetName()))->Fill(vertex_e1much[0],vertex_e1much[1]);

    ((TH1F*)gDirectory->FindObject(h_e1e2hodtimediff->GetName()))        ->Fill( el1trkhodtime-el2trkhodtime );
    ((TH1F*)gDirectory->FindObject(h_e1muhodtimediff->GetName()))        ->Fill( el1trkhodtime-muchtrkhodtime );
    ((TH1F*)gDirectory->FindObject(h_e2muhodtimediff->GetName()))        ->Fill( el2trkhodtime-muchtrkhodtime );




    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"mu_track");
    cutcounter++;
    
    //SO FAR !!!!


    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~
    //"mu_hod_Tmatch"
    if(IS_DATA)    {
        //after final selection this difference has a maximum at -5.5 and sigma 1.5
        if( fabs(muchtrktime-muchtrkhodtime) > 12. )    return 0;  //5 is from Manuel's thesis
    }
    else    {
        if( muchtrkhodtime<135. || muchtrkhodtime>148. )    return 0;
    }
    //
    dirs[cutcounter-1]->cd();


    ((TH1F*)gDirectory->FindObject(h_bx_much->GetName()))            ->Fill(muchp_uncorr[0]);
    ((TH1F*)gDirectory->FindObject(h_by_much->GetName()))            ->Fill(muchp_uncorr[1]);
    ((TH1F*)gDirectory->FindObject(h_bdxdz_much->GetName()))            ->Fill(muchv_uncorr[0]);
    ((TH1F*)gDirectory->FindObject(h_bdydz_much->GetName()))            ->Fill(muchv_uncorr[1]) ;
    ((TH1F*)gDirectory->FindObject(h_bx_vs_by_much->GetName()))        ->Fill(muchp_uncorr[0], muchp_uncorr[1]);
    ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_much->GetName()))->Fill(muchv_uncorr[0], muchv_uncorr[1]);


    file1->cd();
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"mu_hod_Tmatch");
    h_cutflow->Fill(cutcounter);
    cutcounter++;

    


    //------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------
    //                            electron and positron matching
    //------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------


    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~
    //"e1e2_cl_Tmatch")
    //    if( fabs((long double)(el1cltime-el2cltime))>15. )        {cleanup();return 0;}
    if( fabs((long double)(el1cltime-el2cltime))>10. )        {cleanup();return 0;}
    //
    dirs[cutcounter-1]->cd();

    //    ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
    ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );
    
    ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);

    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"e1e2_cl_Tmatch");
    cutcounter++;



    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~
    //"e1e2_DCHinnerR")
    if( pow(elp1_uncorr[0],2.) + pow(elp1_uncorr[1],2.) < 196. )    {cleanup();return 0;} //14*14 /it was 12*12/
    if( pow(elp2_uncorr[0],2.) + pow(elp2_uncorr[1],2.) < 196. )    {cleanup();return 0;}
    //
    dirs[cutcounter-1]->cd();

    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
    ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);

   
    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"e1e2_DCHinnerR");
    cutcounter++;



    
    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~
    //"e1e2_DCHouterR"
    //    if( pow(elp1_uncorr[0],2.) + pow(elp1_uncorr[1],2.) > 13225. )    {cleanup();return 0;} //115*115
    //    if( pow(elp2_uncorr[0],2.) + pow(elp2_uncorr[1],2.) > 13225. )    {cleanup();return 0;}
    if( pow(elp1_uncorr[0],2.) + pow(elp1_uncorr[1],2.) > 12100. )    {cleanup();return 0;} //110*110
    if( pow(elp2_uncorr[0],2.) + pow(elp2_uncorr[1],2.) > 12100. )    {cleanup();return 0;}
    //
    dirs[cutcounter-1]->cd();
    
  
    //    ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
    ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);

   
    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"e1e2_DCHouterR");
    cutcounter++;


    //    cout<<"3"<<endl;


    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~
    //"e1e2_vtx_found"

    if(vertex_OK!=1)    return 0;

    //
    dirs[cutcounter-1]->cd();

    
    ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
    ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);

    //
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

    ((TH1F*)gDirectory->FindObject(h_cda_e1much->GetName()))->Fill(cda_e1much);
    ((TH2F*)gDirectory->FindObject(h_vx_vs_vy_e1much->GetName()))->Fill(vertex_e1much[0],vertex_e1much[1]);
     ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
    ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );


    ((TH1F*)gDirectory->FindObject(h_muchlvE->GetName()))->Fill( (*muchlv).E() );
    ((TH1F*)gDirectory->FindObject(h_muchlvP->GetName()))->Fill( (*muchlv).P() );
    ((TH1F*)gDirectory->FindObject(h_muchlvPt->GetName()))->Fill( (*muchlv).Pt() );
    ((TH1F*)gDirectory->FindObject(h_nu_P->GetName()))->Fill(  (*muvee_neutrino).P() );
    ((TH1F*)gDirectory->FindObject(h_nu_E->GetName()))->Fill(  (*muvee_neutrino).E() );
    ((TH1F*)gDirectory->FindObject(h_nu_Pt->GetName()))->Fill( (*muvee_neutrino).Pt() );
    ((TH1F*)gDirectory->FindObject(h_nu_m->GetName()))->Fill(  (*muvee_neutrino).M() );

    ((TH1F*)gDirectory->FindObject(h_missing_mass->GetName()))->Fill( (*muvee_neutrino).M2() );
    
    

    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"e1e2_vtx_found");
    cutcounter++;



    
    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~
   
    //    if( vertex[2]<-2000. || vertex[2]>8000. )    {cleanup();return 0;}
    if (eelvmass <= 0.140)  {cleanup();return 0;}
    if ( (*el1lv).P() < 3. || (*el2lv).P() < 3. || (*muchlv).P() < 10.) {cleanup();return 0;}
    
    dirs[cutcounter-1]->cd();
    
    ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
    ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );
    ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);
    
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

    ((TH1F*)gDirectory->FindObject(h_cda_e1much->GetName()))->Fill(cda_e1much);
    ((TH2F*)gDirectory->FindObject(h_vx_vs_vy_e1much->GetName()))->Fill(vertex_e1much[0],vertex_e1much[1]);
     ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
    ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );


    ((TH1F*)gDirectory->FindObject(h_muchlvE->GetName()))->Fill( (*muchlv).E() );
    ((TH1F*)gDirectory->FindObject(h_muchlvP->GetName()))->Fill( (*muchlv).P() );
    ((TH1F*)gDirectory->FindObject(h_muchlvPt->GetName()))->Fill( (*muchlv).Pt() );
    ((TH1F*)gDirectory->FindObject(h_nu_P->GetName()))->Fill(  (*muvee_neutrino).P() );
    ((TH1F*)gDirectory->FindObject(h_nu_E->GetName()))->Fill(  (*muvee_neutrino).E() );
    ((TH1F*)gDirectory->FindObject(h_nu_Pt->GetName()))->Fill( (*muvee_neutrino).Pt() );
    ((TH1F*)gDirectory->FindObject(h_nu_m->GetName()))->Fill(  (*muvee_neutrino).M() );
    
    ((TH1F*)gDirectory->FindObject(h_missing_mass->GetName()))->Fill( (*muvee_neutrino).M2() );
    
    
    
 
    
    
    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"ee_mass_140_and_pcut");
    cutcounter++;
    
     //"e1e2_vtx_good" 
    if( vertex[2]<-1500. || vertex[2]>7000. )    {cleanup();return 0;}
    //
    dirs[cutcounter-1]->cd();

   ((TH1I*)gDirectory->FindObject(h_mutrackq->GetName()))->Fill(muchtrkq);
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
    
    ((TH1F*)gDirectory->FindObject(h_muchtrktime->GetName()))        ->Fill(muchtrktime);
    ((TH1F*)gDirectory->FindObject(h_muchtrkhodtime->GetName()))        ->Fill( muchtrkhodtime );
    ((TH1F*)gDirectory->FindObject(h_muchcltime->GetName()))            ->Fill(muchcltime);
    ((TH1F*)gDirectory->FindObject(h_muchtrktimediff->GetName()))    ->Fill( muchtrktime-muchtrkhodtime );
    ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( muchtrktime-muchcltime );
    ((TH1F*)gDirectory->FindObject(h_e1e2hodtimediff->GetName()))        ->Fill( el1trkhodtime-el2trkhodtime );
    ((TH1F*)gDirectory->FindObject(h_e1muhodtimediff->GetName()))        ->Fill( el1trkhodtime-muchtrkhodtime );
    ((TH1F*)gDirectory->FindObject(h_e2muhodtimediff->GetName()))        ->Fill( el2trkhodtime-muchtrkhodtime );

    ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill( elp1_uncorr[0] );
    ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill( elp1_uncorr[1] );
    ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill( elv1_uncorr[0] );
    ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill( elv1_uncorr[1] );
    ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill( elp1_uncorr[0], elp1_uncorr[1] );
    ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill( elv1_uncorr[0], elv1_uncorr[1] );
    

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

    ((TH1F*)gDirectory->FindObject(h_cda_e1much->GetName()))->Fill(cda_e1much);
    ((TH2F*)gDirectory->FindObject(h_vx_vs_vy_e1much->GetName()))->Fill(vertex_e1much[0],vertex_e1much[1]);

    ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( el1lv->P() );
    ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( el2lv->P() );
        ((TH1F*)gDirectory->FindObject(h_e1m->GetName()))->Fill( el1lv->M() );
    ((TH1F*)gDirectory->FindObject(h_e2m->GetName()))->Fill( el2lv->M() );

    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el1lv->E() );
    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el2lv->E() );
    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el1lv->P() );
    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el2lv->P() );
    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el1lv->Pt() );
    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el2lv->Pt() );

    
  
    ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
    ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);
    



    ((TH1F*)gDirectory->FindObject(h_mueeE->GetName()))->Fill( (*muvee).E());
    ((TH1F*)gDirectory->FindObject(h_mueeP->GetName()))->Fill( (*muvee).P());
    ((TH1F*)gDirectory->FindObject(h_mueePt->GetName()))->Fill((*muvee).Pt());
    ((TH1F*)gDirectory->FindObject(h_mueeM->GetName()))->Fill( (*muvee).M());
    ((TH1F*)gDirectory->FindObject(h_3trk_invm->GetName()))->Fill( (*muvee).M2());
    ((TH1F*)gDirectory->FindObject(h_nu_P->GetName()))->Fill(  (*muvee_neutrino).P() );
    ((TH1F*)gDirectory->FindObject(h_nu_E->GetName()))->Fill(  (*muvee_neutrino).E() );
    ((TH1F*)gDirectory->FindObject(h_nu_Pt->GetName()))->Fill( (*muvee_neutrino).Pt() );
    ((TH1F*)gDirectory->FindObject(h_nu_m->GetName()))->Fill(  (*muvee_neutrino).M() );

    ((TH1F*)gDirectory->FindObject(h_missing_mass->GetName()))->Fill( (*muvee_neutrino).M2() );
    
    
    ((TH1F*)gDirectory->FindObject(h_muchlvE->GetName()))->Fill( (*muchlv).E() );
    ((TH1F*)gDirectory->FindObject(h_muchlvP->GetName()))->Fill( (*muchlv).P() );
    ((TH1F*)gDirectory->FindObject(h_muchlvPt->GetName()))->Fill( (*muchlv).Pt() );


    ((TH1F*)gDirectory->FindObject(h_angle_e1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()) );

    
    ((TH1F*)gDirectory->FindObject(h_KE->GetName()))->Fill((*pKaon_simple).E());
    ((TH1F*)gDirectory->FindObject(h_KP->GetName()))->Fill((*pKaon_simple).P());

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


    
    ((TH1F*)gDirectory->FindObject(h_eldistDCH1plane->GetName()))->Fill( eldistDCH1plane );
    ((TH1F*)gDirectory->FindObject(h_eldistLKrplane->GetName()))->Fill( eldistLKrplane );
    
    ((TH1F*)gDirectory->FindObject(h_muche1distDHC1plane->GetName()))->Fill(  muche1distDHC1plane );
    ((TH1F*)gDirectory->FindObject(h_muche2distDHC1plane->GetName()))->Fill(  muche2distDHC1plane );
    ((TH1F*)gDirectory->FindObject(h_muche1distLKrplane->GetName()))->Fill( muche1distLKrplane );
    ((TH1F*)gDirectory->FindObject(h_muche2distLKrplane->GetName()))->Fill( muche2distLKrplane );



    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"e1e2_vtx_good");
    cutcounter++;




    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~
    //"e1e2_dist_DCH1"
    //    if(eldistDCH1plane<2.)        return 0;
    if(eldistDCH1plane<1.)        return 0;
    //
    dirs[cutcounter-1]->cd();

 
    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
    ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);
    ((TH1F*)gDirectory->FindObject(h_nu_P->GetName()))->Fill(  (*muvee_neutrino).P() );
    ((TH1F*)gDirectory->FindObject(h_nu_E->GetName()))->Fill(  (*muvee_neutrino).E() );
    ((TH1F*)gDirectory->FindObject(h_nu_Pt->GetName()))->Fill( (*muvee_neutrino).Pt() );
    ((TH1F*)gDirectory->FindObject(h_nu_m->GetName()))->Fill(  (*muvee_neutrino).M() );
    
    ((TH1F*)gDirectory->FindObject(h_missing_mass->GetName()))->Fill( (*muvee_neutrino).M2() );
    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"e1e2_dist_DCH1");
    cutcounter++;


    //    cout<<"6"<<endl;


    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~
    //"e1e2_dist_LKr"
    //    if( eldistLKrplane<15. )        return 0;
    if( eldistLKrplane<20. )        {cleanup();return 0;}
    //
    dirs[cutcounter-1]->cd();

   
    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
    ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );
    ((TH1F*)gDirectory->FindObject(h_nu_P->GetName()))->Fill(  (*muvee_neutrino).P() );
    ((TH1F*)gDirectory->FindObject(h_nu_E->GetName()))->Fill(  (*muvee_neutrino).E() );
    ((TH1F*)gDirectory->FindObject(h_nu_Pt->GetName()))->Fill( (*muvee_neutrino).Pt() );
    ((TH1F*)gDirectory->FindObject(h_nu_m->GetName()))->Fill(  (*muvee_neutrino).M() );
    ((TH1F*)gDirectory->FindObject(h_missing_mass->GetName()))->Fill( (*muvee_neutrino).M2() );
    ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);
    ((TH1F*)gDirectory->FindObject(h_mueeE->GetName()))->Fill(muvee->E());
    ((TH1F*)gDirectory->FindObject(h_mueeP->GetName()))->Fill(muvee->P());
    ((TH1F*)gDirectory->FindObject(h_mueePt->GetName()))->Fill(muvee->Pt());
    ((TH1F*)gDirectory->FindObject(h_mueeM->GetName()))->Fill (muvee->M());
    ((TH1F*)gDirectory->FindObject(h_mueeM->GetName()))->Fill (muvee->M2());

   

     file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"e1e2_dist_LKr");
    cutcounter++;

				  

    //------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------
    //                        charged muon and electron-positron matching
    //------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------


    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~
    //"mu_ee_Tmatch"
    if( fabs((long double)(el1trktime-muchtrktime)) > 10. )    {cleanup();return 0;}
    if( fabs((long double)(el1cltime-muchcltime))   > 10. )    {cleanup();return 0;}
    //
    dirs[cutcounter-1]->cd();

    ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );
    
    ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);
    ((TH1F*)gDirectory->FindObject(h_nu_P->GetName()))->Fill(  (*muvee_neutrino).P() );
    ((TH1F*)gDirectory->FindObject(h_nu_E->GetName()))->Fill(  (*muvee_neutrino).E() );
    ((TH1F*)gDirectory->FindObject(h_nu_Pt->GetName()))->Fill( (*muvee_neutrino).Pt() );
    ((TH1F*)gDirectory->FindObject(h_nu_m->GetName()))->Fill(  (*muvee_neutrino).M() );

    
    ((TH1F*)gDirectory->FindObject(h_missing_mass->GetName()))->Fill( (*muvee_neutrino).M2() );
     
    //gpchdistLKrplane    = sqrt( pow((pchxLKrplane-gclx),2.)         + pow((pchyLKrplane-gcly),2.) ) ;

    ((TH1F*)gDirectory->FindObject(h_muche1distDHC1plane->GetName()))->Fill( muche1distDHC1plane );
    ((TH1F*)gDirectory->FindObject(h_muche2distDHC1plane->GetName()))->Fill( muche2distDHC1plane );
    ((TH1F*)gDirectory->FindObject(h_muche1distLKrplane->GetName()))->Fill(  muche1distLKrplane );
    ((TH1F*)gDirectory->FindObject(h_muche2distLKrplane->GetName()))->Fill(  muche2distLKrplane );
    //((TH1F*)gDirectory->FindObject(h_gpchdistLKrplane->GetName()))->Fill( gpchdistLKrplane );
    
    ((TH1F*)gDirectory->FindObject(h_mueeE->GetName()))->Fill(muvee->E());
    ((TH1F*)gDirectory->FindObject(h_mueeP->GetName()))->Fill(muvee->P());
    ((TH1F*)gDirectory->FindObject(h_mueePt->GetName()))->Fill(muvee->Pt());
    ((TH1F*)gDirectory->FindObject(h_mueeM->GetName()))->Fill (muvee->M());
    ((TH1F*)gDirectory->FindObject(h_mueeM->GetName()))->Fill (muvee->M2());
        
    ((TH1F*)gDirectory->FindObject(h_muchlvE->GetName()))->Fill( (*muchlv).E() );
    ((TH1F*)gDirectory->FindObject(h_muchlvP->GetName()))->Fill( (*muchlv).P() );
    ((TH1F*)gDirectory->FindObject(h_muchlvPt->GetName()))->Fill( (*muchlv).Pt() );


    

    file1->cd();
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"mu_ee_Tmatch")
    h_cutflow->Fill(cutcounter);
    cutcounter++;



    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~
    //"mu_vtx_good"
    if( !muchee_vertex_OK)    {cleanup();return 0;}
    //    if( fabs((long double)(vertex[2] - vertex_e1pch[2]))  > 1500. )    {cleanup();return 0;} //old value 300
    if( vertex_e1much[2]<-1500.  || vertex_e1much[2]>7000. )    {cleanup();return 0;}
    //
    dirs[cutcounter-1]->cd();



    ((TH1I*)gDirectory->FindObject(h_mutrackq->GetName()))->Fill(muchtrkq);
    ((TH1F*)gDirectory->FindObject(h_el1trktime->GetName()))        ->Fill( el1trktime );
    ((TH1F*)gDirectory->FindObject(h_el1trkhodtime->GetName()))        ->Fill( el1trkhodtime );
    ((TH1F*)gDirectory->FindObject(h_el1trktimediff->GetName()))    ->Fill( el1trktime - el1trkhodtime );
    ((TH1F*)gDirectory->FindObject(h_ecltime1->GetName()))            ->Fill( el1cltime );
    ((TH1F*)gDirectory->FindObject(h_e1trkcltimediff->GetName()))    ->Fill( el1trktime- el1cltime );
    ((TH1F*)gDirectory->FindObject(h_dtrkcl->GetName()))     ->Fill( dtrkcl );
    
    ((TH1F*)gDirectory->FindObject(h_el2trktime->GetName()))        ->Fill( el2trktime );
    ((TH1F*)gDirectory->FindObject(h_el2trkhodtime->GetName()))        ->Fill( el2trkhodtime );
    ((TH1F*)gDirectory->FindObject(h_el2trktimediff->GetName()))    ->Fill( el2trktime - el2trkhodtime );
    ((TH1F*)gDirectory->FindObject(h_ecltime2->GetName()))            ->Fill( el2cltime );
    ((TH1F*)gDirectory->FindObject(h_e2trkcltimediff->GetName()))    ->Fill( el2trktime-el2cltime );
  
    ((TH1F*)gDirectory->FindObject(h_e1e2timediff->GetName()))        ->Fill( el1trktime-el2trktime );
    ((TH1F*)gDirectory->FindObject(h_e1e2cltimediff->GetName()))    ->Fill( el1cltime-el2cltime );
    
    ((TH1F*)gDirectory->FindObject(h_muchtrktime->GetName()))        ->Fill(muchtrktime);
    ((TH1F*)gDirectory->FindObject(h_muchtrkhodtime->GetName()))        ->Fill( muchtrkhodtime );
    ((TH1F*)gDirectory->FindObject(h_muchcltime->GetName()))            ->Fill(muchcltime);
    ((TH1F*)gDirectory->FindObject(h_muchtrktimediff->GetName()))    ->Fill( muchtrktime-muchtrkhodtime );
    ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( muchtrktime-muchcltime );
    ((TH1F*)gDirectory->FindObject(h_e1e2hodtimediff->GetName()))        ->Fill( el1trkhodtime-el2trkhodtime );
    ((TH1F*)gDirectory->FindObject(h_e1muhodtimediff->GetName()))        ->Fill( el1trkhodtime-muchtrkhodtime );
    ((TH1F*)gDirectory->FindObject(h_e2muhodtimediff->GetName()))        ->Fill( el2trkhodtime-muchtrkhodtime );

    ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill( elp1_uncorr[0] );
    ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill( elp1_uncorr[1] );
    ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill( elv1_uncorr[0] );
    ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill( elv1_uncorr[1] );
    ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill( elp1_uncorr[0], elp1_uncorr[1] );
    ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill( elv1_uncorr[0], elv1_uncorr[1] );
    

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

    ((TH1F*)gDirectory->FindObject(h_cda_e1much->GetName()))->Fill(cda_e1much);
    ((TH2F*)gDirectory->FindObject(h_vx_vs_vy_e1much->GetName()))->Fill(vertex_e1much[0],vertex_e1much[1]);

    ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( (*el1lv).P() );
    ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( (*el2lv).P() );
    ((TH1F*)gDirectory->FindObject(h_e1m->GetName()))->Fill( el1lv->M() );
    ((TH1F*)gDirectory->FindObject(h_e2m->GetName()))->Fill( el2lv->M() );

    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill(  (*el1lv).E() );
    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill(  (*el2lv).E() );
    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill(  (*el1lv).P() );
    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill(  (*el2lv).P() );
    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( (*el1lv).Pt() );
    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( (*el2lv).Pt() );

    
  
    ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
    ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);
    



    ((TH1F*)gDirectory->FindObject(h_mueeE->GetName()))->Fill( (*muvee).E());
    ((TH1F*)gDirectory->FindObject(h_mueeP->GetName()))->Fill( (*muvee).P());
    ((TH1F*)gDirectory->FindObject(h_mueePt->GetName()))->Fill((*muvee).Pt());
    ((TH1F*)gDirectory->FindObject(h_mueeM->GetName()))->Fill( (*muvee).M());
    ((TH1F*)gDirectory->FindObject(h_3trk_invm->GetName()))->Fill( (*muvee).M2());

    ((TH1F*)gDirectory->FindObject(h_nu_P->GetName()))->Fill(  (*muvee_neutrino).P() );
    ((TH1F*)gDirectory->FindObject(h_nu_E->GetName()))->Fill(  (*muvee_neutrino).E() );
    ((TH1F*)gDirectory->FindObject(h_nu_Pt->GetName()))->Fill( (*muvee_neutrino).Pt() );
    ((TH1F*)gDirectory->FindObject(h_nu_m->GetName()))->Fill(  (*muvee_neutrino).M() );
    
    ((TH1F*)gDirectory->FindObject(h_missing_mass->GetName()))->Fill( (*muvee_neutrino).M2() );
    
    
    ((TH1F*)gDirectory->FindObject(h_muchlvE->GetName()))->Fill( (*muchlv).E() );
    ((TH1F*)gDirectory->FindObject(h_muchlvP->GetName()))->Fill( (*muchlv).P() );
    ((TH1F*)gDirectory->FindObject(h_muchlvPt->GetName()))->Fill( (*muchlv).Pt() );



    file1->cd();
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"mu_vtx_good")
    h_cutflow->Fill(cutcounter);
    cutcounter++;



    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~
    //"mu_ee_dist_LKr  & hodtime cut"
    if(el1trkhodtime-muchtrkhodtime < -2 || el1trkhodtime-muchtrkhodtime > 2 )    {cleanup();return 0;}
    if(el2trkhodtime-muchtrkhodtime < -2 || el2trkhodtime-muchtrkhodtime > 2 )    {cleanup();return 0;}
    if(muche1distLKrplane<30.)    {cleanup();return 0;}
    if(muche2distLKrplane<30.)    {cleanup();return 0;}
    //
    dirs[cutcounter-1]->cd();


    ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill(elp1_uncorr[0]);
    ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill(elp1_uncorr[1]);
    ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill(elv1_uncorr[0]);
    ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill(elv1_uncorr[1]);
    ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill(elp1_uncorr[0], elp1_uncorr[1] );
    ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill(elv1_uncorr[0], elv1_uncorr[1] );
   
    ((TH1F*)gDirectory->FindObject(h_bx_much->GetName()))            ->Fill(muchp_uncorr[0]);
    ((TH1F*)gDirectory->FindObject(h_by_much->GetName()))            ->Fill(muchp_uncorr[1]);
    ((TH1F*)gDirectory->FindObject(h_bdxdz_much->GetName()))            ->Fill(muchv_uncorr[0]);
    ((TH1F*)gDirectory->FindObject(h_bdydz_much->GetName()))            ->Fill(muchv_uncorr[1]) ;
    ((TH1F*)gDirectory->FindObject(h_bx_vs_by_much->GetName()))        ->Fill(muchp_uncorr[0], muchp_uncorr[1]);
    ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_much->GetName()))->Fill(  muchv_uncorr[0], muchv_uncorr[1]);
    
    ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( (*el1lv).P() );
    ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( (*el2lv).P() );
    ((TH1F*)gDirectory->FindObject(h_e1m->GetName()))->Fill( el1lv->M() );
    ((TH1F*)gDirectory->FindObject(h_e2m->GetName()))->Fill( el2lv->M() );

    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill(  (*el1lv).E() );
    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill(  (*el2lv).E() );
    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill(  (*el1lv).P() );
    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill(  (*el2lv).P() );
    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( (*el1lv).Pt() );
    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( (*el2lv).Pt() );

      
    ((TH1F*)gDirectory->FindObject(h_muchlvE->GetName()))->Fill( (*muchlv).E() );
    ((TH1F*)gDirectory->FindObject(h_muchlvP->GetName()))->Fill( (*muchlv).P() );
    ((TH1F*)gDirectory->FindObject(h_muchlvPt->GetName()))->Fill( (*muchlv).Pt() );

    ((TH1I*)gDirectory->FindObject(h_evttypebytrks->GetName()))->Fill(evttypebytrks);
    ((TH1I*)gDirectory->FindObject(h_nelectrack->GetName()))->Fill(nelectrack);
    ((TH1I*)gDirectory->FindObject(h_mutrack->GetName()))->Fill(nmutrack);
    ((TH1I*)gDirectory->FindObject(h_mutrackq->GetName()))->Fill(muchtrkq);
    ((TH2I*)gDirectory->FindObject(h_neltrk_mutrk->GetName()))->Fill(nelectrack, nmutrack);

    ((TH1F*)gDirectory->FindObject(h_eop_mu->GetName()))->Fill(sevt->cluster[sevt->track[muchtrkindx].iClus].energy/sevt->track[muchtrkindx].p);
    ((TH2F*)gDirectory->FindObject(h_eop_mu_vs_p->GetName()))->Fill(sevt->cluster[sevt->track[muchtrkindx].iClus].energy/sevt->track[muchtrkindx].p, sevt->track[muchtrkindx].p);

    ((TH1F*)gDirectory->FindObject(h_eop_el->GetName()))->Fill(sevt->cluster[sevt->track[eltrkindx1].iClus].energy/sevt->track[eltrkindx1].p);
    ((TH2F*)gDirectory->FindObject(h_eop_el_vs_p->GetName()))->Fill(sevt->cluster[sevt->track[eltrkindx1].iClus].energy/sevt->track[eltrkindx1].p, sevt->track[eltrkindx1].p);

    ((TH1F*)gDirectory->FindObject(h_eop_el->GetName()))->Fill(sevt->cluster[sevt->track[eltrkindx2].iClus].energy/sevt->track[eltrkindx2].p);
    ((TH2F*)gDirectory->FindObject(h_eop_el_vs_p->GetName()))->Fill(sevt->cluster[sevt->track[eltrkindx2].iClus].energy/sevt->track[eltrkindx2].p, sevt->track[eltrkindx2].p);
	

  
    ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
    ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);
    ((TH1F*)gDirectory->FindObject(h_nu_P->GetName()))->Fill(  (*muvee_neutrino).P() );
    ((TH1F*)gDirectory->FindObject(h_nu_E->GetName()))->Fill(  (*muvee_neutrino).E() );
    ((TH1F*)gDirectory->FindObject(h_nu_Pt->GetName()))->Fill( (*muvee_neutrino).Pt() );
    ((TH1F*)gDirectory->FindObject(h_nu_m->GetName()))->Fill(  (*muvee_neutrino).M() );

    ((TH1F*)gDirectory->FindObject(h_missing_mass->GetName()))->Fill( (*muvee_neutrino).M2() );
    ((TH1F*)gDirectory->FindObject(h_mueeE->GetName()))->Fill((*muvee).E());
    ((TH1F*)gDirectory->FindObject(h_mueeP->GetName()))->Fill((*muvee).P());
    ((TH1F*)gDirectory->FindObject(h_mueePt->GetName()))->Fill((*muvee).Pt());
    ((TH1F*)gDirectory->FindObject(h_mueeM->GetName()))->Fill ((*muvee).M());
    ((TH1F*)gDirectory->FindObject(h_3trk_invm->GetName()))->Fill( (*muvee).M2());

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

    ((TH1F*)gDirectory->FindObject(h_muchtrktime->GetName()))        ->Fill(muchtrktime);
    ((TH1F*)gDirectory->FindObject(h_muchtrkhodtime->GetName()))        ->Fill( muchtrkhodtime );
    ((TH1F*)gDirectory->FindObject(h_muchcltime->GetName()))            ->Fill(muchcltime);
    ((TH1F*)gDirectory->FindObject(h_muchtrktimediff->GetName()))    ->Fill( muchtrktime-muchtrkhodtime );
    ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( muchtrktime-muchcltime );
    ((TH1F*)gDirectory->FindObject(h_e1e2hodtimediff->GetName()))        ->Fill( el1trkhodtime-el2trkhodtime );
    ((TH1F*)gDirectory->FindObject(h_e1muhodtimediff->GetName()))        ->Fill( el1trkhodtime-muchtrkhodtime );
    ((TH1F*)gDirectory->FindObject(h_e2muhodtimediff->GetName()))        ->Fill( el2trkhodtime-muchtrkhodtime );


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

    ((TH1F*)gDirectory->FindObject(h_cda_e1much->GetName()))->Fill(cda_e1much);
    ((TH2F*)gDirectory->FindObject(h_vx_vs_vy_e1much->GetName()))->Fill(vertex_e1much[0],vertex_e1much[1]);

    
    ((TH1F*)gDirectory->FindObject(h_muche1distDHC1plane->GetName()))->Fill( muche1distDHC1plane );
    ((TH1F*)gDirectory->FindObject(h_muche2distDHC1plane->GetName()))->Fill( muche2distDHC1plane );
    
    
    ((TH1F*)gDirectory->FindObject(h_muche1distLKrplane->GetName()))->Fill(  muche1distLKrplane );
    ((TH1F*)gDirectory->FindObject(h_muche2distLKrplane->GetName()))->Fill(  muche2distLKrplane );
    
    ((TH1F*)gDirectory->FindObject(h_eldistDCH1plane->GetName()))->Fill( eldistDCH1plane );
    ((TH1F*)gDirectory->FindObject(h_eldistLKrplane->GetName()))->Fill( eldistLKrplane );
 

    //    ((TH1F*)gDirectory->FindObject(h_gE->GetName()))->Fill( (*glv).E() );
    //    ((TH1F*)gDirectory->FindObject(h_gP->GetName()))->Fill( (*glv).P() );
    //    ((TH1F*)gDirectory->FindObject(h_gPt->GetName()))->Fill( (*glv).Pt() );
    //    ((TH1F*)gDirectory->FindObject(h_gtrkRadDHC1->GetName()))->Fill( gtrkRadDHC1 );
    //
    //    ((TH1F*)gDirectory->FindObject(h_pchlvE->GetName()))->Fill( (*pchlv).E() );
    //    ((TH1F*)gDirectory->FindObject(h_pchlvP->GetName()))->Fill( (*pchlv).P() );
    //    ((TH1F*)gDirectory->FindObject(h_pchlvPt->GetName()))->Fill( (*pchlv).Pt() );
    //
    //
    //    ((TH1F*)gDirectory->FindObject(h_angle_e1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_angle_ph_e1e2->GetName()))->Fill( (*eelv).Vect().Angle(glv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_angle_pch_pi0->GetName()))->Fill( (*pi0lv).Vect().Angle(pchlv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_a_e1e2_phe1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()), (*eelv).Vect().Angle(glv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_a_e1e2_pchpi0->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()), (*pi0lv).Vect().Angle(pchlv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_a_phe1e2_pchpi0->GetName()))->Fill( (*eelv).Vect().Angle(glv->Vect()), (*pi0lv).Vect().Angle(pchlv->Vect()) );
    //
    //    ((TH1F*)gDirectory->FindObject(h_cda_ee->GetName()))->Fill(cda);
    //    ((TH1F*)gDirectory->FindObject(h_vx->GetName()))->Fill(vertex[0]);
    //    ((TH1F*)gDirectory->FindObject(h_vy->GetName()))->Fill(vertex[1]);
    //    ((TH1F*)gDirectory->FindObject(h_vz->GetName()))->Fill(vertex[2]);
    //    ((TH2F*)gDirectory->FindObject(h_vx_vs_vy->GetName()))->Fill(vertex[0],vertex[1]);
    //
    //    ((TH1F*)gDirectory->FindObject(h_cda_e1pch->GetName()))->Fill(cda_e1pch);
    //    ((TH1F*)gDirectory->FindObject(h_vx_e1pch->GetName()))->Fill(vertex_e1pch[0]);
    //    ((TH1F*)gDirectory->FindObject(h_vy_e1pch->GetName()))->Fill(vertex_e1pch[1]);
    //    ((TH1F*)gDirectory->FindObject(h_vz_e1pch->GetName()))->Fill(vertex_e1pch[2]);
    //    ((TH2F*)gDirectory->FindObject(h_vx_vs_vy_e1pch->GetName()))->Fill(vertex_e1pch[0],vertex_e1pch[1]);
    //
    //    ((TH1F*)gDirectory->FindObject(h_vz_e1pch_diff->GetName()))->Fill( vertex[2] - vertex_e1pch[2] );

    //    ((TH1F*)gDirectory->FindObject(h_pche1distDHC1plane->GetName()))->Fill( pche1distDHC1plane );
    //    ((TH1F*)gDirectory->FindObject(h_pche2distDHC1plane->GetName()))->Fill( pche2distDHC1plane );
    //    ((TH1F*)gDirectory->FindObject(h_pche1distLKrplane->GetName()))->Fill( pche1distLKrplane );
    //    ((TH1F*)gDirectory->FindObject(h_pche2distLKrplane->GetName()))->Fill( pche2distLKrplane );
    //    ((TH1F*)gDirectory->FindObject(h_gpchdistLKrplane->GetName()))->Fill( gpchdistLKrplane );

    //
    //((TH1F*)gDirectory->FindObject(h_pi0E->GetName()))->Fill(  pi0lv->E() );
    //((TH1F*)gDirectory->FindObject(h_pi0P->GetName()))->Fill(  pi0lv->P() );
    //((TH1F*)gDirectory->FindObject(h_pi0Pt->GetName()))->Fill( pi0lv->Pt() );
    //((TH1F*)gDirectory->FindObject(h_mofpi0->GetName()))->Fill( (*pi0lv).M() );
    //
    //((TH1F*)gDirectory->FindObject(h_pi0Pcm->GetName()))->Fill( pi0lvcm->P() );
    //((TH1F*)gDirectory->FindObject(h_mofpi0diff->GetName()))->Fill( (*pi0lvcm).M()-(*pi0lv).M() );
    //
    //
    //((TH1F*)gDirectory->FindObject(h_gmueeE->GetName()))->Fill(  ((*glv)+(*muvee)).E() );
    //((TH1F*)gDirectory->FindObject(h_gmueeP->GetName()))->Fill(  ((*glv)+(*muvee)).P() );
    //((TH1F*)gDirectory->FindObject(h_gmueePt->GetName()))->Fill( ((*glv)+(*muvee)).Pt());
    //((TH1F*)gDirectory->FindObject(h_gmueeM->GetName()))->Fill(  ((*glv)+(*muvee)).M() );

    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"mu_ee_dist_LKr");
    cutcounter++;

    //"muee energy > 70GeV test"
    if( (*muvee).E() < 70 )        {cleanup();return 0;}
    //
    dirs[cutcounter-1]->cd();
    
    
    ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill(elp1_uncorr[0]);
    ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill(elp1_uncorr[1]);
    ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill(elv1_uncorr[0]);
    ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill(elv1_uncorr[1]);
    ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill(elp1_uncorr[0], elp1_uncorr[1] );
    ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill(elv1_uncorr[0], elv1_uncorr[1] );
   
    ((TH1F*)gDirectory->FindObject(h_bx_much->GetName()))            ->Fill(muchp_uncorr[0]);
    ((TH1F*)gDirectory->FindObject(h_by_much->GetName()))            ->Fill(muchp_uncorr[1]);
    ((TH1F*)gDirectory->FindObject(h_bdxdz_much->GetName()))            ->Fill(muchv_uncorr[0]);
    ((TH1F*)gDirectory->FindObject(h_bdydz_much->GetName()))            ->Fill(muchv_uncorr[1]) ;
    ((TH1F*)gDirectory->FindObject(h_bx_vs_by_much->GetName()))        ->Fill(muchp_uncorr[0], muchp_uncorr[1]);
    ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_much->GetName()))->Fill(  muchv_uncorr[0], muchv_uncorr[1]);
    
    ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( (*el1lv).P() );
    ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( (*el2lv).P() );
    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill(  (*el1lv).E() );
    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill(  (*el2lv).E() );
    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill(  (*el1lv).P() );
    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill(  (*el2lv).P() );
    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( (*el1lv).Pt() );
    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( (*el2lv).Pt() );

      
    ((TH1F*)gDirectory->FindObject(h_muchlvE->GetName()))->Fill( (*muchlv).E() );
    ((TH1F*)gDirectory->FindObject(h_muchlvP->GetName()))->Fill( (*muchlv).P() );
    ((TH1F*)gDirectory->FindObject(h_muchlvPt->GetName()))->Fill( (*muchlv).Pt() );

    ((TH1I*)gDirectory->FindObject(h_evttypebytrks->GetName()))->Fill(evttypebytrks);
    ((TH1I*)gDirectory->FindObject(h_nelectrack->GetName()))->Fill(nelectrack);
    ((TH1I*)gDirectory->FindObject(h_mutrack->GetName()))->Fill(nmutrack);
    ((TH1I*)gDirectory->FindObject(h_mutrackq->GetName()))->Fill(muchtrkq);
    ((TH2I*)gDirectory->FindObject(h_neltrk_mutrk->GetName()))->Fill(nelectrack, nmutrack);

    
    ((TH1F*)gDirectory->FindObject(h_eop_mu->GetName()))->Fill(sevt->cluster[sevt->track[muchtrkindx].iClus].energy/sevt->track[muchtrkindx].p);
    ((TH2F*)gDirectory->FindObject(h_eop_mu_vs_p->GetName()))->Fill(sevt->cluster[sevt->track[muchtrkindx].iClus].energy/sevt->track[muchtrkindx].p, sevt->track[muchtrkindx].p);

     ((TH1F*)gDirectory->FindObject(h_eop_el->GetName()))->Fill(sevt->cluster[sevt->track[eltrkindx1].iClus].energy/sevt->track[eltrkindx1].p);
    ((TH2F*)gDirectory->FindObject(h_eop_el_vs_p->GetName()))->Fill(sevt->cluster[sevt->track[eltrkindx1].iClus].energy/sevt->track[eltrkindx1].p, sevt->track[eltrkindx1].p);

     ((TH1F*)gDirectory->FindObject(h_eop_el->GetName()))->Fill(sevt->cluster[sevt->track[eltrkindx2].iClus].energy/sevt->track[eltrkindx2].p);
    ((TH2F*)gDirectory->FindObject(h_eop_el_vs_p->GetName()))->Fill(sevt->cluster[sevt->track[eltrkindx2].iClus].energy/sevt->track[eltrkindx2].p, sevt->track[eltrkindx2].p);
    
  
    ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
    ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);
    ((TH1F*)gDirectory->FindObject(h_nu_P->GetName()))->Fill(  (*muvee_neutrino).P() );
    ((TH1F*)gDirectory->FindObject(h_nu_E->GetName()))->Fill(  (*muvee_neutrino).E() );
    ((TH1F*)gDirectory->FindObject(h_nu_Pt->GetName()))->Fill( (*muvee_neutrino).Pt() );
    ((TH1F*)gDirectory->FindObject(h_nu_m->GetName()))->Fill(  (*muvee_neutrino).M() );

    ((TH1F*)gDirectory->FindObject(h_missing_mass->GetName()))->Fill( (*muvee_neutrino).M2() );
    ((TH1F*)gDirectory->FindObject(h_mueeE->GetName()))->Fill((*muvee).E());
    ((TH1F*)gDirectory->FindObject(h_mueeP->GetName()))->Fill((*muvee).P());
    ((TH1F*)gDirectory->FindObject(h_mueePt->GetName()))->Fill((*muvee).Pt());
    ((TH1F*)gDirectory->FindObject(h_mueeM->GetName()))->Fill ((*muvee).M());
    ((TH1F*)gDirectory->FindObject(h_3trk_invm->GetName()))->Fill( (*muvee).M2());

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

    ((TH1F*)gDirectory->FindObject(h_muchtrktime->GetName()))        ->Fill(muchtrktime);
    ((TH1F*)gDirectory->FindObject(h_muchtrkhodtime->GetName()))        ->Fill( muchtrkhodtime );
    ((TH1F*)gDirectory->FindObject(h_muchcltime->GetName()))            ->Fill(muchcltime);
    ((TH1F*)gDirectory->FindObject(h_muchtrktimediff->GetName()))    ->Fill( muchtrktime-muchtrkhodtime );
    ((TH1F*)gDirectory->FindObject(h_muchtrkcltimediff->GetName()))    ->Fill( muchtrktime-muchcltime );
    ((TH1F*)gDirectory->FindObject(h_e1e2hodtimediff->GetName()))        ->Fill( el1trkhodtime-el2trkhodtime );
    ((TH1F*)gDirectory->FindObject(h_e1muhodtimediff->GetName()))        ->Fill( el1trkhodtime-muchtrkhodtime );
    ((TH1F*)gDirectory->FindObject(h_e2muhodtimediff->GetName()))        ->Fill( el2trkhodtime-muchtrkhodtime );


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
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str()); //"muee energy > 70GeV test"
    cutcounter++;

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
    if(pKaon_simple)      delete pKaon_simple; pKaon_simple=0;
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