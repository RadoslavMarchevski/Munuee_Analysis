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
TLorentzVector *pchlv=NULL;
TLorentzVector *pieelv=NULL;
TLorentzVector *glv =NULL;
TLorentzVector *pi0lv=NULL;
TLorentzVector *el1lvcm=NULL, *el2lvcm=NULL,  *glvcm=NULL, *pi0lvcm=NULL;
TLorentzVector *Klv=NULL, *pi0lvKcm=NULL, *pchlvKcm=NULL, *KlvKcm=NULL;


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
    static double massMuon = 0.105658369;
    static double massPionC = 0.13957018;
    static double massKaonC = 0.493677;
    //##########################################################
    // End of Cut definitions
    //##########################################################

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
    int    nphcluster=0;
    int    evttypebytrks=-1;
    double el1trkhodtime=-1e+3;
    double el2trkhodtime=-1e+3;
    double el1trktime=-1e+3;
    double el2trktime=-1e+3;

    double norm1, norm2;
    double p1x,p1y, p1z, p2x,p2y, p2z, el1E, el2E;

    int    gclindx=-1;
    double gclenergy=0.;
    double gcltime=0.;
    double gclx=0.;
    double gcly=0.;
    double gclz=0.;
    double gtrkDHC1x, gtrkDHC1y, gtrkRadDHC1 ;


    double pchtrkhodtime, pchE;
    
	double eldistDCH1plane, el1xLKrplane, el1yLKrplane, el2xLKrplane, el2yLKrplane, eldistLKrplane;
    double gel1distLKrplane, gel2distLKrplane;
    double pche1distDHC1plane, pche2distDHC1plane, pche1distLKrplane, pche2distLKrplane;

    double pchxLKrplane, pchyLKrplane, gpchdistLKrplane;

    double Dgvtx_x,Dgvtx_y, Dgvtx_z, gnorm,  gpx, gpy, gpz;

    unsigned int vertex_OK, pchee_vertex_OK;
    double vertex[3], vertex_e2pch[3], vertex_e1pch[3];
    double cda, cda_e2pch, cda_e1pch;

    double elp1_uncorr[3], elp1[3], elv1_uncorr[3], elv1[3], elp2_uncorr[3],elp2[3], elv2_uncorr[3],elv2[3], pchp_uncorr[3],pchp[3], pchv_uncorr[3], pchv[3];

    int    pchtrkindx=-1;
    double pchtrktime=-1e+3;
    int    pchtrkq=0;
    double pchcltime=-1e+3;
    double pcheopopt=0;
    double pcheopposition=0;
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

    dirs[cutcounter-1]->cd();

    h_ntrack->Fill(ntrack);
    h_nclust->Fill(nclust);

    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"no_cuts");
    cutcounter++;



    /// ************* tracks  &  clusters  *******************************

    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~

    if(ntrack!=3 || nclust!=4) {cleanup();return 0;}
    //
    dirs[cutcounter-1]->cd();

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
        if (imu_track>-1)             { trkCorrSlopesV.push_back(NULL); trkCorrMidPointsV.push_back(NULL); continue; }
        if (ddcell_track<2.)        { trkCorrSlopesV.push_back(NULL); trkCorrMidPointsV.push_back(NULL); continue; }
        if (p_track<5.)                { trkCorrSlopesV.push_back(NULL); trkCorrMidPointsV.push_back(NULL); continue; }
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

    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel( cutcounter, cuts[cutcounter-1].c_str() );// "3tracks_4clusters");
    cutcounter++;




    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~

    if(ngoodtrack!=3) {cleanup();return 0;}
    //
    dirs[cutcounter-1]->cd();


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
        if(eop_track>c_eop_e)    {
            //            if(eop_track>c_eop_e  &&  eop_track<c_eop_e_up)    {

            //            cout<< "electron ";
            if(eltrkindx1<0)    {
                eltrkindx1=i;
                eltrkq1 = sevt->track[i].q;
                //                cout<< "first "<<eltrkindx1<<"    icl_track="<<icl_track;

            }
            else if(eltrkindx2<0)    {
                eltrkindx2=i;
                eltrkq1 = sevt->track[i].q;
                //                cout<< "second    "<<eltrkindx2<<"    icl_track="<<icl_track;
            }
            else { ; }

            nelectrack++;
        }
        // pion track
        else if(eop_track<c_eop_e)    {

            //            cout<< "pion ";
            if(pchtrkindx<0)    {
                pchtrkindx=i;
                pchtrkq=sevt->track[i].q;
                //                cout<<pchtrkindx<<"    icl_track="<<icl_track;
            }

            npitrack++;
        }
        //track with too high e/p. skip it.
        else    {
            //            cout<< "?? ";;
        }
        //        cout<< " "<<endl;

    }

    if(nelectrack+npitrack==3)    {
        if(pchtrkq>0.)    {
            if(eltrkq1 != eltrkq2)    {// pi+, e+, e-
                evttypebytrks=0;
            }
            else if(eltrkq1>0)    {    // pi+, e+, e+
                evttypebytrks=1;
            }
            else    {                // pi+, e-, e-
                evttypebytrks=2;
            }
        }
        else    {
            if(eltrkq1 != eltrkq2)    {// pi-, e+, e-
                evttypebytrks=3;
            }
            else if(eltrkq1>0)    {    // pi-, e+, e+
                evttypebytrks=4;
            }
            else    {                // pi-, e-, e-
                evttypebytrks=5;
            }
        }
    }


    ((TH1I*)gDirectory->FindObject(h_evttypebytrks->GetName()))->Fill(evttypebytrks);
    ((TH1I*)gDirectory->FindObject(h_nelectrack->GetName()))->Fill(nelectrack);
    ((TH1I*)gDirectory->FindObject(h_pitrack->GetName()))->Fill(npitrack);
    ((TH2I*)gDirectory->FindObject(h_neltrk_pitrk->GetName()))->Fill(nelectrack, npitrack);
    ((TH1F*)gDirectory->FindObject(h_dtrkcl->GetName()))     ->Fill( dtrkcl );

    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel( cutcounter, cuts[cutcounter-1].c_str() );// "3goodtracks");
    cutcounter++;





    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~

    if(dtrkcl<20.)        {cleanup();return 0;}
    //
    dirs[cutcounter-1]->cd();

    ((TH1I*)gDirectory->FindObject(h_evttypebytrks->GetName()))->Fill(evttypebytrks);
    ((TH1I*)gDirectory->FindObject(h_nelectrack->GetName()))->Fill(nelectrack);
    ((TH1I*)gDirectory->FindObject(h_pitrack->GetName()))->Fill(npitrack);
    ((TH2I*)gDirectory->FindObject(h_neltrk_pitrk->GetName()))->Fill(nelectrack, npitrack);

    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel( cutcounter, cuts[cutcounter-1].c_str() );// "cluster_isolation");
    cutcounter++;





    /// ************* electron / positron

    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~

    if (nelectrack!=2)        {cleanup();return 0;}
    if (eltrkq1==eltrkq2)    {cleanup();return 0;}
    //
    dirs[cutcounter-1]->cd();

    el1trktime      = sevt->track[eltrkindx1].time ;
    el1trkhodtime    = sevt->track[eltrkindx1].hodTime ;
    icl_track       = sevt->track[eltrkindx1].iClus;
    el1cltime        = sevt->cluster[icl_track].time;
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

    //    cout<< "    second "<<eltrkindx2<<"    icl_track="<<icl_track<<"    Dt="<<el2trktime - el2trkhodtime<<endl;

    //e1e2
    ((TH1F*)gDirectory->FindObject(h_e1e2timediff->GetName()))        ->Fill( el1trktime-el2trktime );
    ((TH1F*)gDirectory->FindObject(h_e1e2cltimediff->GetName()))    ->Fill( el1cltime-el2cltime );

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
    el1E = sqrt( p1x*p1x+p1y*p1y+p1z*p1z + massElec*massElec);
    el2E = sqrt( p2x*p2x+p2y*p2y+p2z*p2z + massElec*massElec);

    el1lv = new TLorentzVector(p1x, p1y, p1z, el1E);
    el2lv = new TLorentzVector(p2x, p2y, p2z, el2E);
    eelv  = new TLorentzVector( (*el1lv)+(*el2lv) );

    ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( el1lv->P() );
    ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( el2lv->P() );
    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el1lv->E() );
    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el2lv->E() );
    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el1lv->P() );
    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el2lv->P() );
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


    if ( closap_double_(elp1, elp2, elv1, elv2, &cda, vertex) )
        vertex_OK = TRUE;
    else
        vertex_OK = FALSE;


    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"electron_positron");
    cutcounter++;



    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~

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

    //    ((TH1F*)gDirectory->FindObject(h_el1trktime->GetName()))        ->Fill( el1trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trkhodtime->GetName()))        ->Fill( el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trktimediff->GetName()))    ->Fill( el1trktime - el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime1->GetName()))            ->Fill( el1cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e1trkcltimediff->GetName()))    ->Fill( el1trktime- el1cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_el2trktime->GetName()))        ->Fill( el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trkhodtime->GetName()))        ->Fill( el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trktimediff->GetName()))    ->Fill( el2trktime - el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime2->GetName()))            ->Fill( el2cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e2trkcltimediff->GetName()))    ->Fill( el2trktime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill(elp1_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill(elp1_uncorr[1]);
    //    ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill(elv1_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill(elv1_uncorr[1]);
    //    ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill(elp1_uncorr[0], elp1_uncorr[1] );
    //    ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill(elv1_uncorr[0], elv1_uncorr[1] );
    //
    //    ((TH1F*)gDirectory->FindObject(h_e1e2timediff->GetName()))        ->Fill( el1trktime-el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_e1e2cltimediff->GetName()))    ->Fill( el1cltime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( el1lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( el2lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el1lv->E() );
    //    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el2lv->E() );
    //    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el1lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el2lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el1lv->Pt() );
    //    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el2lv->Pt() );
    //
    //    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).E() );
    //    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
    //    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).Pt() );
    //    ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);
    //
    //    ((TH1F*)gDirectory->FindObject(h_angle_e1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()) );

    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"e1e2_hod_Tmatch");
    cutcounter++;




    /// **********  photon


    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~

    if(gclindx<0)     {cleanup();return 0;}
    //
    dirs[cutcounter-1]->cd();

    //    ((TH1F*)gDirectory->FindObject(h_el1trktime->GetName()))        ->Fill( el1trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trkhodtime->GetName()))        ->Fill( el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trktimediff->GetName()))    ->Fill( el1trktime - el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime1->GetName()))            ->Fill( el1cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e1trkcltimediff->GetName()))    ->Fill( el1trktime- el1cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_el2trktime->GetName()))        ->Fill( el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trkhodtime->GetName()))        ->Fill( el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trktimediff->GetName()))    ->Fill( el2trktime - el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime2->GetName()))            ->Fill( el2cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e2trkcltimediff->GetName()))    ->Fill( el2trktime-el2cltime );

    double gRadLKr = sqrt(gclx*gclx + gcly*gcly);

    Dgvtx_x = gclx - vertex[0] ;
    Dgvtx_y = gcly - vertex[1] ;
    Dgvtx_z = gclz - vertex[2] ;
    //    double Dgvtx_z = LKrz - vertex[2] ;
    gnorm = 1./sqrt(Dgvtx_x*Dgvtx_x + Dgvtx_y*Dgvtx_y + Dgvtx_z*Dgvtx_z);
    gpx = gclenergy * gnorm * Dgvtx_x ;
    gpy = gclenergy * gnorm * Dgvtx_y ;
    gpz = gclenergy * gnorm * Dgvtx_z ;

    glv = new TLorentzVector(gpx,gpy,gpz,gclenergy);

    ((TH1F*)gDirectory->FindObject(h_gE->GetName()))->Fill( (*glv).E() );
    ((TH1F*)gDirectory->FindObject(h_gP->GetName()))->Fill( (*glv).P() );
    ((TH1F*)gDirectory->FindObject(h_gPt->GetName()))->Fill((*glv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_angle_e1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_angle_ph_e1e2->GetName()))->Fill( (*eelv).Vect().Angle(glv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_a_e1e2_phe1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()), (*eelv).Vect().Angle(glv->Vect()) );

    gtrkDHC1x = vertex[0] + (DCHbz-vertex[2])/Dgvtx_z * Dgvtx_x;
    gtrkDHC1y = vertex[1] + (DCHbz-vertex[2])/Dgvtx_z * Dgvtx_y;
    gtrkRadDHC1 = sqrt(gtrkDHC1x*gtrkDHC1x + gtrkDHC1y*gtrkDHC1y);

    ((TH1F*)gDirectory->FindObject(h_gtrkRadDHC1->GetName()))->Fill( gtrkRadDHC1 );

    //
    //pi0 in lab frame
    pi0lv = new TLorentzVector( (*el1lv)+(*el2lv)+(*glv) );
    TVector3 pi0bv = pi0lv->BoostVector();
    el1lvcm = new TLorentzVector( *el1lv );
    el2lvcm = new TLorentzVector( *el2lv );
    glvcm   = new TLorentzVector( *glv );

    //
    //pi0 in its own cm frame
    el1lvcm->Boost(-pi0bv);
    el2lvcm->Boost(-pi0bv);
    glvcm->Boost(  -pi0bv);
    pi0lvcm = new TLorentzVector( (*el1lvcm)+(*el2lvcm)+(*glvcm) );

    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());// "ph_cluster");
    cutcounter++;



    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~

    if(gclenergy<3.)     {cleanup();return 0;}
    //
    dirs[cutcounter-1]->cd();

    //    ((TH1F*)gDirectory->FindObject(h_el1trktime->GetName()))        ->Fill( el1trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trkhodtime->GetName()))        ->Fill( el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trktimediff->GetName()))    ->Fill( el1trktime - el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime1->GetName()))            ->Fill( el1cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e1trkcltimediff->GetName()))    ->Fill( el1trktime- el1cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_el2trktime->GetName()))        ->Fill( el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trkhodtime->GetName()))        ->Fill( el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trktimediff->GetName()))    ->Fill( el2trktime - el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime2->GetName()))            ->Fill( el2cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e2trkcltimediff->GetName()))    ->Fill( el2trktime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_gcltime->GetName()))->Fill( gcltime );
    //    ((TH1F*)gDirectory->FindObject(h_e1gtimediff->GetName()))->Fill( el1cltime-gcltime );
    //    ((TH1F*)gDirectory->FindObject(h_e2gtimediff->GetName()))->Fill( el2cltime-gcltime );

    ((TH1F*)gDirectory->FindObject(h_gE->GetName()))->Fill( (*glv).E() );
    ((TH1F*)gDirectory->FindObject(h_gP->GetName()))->Fill( (*glv).P() );
    ((TH1F*)gDirectory->FindObject(h_gPt->GetName()))->Fill( (*glv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_angle_e1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_angle_ph_e1e2->GetName()))->Fill( (*eelv).Vect().Angle(glv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_a_e1e2_phe1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()), (*eelv).Vect().Angle(glv->Vect()) );

    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"ph_energy");
    cutcounter++;



    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~

    if(gtrkRadDHC1<14.)     {cleanup();return 0;}
    //
    dirs[cutcounter-1]->cd();

    //    ((TH1F*)gDirectory->FindObject(h_el1trktime->GetName()))        ->Fill( el1trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trkhodtime->GetName()))        ->Fill( el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trktimediff->GetName()))    ->Fill( el1trktime - el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime1->GetName()))            ->Fill( el1cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e1trkcltimediff->GetName()))    ->Fill( el1trktime- el1cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_el2trktime->GetName()))        ->Fill( el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trkhodtime->GetName()))        ->Fill( el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trktimediff->GetName()))    ->Fill( el2trktime - el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime2->GetName()))            ->Fill( el2cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e2trkcltimediff->GetName()))    ->Fill( el2trktime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_gcltime->GetName()))->Fill( gcltime );
    //    ((TH1F*)gDirectory->FindObject(h_e1gtimediff->GetName()))->Fill( el1cltime-gcltime );
    //    ((TH1F*)gDirectory->FindObject(h_e2gtimediff->GetName()))->Fill( el2cltime-gcltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_gE->GetName()))->Fill( (*glv).E() );
    //    ((TH1F*)gDirectory->FindObject(h_gP->GetName()))->Fill( (*glv).P() );
    //    ((TH1F*)gDirectory->FindObject(h_gPt->GetName()))->Fill( (*glv).Pt() );
    //
    //    ((TH1F*)gDirectory->FindObject(h_angle_e1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_angle_ph_e1e2->GetName()))->Fill( (*eelv).Vect().Angle(glv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_a_e1e2_phe1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()), (*eelv).Vect().Angle(glv->Vect()) );

    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"ph_RatDHC1");
    cutcounter++;




    /// **********  charged pion


    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~

    if (pchtrkindx<0)    {cleanup();return 0;}
    //
    dirs[cutcounter-1]->cd();

    //    ((TH1F*)gDirectory->FindObject(h_el1trktime->GetName()))        ->Fill( el1trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trkhodtime->GetName()))        ->Fill( el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trktimediff->GetName()))    ->Fill( el1trktime - el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime1->GetName()))            ->Fill( el1cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e1trkcltimediff->GetName()))    ->Fill( el1trktime- el1cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_el2trktime->GetName()))        ->Fill( el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trkhodtime->GetName()))        ->Fill( el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trktimediff->GetName()))    ->Fill( el2trktime - el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime2->GetName()))            ->Fill( el2cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e2trkcltimediff->GetName()))    ->Fill( el2trktime-el2cltime );


    pchtrktime         = sevt->track[pchtrkindx].time;
    pchcltime         = sevt->cluster[ sevt->track[pchtrkindx].iClus ].time;
    pchtrkhodtime     = sevt->track[pchtrkindx].hodTime;
    icl_track       = sevt->track[pchtrkindx].iClus;

    ((TH1F*)gDirectory->FindObject(h_pchtrktime->GetName()))        ->Fill(pchtrktime);
    ((TH1F*)gDirectory->FindObject(h_pchtrkhodtime->GetName()))        ->Fill( pchtrkhodtime );
    ((TH1F*)gDirectory->FindObject(h_pchcltime->GetName()))            ->Fill(pchcltime);
    ((TH1F*)gDirectory->FindObject(h_pchtrktimediff->GetName()))    ->Fill( pchtrktime-pchtrkhodtime );
    ((TH1F*)gDirectory->FindObject(h_pchtrkcltimediff->GetName()))    ->Fill( pchtrktime-pchcltime );

    //    cout<< "    pppion "<<pchtrkindx<<"    icl_track="<<icl_track<<"    Dt="<<pchtrktime-pchtrkhodtime<<endl;

    pchv_uncorr[0] = sevt->track[pchtrkindx].bdxdz;
    pchv_uncorr[1] = sevt->track[pchtrkindx].bdydz;
    pchv_uncorr[2] = 1.;
    pchp_uncorr[0] = sevt->track[pchtrkindx].bx;
    pchp_uncorr[1] = sevt->track[pchtrkindx].by;
    pchp_uncorr[2] = DCHbz;

    norm1 = 1./sqrt(pchv_uncorr[0]*pchv_uncorr[0] + pchv_uncorr[1]*pchv_uncorr[1] + pchv_uncorr[2]*pchv_uncorr[2] );
    p1x = norm1 * pchv_uncorr[0] * sevt->track[pchtrkindx].p ;
    p1y = norm1 * pchv_uncorr[1] * sevt->track[pchtrkindx].p ;
    p1z = norm1 *                  sevt->track[pchtrkindx].p ;
    pchE = sqrt( p1x*p1x+p1y*p1y+p1z*p1z + massPionC*massPionC);

    pchlv = new TLorentzVector(p1x, p1y, p1z, pchE);


    // find vertex of pi and ee
    pchv[0] = trkCorrSlopesV[pchtrkindx]->x();
    pchv[1] = trkCorrSlopesV[pchtrkindx]->y();
    pchv[2] = 1.;
    pchp[0] = trkCorrMidPointsV[pchtrkindx]->x();
    pchp[1] = trkCorrMidPointsV[pchtrkindx]->y();
    pchp[2] = trkCorrMidPointsV[pchtrkindx]->z();

    if ( closap_double_(elp1, pchp, elv1, pchv, &cda_e1pch, vertex_e1pch) )
        pchee_vertex_OK = TRUE;
    else
        pchee_vertex_OK = FALSE;

    ((TH1F*)gDirectory->FindObject(h_pchlvE->GetName()))->Fill( (*pchlv).E() );
    ((TH1F*)gDirectory->FindObject(h_pchlvP->GetName()))->Fill( (*pchlv).P() );
    ((TH1F*)gDirectory->FindObject(h_pchlvPt->GetName()))->Fill((*pchlv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_angle_e1e2->GetName()))->Fill(     (*el1lv).Vect().Angle(el2lv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_angle_ph_e1e2->GetName()))->Fill(  (*eelv).Vect().Angle(glv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_angle_pch_pi0->GetName()))->Fill(  (*pi0lv).Vect().Angle(pchlv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_a_e1e2_phe1e2->GetName()))->Fill(  (*el1lv).Vect().Angle(el2lv->Vect()),     (*eelv).Vect().Angle(glv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_a_e1e2_pchpi0->GetName()))->Fill(     (*el1lv).Vect().Angle(el2lv->Vect()),     (*pi0lv).Vect().Angle(pchlv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_a_phe1e2_pchpi0->GetName()))->Fill((*eelv).Vect().Angle(glv->Vect()),         (*pi0lv).Vect().Angle(pchlv->Vect()) );

    ((TH1F*)gDirectory->FindObject(h_bx_pch->GetName()))            ->Fill(pchp_uncorr[0]);
    ((TH1F*)gDirectory->FindObject(h_by_pch->GetName()))            ->Fill(pchp_uncorr[1]);
    ((TH1F*)gDirectory->FindObject(h_bdxdz_pch->GetName()))            ->Fill(pchv_uncorr[0]);
    ((TH1F*)gDirectory->FindObject(h_bdydz_pch->GetName()))            ->Fill(pchv_uncorr[1]) ;
    ((TH1F*)gDirectory->FindObject(h_bx_vs_by_pch->GetName()))        ->Fill(pchp_uncorr[0], pchp_uncorr[1]);
    ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_pch->GetName()))->Fill(pchv_uncorr[0], pchv_uncorr[1]);

    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"pi_track");
    cutcounter++;




    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~

    if(IS_DATA)    {
        //after final selection this difference has a maximum at -5.5 and sigma 1.5
        if( fabs(pchtrktime-pchtrkhodtime) > 12. )    return 0;  //5 is from Manuel's thesis
    }
    else    {
        if( pchtrkhodtime<135. || pchtrkhodtime>148. )    return 0;
    }
    //
    dirs[cutcounter-1]->cd();

    //    ((TH1F*)gDirectory->FindObject(h_el1trktime->GetName()))        ->Fill( el1trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trkhodtime->GetName()))        ->Fill( el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trktimediff->GetName()))    ->Fill( el1trktime - el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime1->GetName()))            ->Fill( el1cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e1trkcltimediff->GetName()))    ->Fill( el1trktime- el1cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_el2trktime->GetName()))        ->Fill( el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trkhodtime->GetName()))        ->Fill( el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trktimediff->GetName()))    ->Fill( el2trktime - el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime2->GetName()))            ->Fill( el2cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e2trkcltimediff->GetName()))    ->Fill( el2trktime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_e1e2timediff->GetName()))        ->Fill( el1trktime-el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_e1e2cltimediff->GetName()))    ->Fill( el1cltime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_pchtrktime->GetName()))        ->Fill(pchtrktime);
    //    ((TH1F*)gDirectory->FindObject(h_pchtrkhodtime->GetName()))        ->Fill( pchtrkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_pchcltime->GetName()))            ->Fill(pchcltime);
    //    ((TH1F*)gDirectory->FindObject(h_pchtrktimediff->GetName()))    ->Fill( pchtrktime-pchtrkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrkcltimediff->GetName()))    ->Fill( pchtrktime-pchcltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_pchlvE->GetName()))->Fill( (*pchlv).E() );
    //    ((TH1F*)gDirectory->FindObject(h_pchlvP->GetName()))->Fill( (*pchlv).P() );
    //    ((TH1F*)gDirectory->FindObject(h_pchlvPt->GetName()))->Fill( (*pchlv).Pt() );
    //
    //    ((TH1F*)gDirectory->FindObject(h_angle_e1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_angle_ph_e1e2->GetName()))->Fill( (*eelv).Vect().Angle(glv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_angle_pch_pi0->GetName()))->Fill( (*pi0lv).Vect().Angle(pchlv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_a_e1e2_phe1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()), (*eelv).Vect().Angle(glv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_a_e1e2_pchpi0->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()), (*pi0lv).Vect().Angle(pchlv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_a_phe1e2_pchpi0->GetName()))->Fill( (*eelv).Vect().Angle(glv->Vect()), (*pi0lv).Vect().Angle(pchlv->Vect()) );

    ((TH1F*)gDirectory->FindObject(h_bx_pch->GetName()))            ->Fill(pchp_uncorr[0]);
    ((TH1F*)gDirectory->FindObject(h_by_pch->GetName()))            ->Fill(pchp_uncorr[1]);
    ((TH1F*)gDirectory->FindObject(h_bdxdz_pch->GetName()))            ->Fill(pchv_uncorr[0]);
    ((TH1F*)gDirectory->FindObject(h_bdydz_pch->GetName()))            ->Fill(pchv_uncorr[1]) ;
    ((TH1F*)gDirectory->FindObject(h_bx_vs_by_pch->GetName()))        ->Fill(pchp_uncorr[0], pchp_uncorr[1]);
    ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_pch->GetName()))->Fill(pchv_uncorr[0], pchv_uncorr[1]);

    file1->cd();
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"pi_hod_Tmatch");
    h_cutflow->Fill(cutcounter);
    cutcounter++;




    //------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------
    //                            electron and positron matching
    //------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------


    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~

    //    if( fabs((long double)(el1cltime-el2cltime))>15. )        {cleanup();return 0;}
    if( fabs((long double)(el1cltime-el2cltime))>10. )        {cleanup();return 0;}
    //
    dirs[cutcounter-1]->cd();

    //    cout<<"trkCorrSlopesV[eltrkindx1]->x()    = "<<trkCorrSlopesV[eltrkindx1]->x()     <<"        sevt->track[eltrkindx1].bdxdz="<<sevt->track[eltrkindx1].bdxdz<<endl;
    //    cout<<"trkCorrMidPointsV[eltrkindx1]->x() = "<<trkCorrMidPointsV[eltrkindx1]->x()<<"        sevt->track[eltrkindx1].bx   ="<<sevt->track[eltrkindx1].bx<<endl;

    //    ((TH1F*)gDirectory->FindObject(h_el1trktime->GetName()))        ->Fill( el1trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trkhodtime->GetName()))        ->Fill( el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trktimediff->GetName()))    ->Fill( el1trktime - el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime1->GetName()))            ->Fill( el1cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e1trkcltimediff->GetName()))    ->Fill( el1trktime- el1cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_el2trktime->GetName()))        ->Fill( el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trkhodtime->GetName()))        ->Fill( el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trktimediff->GetName()))    ->Fill( el2trktime - el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime2->GetName()))            ->Fill( el2cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e2trkcltimediff->GetName()))    ->Fill( el2trktime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_e1e2timediff->GetName()))        ->Fill( el1trktime-el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_e1e2cltimediff->GetName()))    ->Fill( el1cltime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_pchtrktime->GetName()))        ->Fill( pchtrktime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrkhodtime->GetName()))        ->Fill( pchtrkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_pchcltime->GetName()))            ->Fill( pchcltime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrktimediff->GetName()))    ->Fill( pchtrktime-pchtrkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrkcltimediff->GetName()))    ->Fill( pchtrktime-pchcltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill(elp1_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill(elp1_uncorr[1]);
    //    ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill(elv1_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill(elv1_uncorr[1]);
    //    ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill(elp1_uncorr[0], elp1_uncorr[1] );
    //    ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill(elv1_uncorr[0], elv1_uncorr[1] );
    //
    //    ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( el1lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( el2lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el1lv->E() );
    //    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el2lv->E() );
    //    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el1lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el2lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el1lv->Pt() );
    //    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el2lv->Pt() );
    //
    //    ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
    ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);

    //    ((TH1F*)gDirectory->FindObject(h_gE->GetName()))->Fill( (*glv).E() );
    //    ((TH1F*)gDirectory->FindObject(h_gP->GetName()))->Fill( (*glv).P() );
    //    ((TH1F*)gDirectory->FindObject(h_gPt->GetName()))->Fill( (*glv).Pt() );
    //
    //    ((TH1F*)gDirectory->FindObject(h_pchlvE->GetName()))->Fill( (*pchlv).E() );
    //    ((TH1F*)gDirectory->FindObject(h_pchlvP->GetName()))->Fill( (*pchlv).P() );
    //    ((TH1F*)gDirectory->FindObject(h_pchlvPt->GetName()))->Fill( (*pchlv).Pt() );
    //
    //    ((TH1F*)gDirectory->FindObject(h_angle_e1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_angle_ph_e1e2->GetName()))->Fill( (*eelv).Vect().Angle(glv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_angle_pch_pi0->GetName()))->Fill( (*pi0lv).Vect().Angle(pchlv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_a_e1e2_phe1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()), (*eelv).Vect().Angle(glv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_a_e1e2_pchpi0->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()), (*pi0lv).Vect().Angle(pchlv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_a_phe1e2_pchpi0->GetName()))->Fill( (*eelv).Vect().Angle(glv->Vect()), (*pi0lv).Vect().Angle(pchlv->Vect()) );

    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"e1e2_cl_Tmatch");
    cutcounter++;



    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~

    if( pow(elp1_uncorr[0],2.) + pow(elp1_uncorr[1],2.) < 196. )    {cleanup();return 0;} //14*14 /it was 12*12/
    if( pow(elp2_uncorr[0],2.) + pow(elp2_uncorr[1],2.) < 196. )    {cleanup();return 0;}
    //
    dirs[cutcounter-1]->cd();

    //    ((TH1F*)gDirectory->FindObject(h_el1trktime->GetName()))        ->Fill( el1trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trkhodtime->GetName()))        ->Fill( el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trktimediff->GetName()))    ->Fill( el1trktime - el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime1->GetName()))            ->Fill( el1cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e1trkcltimediff->GetName()))    ->Fill( el1trktime- el1cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_el2trktime->GetName()))        ->Fill( el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trkhodtime->GetName()))        ->Fill( el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trktimediff->GetName()))    ->Fill( el2trktime - el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime2->GetName()))            ->Fill( el2cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e2trkcltimediff->GetName()))    ->Fill( el2trktime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_e1e2timediff->GetName()))        ->Fill( el1trktime-el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_e1e2cltimediff->GetName()))    ->Fill( el1cltime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_pchtrktime->GetName()))        ->Fill( pchtrktime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrkhodtime->GetName()))        ->Fill( pchtrkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_pchcltime->GetName()))            ->Fill( pchcltime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrktimediff->GetName()))    ->Fill( pchtrktime-pchtrkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrkcltimediff->GetName()))    ->Fill( pchtrktime-pchcltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill(elp1_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill(elp1_uncorr[1]);
    //    ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill(elv1_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill(elv1_uncorr[1]);
    //    ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill(elp1_uncorr[0], elp1_uncorr[1] );
    //    ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill(elv1_uncorr[0], elv1_uncorr[1] );
    //
    //    ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( el1lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( el2lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el1lv->E() );
    //    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el2lv->E() );
    //    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el1lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el2lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el1lv->Pt() );
    //    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el2lv->Pt() );
    //
    //    ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
    ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);

    //    ((TH1F*)gDirectory->FindObject(h_gE->GetName()))->Fill( (*glv).E() );
    //    ((TH1F*)gDirectory->FindObject(h_gP->GetName()))->Fill( (*glv).P() );
    //    ((TH1F*)gDirectory->FindObject(h_gPt->GetName()))->Fill( (*glv).Pt() );
    //
    //    ((TH1F*)gDirectory->FindObject(h_pchlvE->GetName()))->Fill( (*pchlv).E() );
    //    ((TH1F*)gDirectory->FindObject(h_pchlvP->GetName()))->Fill( (*pchlv).P() );
    //    ((TH1F*)gDirectory->FindObject(h_pchlvPt->GetName()))->Fill( (*pchlv).Pt() );
    //
    //    ((TH1F*)gDirectory->FindObject(h_angle_e1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_angle_ph_e1e2->GetName()))->Fill( (*eelv).Vect().Angle(glv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_angle_pch_pi0->GetName()))->Fill( (*pi0lv).Vect().Angle(pchlv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_a_e1e2_phe1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()), (*eelv).Vect().Angle(glv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_a_e1e2_pchpi0->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()), (*pi0lv).Vect().Angle(pchlv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_a_phe1e2_pchpi0->GetName()))->Fill( (*eelv).Vect().Angle(glv->Vect()), (*pi0lv).Vect().Angle(pchlv->Vect()) );

    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"e1e2_DCHinnerR");
    cutcounter++;


    //    cout<<"2"<<endl;


    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~

    //    if( pow(elp1_uncorr[0],2.) + pow(elp1_uncorr[1],2.) > 13225. )    {cleanup();return 0;} //115*115
    //    if( pow(elp2_uncorr[0],2.) + pow(elp2_uncorr[1],2.) > 13225. )    {cleanup();return 0;}
    if( pow(elp1_uncorr[0],2.) + pow(elp1_uncorr[1],2.) > 12100. )    {cleanup();return 0;} //110*110
    if( pow(elp2_uncorr[0],2.) + pow(elp2_uncorr[1],2.) > 12100. )    {cleanup();return 0;}
    //
    dirs[cutcounter-1]->cd();

    //    ((TH1F*)gDirectory->FindObject(h_el1trktime->GetName()))        ->Fill( el1trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trkhodtime->GetName()))        ->Fill( el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trktimediff->GetName()))    ->Fill( el1trktime - el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime1->GetName()))            ->Fill( el1cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e1trkcltimediff->GetName()))    ->Fill( el1trktime- el1cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_el2trktime->GetName()))        ->Fill( el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trkhodtime->GetName()))        ->Fill( el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trktimediff->GetName()))    ->Fill( el2trktime - el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime2->GetName()))            ->Fill( el2cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e2trkcltimediff->GetName()))    ->Fill( el2trktime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_e1e2timediff->GetName()))        ->Fill( el1trktime-el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_e1e2cltimediff->GetName()))    ->Fill( el1cltime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_pchtrktime->GetName()))        ->Fill( pchtrktime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrkhodtime->GetName()))        ->Fill( pchtrkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_pchcltime->GetName()))            ->Fill( pchcltime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrktimediff->GetName()))    ->Fill( pchtrktime-pchtrkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrkcltimediff->GetName()))    ->Fill( pchtrktime-pchcltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill(elp1_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill(elp1_uncorr[1]);
    //    ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill(elv1_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill(elv1_uncorr[1]);
    //    ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill(elp1_uncorr[0], elp1_uncorr[1] );
    //    ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill(elv1_uncorr[0], elv1_uncorr[1] );
    //
    //    ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( el1lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( el2lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el1lv->E() );
    //    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el2lv->E() );
    //    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el1lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el2lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el1lv->Pt() );
    //    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el2lv->Pt() );
    //
    //    ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
    ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);

    //    ((TH1F*)gDirectory->FindObject(h_gE->GetName()))->Fill( (*glv).E() );
    //    ((TH1F*)gDirectory->FindObject(h_gP->GetName()))->Fill( (*glv).P() );
    //    ((TH1F*)gDirectory->FindObject(h_gPt->GetName()))->Fill( (*glv).Pt() );
    //
    //    ((TH1F*)gDirectory->FindObject(h_pchlvE->GetName()))->Fill( (*pchlv).E() );
    //    ((TH1F*)gDirectory->FindObject(h_pchlvP->GetName()))->Fill( (*pchlv).P() );
    //    ((TH1F*)gDirectory->FindObject(h_pchlvPt->GetName()))->Fill( (*pchlv).Pt() );
    //
    //    ((TH1F*)gDirectory->FindObject(h_angle_e1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_angle_ph_e1e2->GetName()))->Fill( (*eelv).Vect().Angle(glv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_angle_pch_pi0->GetName()))->Fill( (*pi0lv).Vect().Angle(pchlv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_a_e1e2_phe1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()), (*eelv).Vect().Angle(glv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_a_e1e2_pchpi0->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()), (*pi0lv).Vect().Angle(pchlv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_a_phe1e2_pchpi0->GetName()))->Fill( (*eelv).Vect().Angle(glv->Vect()), (*pi0lv).Vect().Angle(pchlv->Vect()) );

    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"e1e2_DCHouterR");
    cutcounter++;


    //    cout<<"3"<<endl;


    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~

    if(vertex_OK!=1)    return 0;
    //
    dirs[cutcounter-1]->cd();

    //    ((TH1F*)gDirectory->FindObject(h_el1trktime->GetName()))        ->Fill( el1trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trkhodtime->GetName()))        ->Fill( el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trktimediff->GetName()))    ->Fill( el1trktime - el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime1->GetName()))            ->Fill( el1cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e1trkcltimediff->GetName()))    ->Fill( el1trktime- el1cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_el2trktime->GetName()))        ->Fill( el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trkhodtime->GetName()))        ->Fill( el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trktimediff->GetName()))    ->Fill( el2trktime - el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime2->GetName()))            ->Fill( el2cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e2trkcltimediff->GetName()))    ->Fill( el2trktime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_e1e2timediff->GetName()))        ->Fill( el1trktime-el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_e1e2cltimediff->GetName()))    ->Fill( el1cltime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_pchtrktime->GetName()))        ->Fill( pchtrktime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrkhodtime->GetName()))        ->Fill( pchtrkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_pchcltime->GetName()))            ->Fill( pchcltime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrktimediff->GetName()))    ->Fill( pchtrktime-pchtrkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrkcltimediff->GetName()))    ->Fill( pchtrktime-pchcltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill(elp1_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill(elp1_uncorr[1]);
    //    ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill(elv1_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill(elv1_uncorr[1]);
    //    ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill(elp1_uncorr[0], elp1_uncorr[1] );
    //    ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill(elv1_uncorr[0], elv1_uncorr[1] );
    //
    //    ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( el1lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( el2lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el1lv->E() );
    //    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el2lv->E() );
    //    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el1lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el2lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el1lv->Pt() );
    //    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el2lv->Pt() );
    //
    //    ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
    ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);
    //
    //    ((TH1F*)gDirectory->FindObject(h_gE->GetName()))->Fill( (*glv).E() );
    //    ((TH1F*)gDirectory->FindObject(h_gP->GetName()))->Fill( (*glv).P() );
    //    ((TH1F*)gDirectory->FindObject(h_gPt->GetName()))->Fill( (*glv).Pt() );
    //
    //    ((TH1F*)gDirectory->FindObject(h_pchlvE->GetName()))->Fill( (*pchlv).E() );
    //    ((TH1F*)gDirectory->FindObject(h_pchlvP->GetName()))->Fill( (*pchlv).P() );
    //    ((TH1F*)gDirectory->FindObject(h_pchlvPt->GetName()))->Fill( (*pchlv).Pt() );
    //
    //    ((TH1F*)gDirectory->FindObject(h_angle_e1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_angle_ph_e1e2->GetName()))->Fill( (*eelv).Vect().Angle(glv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_angle_pch_pi0->GetName()))->Fill( (*pi0lv).Vect().Angle(pchlv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_a_e1e2_phe1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()), (*eelv).Vect().Angle(glv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_a_e1e2_pchpi0->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()), (*pi0lv).Vect().Angle(pchlv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_a_phe1e2_pchpi0->GetName()))->Fill( (*eelv).Vect().Angle(glv->Vect()), (*pi0lv).Vect().Angle(pchlv->Vect()) );
    //
    ((TH1F*)gDirectory->FindObject(h_cda_ee->GetName()))->Fill(cda);
    ((TH1F*)gDirectory->FindObject(h_vx->GetName()))->Fill(vertex[0]);
    ((TH1F*)gDirectory->FindObject(h_vy->GetName()))->Fill(vertex[1]);
    ((TH1F*)gDirectory->FindObject(h_vz->GetName()))->Fill(vertex[2]);
    ((TH2F*)gDirectory->FindObject(h_vx_vs_vy->GetName()))->Fill(vertex[0],vertex[1]);

    ((TH1F*)gDirectory->FindObject(h_vx_e1pch->GetName()))->Fill(vertex_e1pch[0]);
    ((TH1F*)gDirectory->FindObject(h_vy_e1pch->GetName()))->Fill(vertex_e1pch[1]);
    ((TH1F*)gDirectory->FindObject(h_vz_e1pch->GetName()))->Fill(vertex_e1pch[2]);
    ((TH1F*)gDirectory->FindObject(h_cda_e1pch->GetName()))->Fill(cda_e1pch);
    ((TH1F*)gDirectory->FindObject(h_vz_e1pch_diff->GetName()))->Fill( vertex[2] - vertex_e1pch[2] );

    ((TH1F*)gDirectory->FindObject(h_cda_e1pch->GetName()))->Fill(cda_e1pch);
    ((TH2F*)gDirectory->FindObject(h_vx_vs_vy_e1pch->GetName()))->Fill(vertex_e1pch[0],vertex_e1pch[1]);

    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"e1e2_vtx_found");
    cutcounter++;




    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~

    //    if( vertex[2]<-2000. || vertex[2]>8000. )    {cleanup();return 0;}
    if( vertex[2]<-1500. || vertex[2]>7000. )    {cleanup();return 0;}
    //
    dirs[cutcounter-1]->cd();

    //    ((TH1F*)gDirectory->FindObject(h_el1trktime->GetName()))        ->Fill( el1trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trkhodtime->GetName()))        ->Fill( el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trktimediff->GetName()))    ->Fill( el1trktime - el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime1->GetName()))            ->Fill( el1cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e1trkcltimediff->GetName()))    ->Fill( el1trktime- el1cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_el2trktime->GetName()))        ->Fill( el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trkhodtime->GetName()))        ->Fill( el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trktimediff->GetName()))    ->Fill( el2trktime - el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime2->GetName()))            ->Fill( el2cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e2trkcltimediff->GetName()))    ->Fill( el2trktime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_e1e2timediff->GetName()))        ->Fill( el1trktime-el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_e1e2cltimediff->GetName()))    ->Fill( el1cltime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_pchtrktime->GetName()))        ->Fill( pchtrktime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrkhodtime->GetName()))        ->Fill( pchtrkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_pchcltime->GetName()))            ->Fill( pchcltime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrktimediff->GetName()))    ->Fill( pchtrktime-pchtrkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrkcltimediff->GetName()))    ->Fill( pchtrktime-pchcltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill(elp1_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill(elp1_uncorr[1]);
    //    ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill(elv1_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill(elv1_uncorr[1]);
    //    ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill(elp1_uncorr[0], elp1_uncorr[1] );
    //    ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill(elv1_uncorr[0], elv1_uncorr[1] );
    //
    //    ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( el1lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( el2lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el1lv->E() );
    //    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el2lv->E() );
    //    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el1lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el2lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el1lv->Pt() );
    //    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el2lv->Pt() );
    //
    //    ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
    ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);

    ((TH1F*)gDirectory->FindObject(h_gE->GetName()))->Fill( (*glv).E() );
    ((TH1F*)gDirectory->FindObject(h_gP->GetName()))->Fill( (*glv).P() );
    ((TH1F*)gDirectory->FindObject(h_gPt->GetName()))->Fill( (*glv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_pchlvE->GetName()))->Fill( (*pchlv).E() );
    ((TH1F*)gDirectory->FindObject(h_pchlvP->GetName()))->Fill( (*pchlv).P() );
    ((TH1F*)gDirectory->FindObject(h_pchlvPt->GetName()))->Fill( (*pchlv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_angle_e1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_angle_ph_e1e2->GetName()))->Fill( (*eelv).Vect().Angle(glv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_angle_pch_pi0->GetName()))->Fill( (*pi0lv).Vect().Angle(pchlv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_a_e1e2_phe1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()), (*eelv).Vect().Angle(glv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_a_e1e2_pchpi0->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()), (*pi0lv).Vect().Angle(pchlv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_a_phe1e2_pchpi0->GetName()))->Fill( (*eelv).Vect().Angle(glv->Vect()), (*pi0lv).Vect().Angle(pchlv->Vect()) );

    ((TH1F*)gDirectory->FindObject(h_cda_ee->GetName()))->Fill(cda);
    ((TH1F*)gDirectory->FindObject(h_vx->GetName()))->Fill(vertex[0]);
    ((TH1F*)gDirectory->FindObject(h_vy->GetName()))->Fill(vertex[1]);
    ((TH1F*)gDirectory->FindObject(h_vz->GetName()))->Fill(vertex[2]);
    ((TH2F*)gDirectory->FindObject(h_vx_vs_vy->GetName()))->Fill(vertex[0],vertex[1]);

    ((TH1F*)gDirectory->FindObject(h_vx_e1pch->GetName()))->Fill(vertex_e1pch[0]);
    ((TH1F*)gDirectory->FindObject(h_vy_e1pch->GetName()))->Fill(vertex_e1pch[1]);
    ((TH1F*)gDirectory->FindObject(h_vz_e1pch->GetName()))->Fill(vertex_e1pch[2]);
    ((TH1F*)gDirectory->FindObject(h_cda_e1pch->GetName()))->Fill(cda_e1pch);
    ((TH1F*)gDirectory->FindObject(h_vz_e1pch_diff->GetName()))->Fill( vertex[2] - vertex_e1pch[2] );

    ((TH1F*)gDirectory->FindObject(h_cda_e1pch->GetName()))->Fill(cda_e1pch);
    ((TH2F*)gDirectory->FindObject(h_vx_vs_vy_e1pch->GetName()))->Fill(vertex_e1pch[0],vertex_e1pch[1]);

    eldistDCH1plane = sqrt( pow((elp1_uncorr[0]-elp2_uncorr[0]),2.) + pow((elp1_uncorr[1]-elp2_uncorr[1]),2.) ) ;

    el1xLKrplane = sevt->track[eltrkindx1].x + (LKrz-DCHz) * sevt->track[eltrkindx1].dxdz ;
    el1yLKrplane = sevt->track[eltrkindx1].y + (LKrz-DCHz) * sevt->track[eltrkindx1].dydz ;
    el2xLKrplane = sevt->track[eltrkindx2].x + (LKrz-DCHz) * sevt->track[eltrkindx2].dxdz ;
    el2yLKrplane = sevt->track[eltrkindx2].y + (LKrz-DCHz) * sevt->track[eltrkindx2].dydz ;
    eldistLKrplane = sqrt( pow((el1xLKrplane-el2xLKrplane),2.) + pow((el1yLKrplane-el2yLKrplane),2.) ) ;

    ((TH1F*)gDirectory->FindObject(h_eldistDCH1plane->GetName()))->Fill( eldistDCH1plane );
    //    ((TH2F*)gDirectory->FindObject(h_eldistDCH1_ee_vz->GetName()))->Fill( eldistDCH1plane, vertex[2] );
    //    ((TH2F*)gDirectory->FindObject(h_eldistDCH1_pch_vz->GetName()))->Fill( eldistDCH1plane, vertex_e1pch[2] );

    ((TH1F*)gDirectory->FindObject(h_eldistLKrplane->GetName()))->Fill( eldistLKrplane );

    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"e1e2_vtx_good");
    cutcounter++;




    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~

    //    if(eldistDCH1plane<2.)        return 0;
    if(eldistDCH1plane<1.)        return 0;
    //
    dirs[cutcounter-1]->cd();

    //    ((TH1F*)gDirectory->FindObject(h_el1trktime->GetName()))        ->Fill( el1trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trkhodtime->GetName()))        ->Fill( el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trktimediff->GetName()))    ->Fill( el1trktime - el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime1->GetName()))            ->Fill( el1cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e1trkcltimediff->GetName()))    ->Fill( el1trktime- el1cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_el2trktime->GetName()))        ->Fill( el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trkhodtime->GetName()))        ->Fill( el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trktimediff->GetName()))    ->Fill( el2trktime - el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime2->GetName()))            ->Fill( el2cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e2trkcltimediff->GetName()))    ->Fill( el2trktime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_e1e2timediff->GetName()))        ->Fill( el1trktime-el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_e1e2cltimediff->GetName()))    ->Fill( el1cltime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_pchtrktime->GetName()))        ->Fill( pchtrktime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrkhodtime->GetName()))        ->Fill( pchtrkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_pchcltime->GetName()))            ->Fill( pchcltime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrktimediff->GetName()))    ->Fill( pchtrktime-pchtrkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrkcltimediff->GetName()))    ->Fill( pchtrktime-pchcltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill(elp1_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill(elp1_uncorr[1]);
    //    ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill(elv1_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill(elv1_uncorr[1]);
    //    ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill(elp1_uncorr[0], elp1_uncorr[1] );
    //    ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill(elv1_uncorr[0], elv1_uncorr[1] );
    //
    //    ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( el1lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( el2lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el1lv->E() );
    //    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el2lv->E() );
    //    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el1lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el2lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el1lv->Pt() );
    //    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el2lv->Pt() );
    //
    //    ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
    ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);

    //    ((TH1F*)gDirectory->FindObject(h_gE->GetName()))->Fill( (*glv).E() );
    //    ((TH1F*)gDirectory->FindObject(h_gP->GetName()))->Fill( (*glv).P() );
    //    ((TH1F*)gDirectory->FindObject(h_gPt->GetName()))->Fill( (*glv).Pt() );
    //
    //    ((TH1F*)gDirectory->FindObject(h_pchlvE->GetName()))->Fill( (*pchlv).E() );
    //    ((TH1F*)gDirectory->FindObject(h_pchlvP->GetName()))->Fill( (*pchlv).P() );
    //    ((TH1F*)gDirectory->FindObject(h_pchlvPt->GetName()))->Fill( (*pchlv).Pt() );
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
    //
    //    ((TH1F*)gDirectory->FindObject(h_eldistLKrplane->GetName()))->Fill( eldistLKrplane );

    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"e1e2_dist_DCH1");
    cutcounter++;


    //    cout<<"6"<<endl;


    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~

    //    if( eldistLKrplane<15. )        return 0;
    if( eldistLKrplane<20. )        {cleanup();return 0;}
    //
    dirs[cutcounter-1]->cd();

    //    ((TH1F*)gDirectory->FindObject(h_el1trktime->GetName()))        ->Fill( el1trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trkhodtime->GetName()))        ->Fill( el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trktimediff->GetName()))    ->Fill( el1trktime - el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime1->GetName()))            ->Fill( el1cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e1trkcltimediff->GetName()))    ->Fill( el1trktime- el1cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_el2trktime->GetName()))        ->Fill( el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trkhodtime->GetName()))        ->Fill( el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trktimediff->GetName()))    ->Fill( el2trktime - el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime2->GetName()))            ->Fill( el2cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e2trkcltimediff->GetName()))    ->Fill( el2trktime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_e1e2timediff->GetName()))        ->Fill( el1trktime-el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_e1e2cltimediff->GetName()))    ->Fill( el1cltime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_pchtrktime->GetName()))        ->Fill( pchtrktime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrkhodtime->GetName()))        ->Fill( pchtrkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_pchcltime->GetName()))            ->Fill( pchcltime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrktimediff->GetName()))    ->Fill( pchtrktime-pchtrkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrkcltimediff->GetName()))    ->Fill( pchtrktime-pchcltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill(elp1_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill(elp1_uncorr[1]);
    //    ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill(elv1_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill(elv1_uncorr[1]);
    //    ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill(elp1_uncorr[0], elp1_uncorr[1] );
    //    ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill(elv1_uncorr[0], elv1_uncorr[1] );
    //
    //    ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( el1lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( el2lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el1lv->E() );
    //    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el2lv->E() );
    //    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el1lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el2lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el1lv->Pt() );
    //    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el2lv->Pt() );
    //
    //    ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
    ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);

    //    ((TH1F*)gDirectory->FindObject(h_gE->GetName()))->Fill( (*glv).E() );
    //    ((TH1F*)gDirectory->FindObject(h_gP->GetName()))->Fill( (*glv).P() );
    //    ((TH1F*)gDirectory->FindObject(h_gPt->GetName()))->Fill( (*glv).Pt() );
    //
    //    ((TH1F*)gDirectory->FindObject(h_pchlvE->GetName()))->Fill( (*pchlv).E() );
    //    ((TH1F*)gDirectory->FindObject(h_pchlvP->GetName()))->Fill( (*pchlv).P() );
    //    ((TH1F*)gDirectory->FindObject(h_pchlvPt->GetName()))->Fill( (*pchlv).Pt() );
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
    //     ((TH1F*)gDirectory->FindObject(h_vz_e1pch_diff->GetName()))->Fill( vertex[2] - vertex_e1pch[2] );

     file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"e1e2_dist_LKr");
    cutcounter++;




    //------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------
    //                        photon and electron-positron matching
    //------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------


    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~

    if( fabs((long double)(el1cltime-gcltime)) > 10. )    {cleanup();return 0;}
    if( fabs((long double)(el2cltime-gcltime)) > 10. )    {cleanup();return 0;}
    //
    dirs[cutcounter-1]->cd();

    //    ((TH1F*)gDirectory->FindObject(h_el1trktime->GetName()))        ->Fill( el1trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trkhodtime->GetName()))        ->Fill( el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trktimediff->GetName()))    ->Fill( el1trktime - el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime1->GetName()))            ->Fill( el1cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e1trkcltimediff->GetName()))    ->Fill( el1trktime- el1cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_el2trktime->GetName()))        ->Fill( el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trkhodtime->GetName()))        ->Fill( el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trktimediff->GetName()))    ->Fill( el2trktime - el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime2->GetName()))            ->Fill( el2cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e2trkcltimediff->GetName()))    ->Fill( el2trktime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_e1e2timediff->GetName()))        ->Fill( el1trktime-el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_e1e2cltimediff->GetName()))    ->Fill( el1cltime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_pchtrktime->GetName()))        ->Fill( pchtrktime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrkhodtime->GetName()))        ->Fill( pchtrkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_pchcltime->GetName()))            ->Fill( pchcltime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrktimediff->GetName()))    ->Fill( pchtrktime-pchtrkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrkcltimediff->GetName()))    ->Fill( pchtrktime-pchcltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill(elp1_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill(elp1_uncorr[1]);
    //    ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill(elv1_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill(elv1_uncorr[1]);
    //    ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill(elp1_uncorr[0], elp1_uncorr[1] );
    //    ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill(elv1_uncorr[0], elv1_uncorr[1] );
    //
    //    ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( el1lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( el2lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el1lv->E() );
    //    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el2lv->E() );
    //    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el1lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el2lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el1lv->Pt() );
    //    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el2lv->Pt() );
    //
    //    ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
    ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);

    //    ((TH1F*)gDirectory->FindObject(h_gE->GetName()))->Fill( (*glv).E() );
    //    ((TH1F*)gDirectory->FindObject(h_gP->GetName()))->Fill( (*glv).P() );
    //    ((TH1F*)gDirectory->FindObject(h_gPt->GetName()))->Fill( (*glv).Pt() );
    //    ((TH1F*)gDirectory->FindObject(h_gtrkRadDHC1->GetName()))->Fill( gtrkRadDHC1 );
    //
    //
    //    ((TH1F*)gDirectory->FindObject(h_pchlvE->GetName()))->Fill( (*pchlv).E() );
    //    ((TH1F*)gDirectory->FindObject(h_pchlvP->GetName()))->Fill( (*pchlv).P() );
    //    ((TH1F*)gDirectory->FindObject(h_pchlvPt->GetName()))->Fill( (*pchlv).Pt() );
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

    trkxcl = sevt->track[eltrkindx1].x + (gclz-DCHz) * sevt->track[eltrkindx1].dxdz ;
    trkycl = sevt->track[eltrkindx1].y + (gclz-DCHz) * sevt->track[eltrkindx1].dydz ;
    gel1distLKrplane = sqrt( pow(trkxcl-gclx, 2.) + pow(trkycl-gcly, 2.) );
    trkxcl = sevt->track[eltrkindx2].x + (gclz-DCHz) * sevt->track[eltrkindx2].dxdz ;
    trkycl = sevt->track[eltrkindx2].y + (gclz-DCHz) * sevt->track[eltrkindx2].dydz ;
    gel2distLKrplane = sqrt( pow(trkxcl-gclx, 2.) + pow(trkycl-gcly, 2.) );

    //    gel1distLKrplane = sqrt( pow((el1xLKrplane-gclx),2.) + pow((el1yLKrplane-gcly),2.) ) ;
    //    gel2distLKrplane = sqrt( pow((el2xLKrplane-gclx),2.) + pow((el2yLKrplane-gcly),2.) ) ;

    ((TH1F*)gDirectory->FindObject(h_gel1distLKrplane->GetName()))->Fill( gel1distLKrplane );
    ((TH1F*)gDirectory->FindObject(h_gel2distLKrplane->GetName()))->Fill( gel2distLKrplane );

    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"ph_ee_Tmatch");
    cutcounter++;



    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~

    if(gel1distLKrplane<20.)    {cleanup();return 0;}
    if(gel2distLKrplane<20.)    {cleanup();return 0;}
    //
    dirs[cutcounter-1]->cd();

    //    ((TH1F*)gDirectory->FindObject(h_el1trktime->GetName()))        ->Fill( el1trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trkhodtime->GetName()))        ->Fill( el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trktimediff->GetName()))    ->Fill( el1trktime - el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime1->GetName()))            ->Fill( el1cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e1trkcltimediff->GetName()))    ->Fill( el1trktime- el1cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_el2trktime->GetName()))        ->Fill( el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trkhodtime->GetName()))        ->Fill( el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trktimediff->GetName()))    ->Fill( el2trktime - el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime2->GetName()))            ->Fill( el2cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e2trkcltimediff->GetName()))    ->Fill( el2trktime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_e1e2timediff->GetName()))        ->Fill( el1trktime-el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_e1e2cltimediff->GetName()))    ->Fill( el1cltime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_pchtrktime->GetName()))        ->Fill( pchtrktime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrkhodtime->GetName()))        ->Fill( pchtrkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_pchcltime->GetName()))            ->Fill( pchcltime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrktimediff->GetName()))    ->Fill( pchtrktime-pchtrkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrkcltimediff->GetName()))    ->Fill( pchtrktime-pchcltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill(elp1_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill(elp1_uncorr[1]);
    //    ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill(elv1_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill(elv1_uncorr[1]);
    //    ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill(elp1_uncorr[0], elp1_uncorr[1] );
    //    ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill(elv1_uncorr[0], elv1_uncorr[1] );
    //
    //    ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( el1lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( el2lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el1lv->E() );
    //    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el2lv->E() );
    //    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el1lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el2lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el1lv->Pt() );
    //    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el2lv->Pt() );
    //
    //    ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
    ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);

    //    ((TH1F*)gDirectory->FindObject(h_gE->GetName()))->Fill( (*glv).E() );
    //    ((TH1F*)gDirectory->FindObject(h_gP->GetName()))->Fill( (*glv).P() );
    //    ((TH1F*)gDirectory->FindObject(h_gPt->GetName()))->Fill( (*glv).Pt() );
    //    ((TH1F*)gDirectory->FindObject(h_gtrkRadDHC1->GetName()))->Fill( gtrkRadDHC1 );
    //
    //    ((TH1F*)gDirectory->FindObject(h_pchlvE->GetName()))->Fill( (*pchlv).E() );
    //    ((TH1F*)gDirectory->FindObject(h_pchlvP->GetName()))->Fill( (*pchlv).P() );
    //    ((TH1F*)gDirectory->FindObject(h_pchlvPt->GetName()))->Fill( (*pchlv).Pt() );
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

    //
    ((TH1F*)gDirectory->FindObject(h_pi0E->GetName()))->Fill(  pi0lv->E() );
    ((TH1F*)gDirectory->FindObject(h_pi0P->GetName()))->Fill(  pi0lv->P() );
    ((TH1F*)gDirectory->FindObject(h_pi0Pt->GetName()))->Fill( pi0lv->Pt() );
    ((TH1F*)gDirectory->FindObject(h_mofpi0->GetName()))->Fill( (*pi0lv).M() );

    ((TH1F*)gDirectory->FindObject(h_pi0Pcm->GetName()))->Fill( pi0lvcm->P() );
    ((TH1F*)gDirectory->FindObject(h_mofpi0diff->GetName()))->Fill( (*pi0lvcm).M()-(*pi0lv).M() );

    //
    //todo investigate what's wrong here
    double msq12 = 2 * pow(massElec, 2.) + 2* (*el1lvcm)*(*el2lvcm);
    double msq23 =     pow(massElec, 2.) + 2* (*el2lvcm)*(*glvcm);
    ((TH2F*)gDirectory->FindObject(h_pi0_dalitz->GetName()))->Fill(msq12, msq23);

    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"ph_ee_dist_LKr");
    cutcounter++;


    //    cout<<"8"<<endl;


    //------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------
    //                        charged pion and electron-positron matching
    //------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------


    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~

    if( fabs((long double)(el1trktime-pchtrktime)) > 10. )    {cleanup();return 0;}
    if( fabs((long double)(el1cltime-pchcltime))   > 10. )    {cleanup();return 0;}
    //
    dirs[cutcounter-1]->cd();

    //    ((TH1F*)gDirectory->FindObject(h_el1trktime->GetName()))        ->Fill( el1trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trkhodtime->GetName()))        ->Fill( el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trktimediff->GetName()))    ->Fill( el1trktime - el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime1->GetName()))            ->Fill( el1cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e1trkcltimediff->GetName()))    ->Fill( el1trktime- el1cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_el2trktime->GetName()))        ->Fill( el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trkhodtime->GetName()))        ->Fill( el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trktimediff->GetName()))    ->Fill( el2trktime - el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime2->GetName()))            ->Fill( el2cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e2trkcltimediff->GetName()))    ->Fill( el2trktime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_e1e2timediff->GetName()))        ->Fill( el1trktime-el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_e1e2cltimediff->GetName()))    ->Fill( el1cltime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_pchtrktime->GetName()))        ->Fill( pchtrktime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrkhodtime->GetName()))        ->Fill( pchtrkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_pchcltime->GetName()))            ->Fill( pchcltime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrktimediff->GetName()))    ->Fill( pchtrktime-pchtrkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrkcltimediff->GetName()))    ->Fill( pchtrktime-pchcltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill(elp1_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill(elp1_uncorr[1]);
    //    ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill(elv1_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill(elv1_uncorr[1]);
    //    ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill(elp1_uncorr[0], elp1_uncorr[1] );
    //    ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill(elv1_uncorr[0], elv1_uncorr[1] );
    //
    //    ((TH1F*)gDirectory->FindObject(h_bx_pch->GetName()))            ->Fill(pchp_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_by_pch->GetName()))            ->Fill(pchp_uncorr[1]);
    //    ((TH1F*)gDirectory->FindObject(h_bdxdz_pch->GetName()))            ->Fill(pchv_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_bdydz_pch->GetName()))            ->Fill(pchv_uncorr[1]) ;
    //    ((TH1F*)gDirectory->FindObject(h_bx_vs_by_pch->GetName()))        ->Fill(pchp_uncorr[0], pchp_uncorr[1]);
    //    ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_pch->GetName()))->Fill(pchv_uncorr[0], pchv_uncorr[1]);
    //
    //    ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( el1lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( el2lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el1lv->E() );
    //    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el2lv->E() );
    //    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el1lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el2lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el1lv->Pt() );
    //    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el2lv->Pt() );
    //
    //    ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
    ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);

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

    ((TH1F*)gDirectory->FindObject(h_pi0E->GetName()))->Fill(  pi0lv->E() );
    ((TH1F*)gDirectory->FindObject(h_pi0P->GetName()))->Fill(  pi0lv->P() );
    ((TH1F*)gDirectory->FindObject(h_pi0Pt->GetName()))->Fill( pi0lv->Pt() );
    ((TH1F*)gDirectory->FindObject(h_mofpi0->GetName()))->Fill( (*pi0lv).M() );

    ((TH1F*)gDirectory->FindObject(h_pi0Pcm->GetName()))->Fill( pi0lvcm->P() );
    ((TH1F*)gDirectory->FindObject(h_mofpi0diff->GetName()))->Fill( (*pi0lvcm).M()-(*pi0lv).M() );

    pche1distDHC1plane = sqrt( pow((elp1_uncorr[0] - sevt->track[pchtrkindx].x),2.) + pow((elp1_uncorr[1] - sevt->track[pchtrkindx].y),2.) ) ;
    pche2distDHC1plane = sqrt( pow((elp2_uncorr[0] - sevt->track[pchtrkindx].x),2.) + pow((elp2_uncorr[1] - sevt->track[pchtrkindx].y),2.) ) ;
    pchxLKrplane         = sevt->track[pchtrkindx].x + (LKrz-DCHz) * sevt->track[pchtrkindx].dxdz ;
    pchyLKrplane         = sevt->track[pchtrkindx].y + (LKrz-DCHz) * sevt->track[pchtrkindx].dydz ;
    pche1distLKrplane     = sqrt( pow((el1xLKrplane-pchxLKrplane),2.) + pow((el1yLKrplane-pchxLKrplane),2.) ) ;
    pche2distLKrplane     = sqrt( pow((el2xLKrplane-pchxLKrplane),2.) + pow((el2yLKrplane-pchxLKrplane),2.) ) ;
    gpchdistLKrplane    = sqrt( pow((pchxLKrplane-gclx),2.)         + pow((pchyLKrplane-gcly),2.) ) ;

    ((TH1F*)gDirectory->FindObject(h_pche1distDHC1plane->GetName()))->Fill( pche1distDHC1plane );
    ((TH1F*)gDirectory->FindObject(h_pche2distDHC1plane->GetName()))->Fill( pche2distDHC1plane );
    ((TH1F*)gDirectory->FindObject(h_pche1distLKrplane->GetName()))->Fill( pche1distLKrplane );
    ((TH1F*)gDirectory->FindObject(h_pche2distLKrplane->GetName()))->Fill( pche2distLKrplane );
    ((TH1F*)gDirectory->FindObject(h_gpchdistLKrplane->GetName()))->Fill( gpchdistLKrplane );

    pieelv = new TLorentzVector( (*el1lv)+(*el2lv)+(*pchlv) );

    ((TH1F*)gDirectory->FindObject(h_pieeE->GetName()))->Fill(pieelv->E());
    ((TH1F*)gDirectory->FindObject(h_pieeP->GetName()))->Fill(pieelv->P());
    ((TH1F*)gDirectory->FindObject(h_pieePt->GetName()))->Fill(pieelv->Pt());
    ((TH1F*)gDirectory->FindObject(h_pieeM->GetName()))->Fill(pieelv->M());

    file1->cd();
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"pi_ee_vtx_match")
    h_cutflow->Fill(cutcounter);
    cutcounter++;



    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~

    if( !pchee_vertex_OK)    {cleanup();return 0;}
    //    if( fabs((long double)(vertex[2] - vertex_e1pch[2]))  > 1500. )    {cleanup();return 0;} //old value 300
    if( vertex_e1pch[2]<-1500.  || vertex_e1pch[2]>7000. )    {cleanup();return 0;}
    //
    dirs[cutcounter-1]->cd();

    //    ((TH1F*)gDirectory->FindObject(h_el1trktime->GetName()))        ->Fill( el1trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trkhodtime->GetName()))        ->Fill( el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trktimediff->GetName()))    ->Fill( el1trktime - el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime1->GetName()))            ->Fill( el1cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e1trkcltimediff->GetName()))    ->Fill( el1trktime- el1cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_el2trktime->GetName()))        ->Fill( el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trkhodtime->GetName()))        ->Fill( el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trktimediff->GetName()))    ->Fill( el2trktime - el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime2->GetName()))            ->Fill( el2cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e2trkcltimediff->GetName()))    ->Fill( el2trktime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_e1e2timediff->GetName()))        ->Fill( el1trktime-el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_e1e2cltimediff->GetName()))    ->Fill( el1cltime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_pchtrktime->GetName()))        ->Fill( pchtrktime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrkhodtime->GetName()))        ->Fill( pchtrkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_pchcltime->GetName()))            ->Fill( pchcltime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrktimediff->GetName()))    ->Fill( pchtrktime-pchtrkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrkcltimediff->GetName()))    ->Fill( pchtrktime-pchcltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill(elp1_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill(elp1_uncorr[1]);
    //    ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill(elv1_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill(elv1_uncorr[1]);
    //    ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill(elp1_uncorr[0], elp1_uncorr[1] );
    //    ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill(elv1_uncorr[0], elv1_uncorr[1] );
    //
    //    ((TH1F*)gDirectory->FindObject(h_bx_pch->GetName()))            ->Fill(pchp_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_by_pch->GetName()))            ->Fill(pchp_uncorr[1]);
    //    ((TH1F*)gDirectory->FindObject(h_bdxdz_pch->GetName()))            ->Fill(pchv_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_bdydz_pch->GetName()))            ->Fill(pchv_uncorr[1]) ;
    //    ((TH1F*)gDirectory->FindObject(h_bx_vs_by_pch->GetName()))        ->Fill(pchp_uncorr[0], pchp_uncorr[1]);
    //    ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_pch->GetName()))->Fill(pchv_uncorr[0], pchv_uncorr[1]);
    //
    //    ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( el1lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( el2lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el1lv->E() );
    //    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el2lv->E() );
    //    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el1lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el2lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el1lv->Pt() );
    //    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el2lv->Pt() );
    //
    //    ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
    ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);

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

    ((TH1F*)gDirectory->FindObject(h_pi0E->GetName()))->Fill(  pi0lv->E() );
    ((TH1F*)gDirectory->FindObject(h_pi0P->GetName()))->Fill(  pi0lv->P() );
    ((TH1F*)gDirectory->FindObject(h_pi0Pt->GetName()))->Fill( pi0lv->Pt() );
    ((TH1F*)gDirectory->FindObject(h_mofpi0->GetName()))->Fill( (*pi0lv).M() );

    ((TH1F*)gDirectory->FindObject(h_pi0Pcm->GetName()))->Fill( pi0lvcm->P() );
    ((TH1F*)gDirectory->FindObject(h_mofpi0diff->GetName()))->Fill( (*pi0lvcm).M()-(*pi0lv).M() );

    //    ((TH1F*)gDirectory->FindObject(h_pche1distDHC1plane->GetName()))->Fill( pche1distDHC1plane );
    //    ((TH1F*)gDirectory->FindObject(h_pche2distDHC1plane->GetName()))->Fill( pche2distDHC1plane );
    //    ((TH1F*)gDirectory->FindObject(h_pche1distLKrplane->GetName()))->Fill( pche1distLKrplane );
    //    ((TH1F*)gDirectory->FindObject(h_pche2distLKrplane->GetName()))->Fill( pche2distLKrplane );
    //    ((TH1F*)gDirectory->FindObject(h_gpchdistLKrplane->GetName()))->Fill( gpchdistLKrplane );

    ((TH1F*)gDirectory->FindObject(h_pieeE->GetName()))->Fill(pieelv->E());
    ((TH1F*)gDirectory->FindObject(h_pieeP->GetName()))->Fill(pieelv->P());
    ((TH1F*)gDirectory->FindObject(h_pieePt->GetName()))->Fill(pieelv->Pt());
    ((TH1F*)gDirectory->FindObject(h_pieeM->GetName()))->Fill(pieelv->M());

    file1->cd();
    //    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"pi_ee_vtx_match")
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"pi_vtx_good")
    h_cutflow->Fill(cutcounter);
    cutcounter++;



    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~

    if(pche1distLKrplane<30.)    {cleanup();return 0;}
    if(pche2distLKrplane<30.)    {cleanup();return 0;}
    //
    dirs[cutcounter-1]->cd();

    //    ((TH1F*)gDirectory->FindObject(h_el1trktime->GetName()))        ->Fill( el1trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trkhodtime->GetName()))        ->Fill( el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trktimediff->GetName()))    ->Fill( el1trktime - el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime1->GetName()))            ->Fill( el1cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e1trkcltimediff->GetName()))    ->Fill( el1trktime- el1cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_el2trktime->GetName()))        ->Fill( el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trkhodtime->GetName()))        ->Fill( el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trktimediff->GetName()))    ->Fill( el2trktime - el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime2->GetName()))            ->Fill( el2cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e2trkcltimediff->GetName()))    ->Fill( el2trktime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_e1e2timediff->GetName()))        ->Fill( el1trktime-el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_e1e2cltimediff->GetName()))    ->Fill( el1cltime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_pchtrktime->GetName()))        ->Fill( pchtrktime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrkhodtime->GetName()))        ->Fill( pchtrkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_pchcltime->GetName()))            ->Fill( pchcltime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrktimediff->GetName()))    ->Fill( pchtrktime-pchtrkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrkcltimediff->GetName()))    ->Fill( pchtrktime-pchcltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill(elp1_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill(elp1_uncorr[1]);
    //    ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill(elv1_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill(elv1_uncorr[1]);
    //    ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill(elp1_uncorr[0], elp1_uncorr[1] );
    //    ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill(elv1_uncorr[0], elv1_uncorr[1] );
    //
    //    ((TH1F*)gDirectory->FindObject(h_bx_pch->GetName()))            ->Fill(pchp_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_by_pch->GetName()))            ->Fill(pchp_uncorr[1]);
    //    ((TH1F*)gDirectory->FindObject(h_bdxdz_pch->GetName()))            ->Fill(pchv_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_bdydz_pch->GetName()))            ->Fill(pchv_uncorr[1]) ;
    //    ((TH1F*)gDirectory->FindObject(h_bx_vs_by_pch->GetName()))        ->Fill(pchp_uncorr[0], pchp_uncorr[1]);
    //    ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_pch->GetName()))->Fill(pchv_uncorr[0], pchv_uncorr[1]);
    //
    //    ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( el1lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( el2lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el1lv->E() );
    //    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el2lv->E() );
    //    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el1lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el2lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el1lv->Pt() );
    //    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el2lv->Pt() );
    //
    //    ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
    ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);

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
    ((TH1F*)gDirectory->FindObject(h_pi0E->GetName()))->Fill(  pi0lv->E() );
    ((TH1F*)gDirectory->FindObject(h_pi0P->GetName()))->Fill(  pi0lv->P() );
    ((TH1F*)gDirectory->FindObject(h_pi0Pt->GetName()))->Fill( pi0lv->Pt() );
    ((TH1F*)gDirectory->FindObject(h_mofpi0->GetName()))->Fill( (*pi0lv).M() );

    ((TH1F*)gDirectory->FindObject(h_pi0Pcm->GetName()))->Fill( pi0lvcm->P() );
    ((TH1F*)gDirectory->FindObject(h_mofpi0diff->GetName()))->Fill( (*pi0lvcm).M()-(*pi0lv).M() );


    ((TH1F*)gDirectory->FindObject(h_gpieeE->GetName()))->Fill(  ((*glv)+(*pieelv)).E() );
    ((TH1F*)gDirectory->FindObject(h_gpieeP->GetName()))->Fill(  ((*glv)+(*pieelv)).P() );
    ((TH1F*)gDirectory->FindObject(h_gpieePt->GetName()))->Fill( ((*glv)+(*pieelv)).Pt());
    ((TH1F*)gDirectory->FindObject(h_gpieeM->GetName()))->Fill(  ((*glv)+(*pieelv)).M() );

    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"pi_ee_dist_LKr");
    cutcounter++;





    //------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------
    //                        charged pion and photon matching
    //------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------


    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~

    if(gpchdistLKrplane<30.)    {cleanup();return 0;}
    //
    dirs[cutcounter-1]->cd();

    //    ((TH1F*)gDirectory->FindObject(h_el1trktime->GetName()))        ->Fill( el1trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trkhodtime->GetName()))        ->Fill( el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trktimediff->GetName()))    ->Fill( el1trktime - el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime1->GetName()))            ->Fill( el1cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e1trkcltimediff->GetName()))    ->Fill( el1trktime- el1cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_el2trktime->GetName()))        ->Fill( el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trkhodtime->GetName()))        ->Fill( el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trktimediff->GetName()))    ->Fill( el2trktime - el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime2->GetName()))            ->Fill( el2cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e2trkcltimediff->GetName()))    ->Fill( el2trktime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_e1e2timediff->GetName()))        ->Fill( el1trktime-el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_e1e2cltimediff->GetName()))    ->Fill( el1cltime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_pchtrktime->GetName()))        ->Fill( pchtrktime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrkhodtime->GetName()))        ->Fill( pchtrkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_pchcltime->GetName()))            ->Fill( pchcltime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrktimediff->GetName()))    ->Fill( pchtrktime-pchtrkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrkcltimediff->GetName()))    ->Fill( pchtrktime-pchcltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill(elp1_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill(elp1_uncorr[1]);
    //    ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill(elv1_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill(elv1_uncorr[1]);
    //    ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill(elp1_uncorr[0], elp1_uncorr[1] );
    //    ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill(elv1_uncorr[0], elv1_uncorr[1] );
    //
    //    ((TH1F*)gDirectory->FindObject(h_bx_pch->GetName()))            ->Fill(pchp_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_by_pch->GetName()))            ->Fill(pchp_uncorr[1]);
    //    ((TH1F*)gDirectory->FindObject(h_bdxdz_pch->GetName()))            ->Fill(pchv_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_bdydz_pch->GetName()))            ->Fill(pchv_uncorr[1]) ;
    //    ((TH1F*)gDirectory->FindObject(h_bx_vs_by_pch->GetName()))        ->Fill(pchp_uncorr[0], pchp_uncorr[1]);
    //    ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_pch->GetName()))->Fill(pchv_uncorr[0], pchv_uncorr[1]);
    //
    //    ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( el1lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( el2lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el1lv->E() );
    //    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el2lv->E() );
    //    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el1lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el2lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el1lv->Pt() );
    //    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el2lv->Pt() );
    //
    //    ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
    ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);

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

    ((TH1F*)gDirectory->FindObject(h_pi0E->GetName()))->Fill(  pi0lv->E() );
    ((TH1F*)gDirectory->FindObject(h_pi0P->GetName()))->Fill(  pi0lv->P() );
    ((TH1F*)gDirectory->FindObject(h_pi0Pt->GetName()))->Fill( pi0lv->Pt() );
    ((TH1F*)gDirectory->FindObject(h_mofpi0->GetName()))->Fill( (*pi0lv).M() );

    ((TH1F*)gDirectory->FindObject(h_pi0Pcm->GetName()))->Fill( pi0lvcm->P() );
    ((TH1F*)gDirectory->FindObject(h_mofpi0diff->GetName()))->Fill( (*pi0lvcm).M()-(*pi0lv).M() );

    ((TH1F*)gDirectory->FindObject(h_gpieeE->GetName()))->Fill(  ((*glv)+(*pieelv)).E() );
    ((TH1F*)gDirectory->FindObject(h_gpieeP->GetName()))->Fill(  ((*glv)+(*pieelv)).P() );
    ((TH1F*)gDirectory->FindObject(h_gpieePt->GetName()))->Fill( ((*glv)+(*pieelv)).Pt());
    ((TH1F*)gDirectory->FindObject(h_gpieeM->GetName()))->Fill(  ((*glv)+(*pieelv)).M() );


    //
    Klv = new TLorentzVector( (*pchlv)+(*pi0lv) );
    //beta correction for Kaon
    //Klv->SetPxPyPzE( (1+beta)*Klv->Px(), (1+beta)*Klv->Py(), (1+beta)*Klv->Pz(), Klv->E() );
    TVector3 Kbv = Klv->BoostVector();
    ((TH1F*)gDirectory->FindObject(h_KE->GetName()))->Fill(Klv->E());
    ((TH1F*)gDirectory->FindObject(h_KP->GetName()))->Fill(Klv->P());
    ((TH1F*)gDirectory->FindObject(h_KPt->GetName()))->Fill(Klv->Pt());
    ((TH1F*)gDirectory->FindObject(h_mofK->GetName()))->Fill(Klv->M());

    pi0lvKcm = new TLorentzVector( *pi0lv );
    pchlvKcm = new TLorentzVector( *pchlv );
    pi0lvKcm->Boost(-Kbv);
    pchlvKcm->Boost(-Kbv);

    //
    KlvKcm = new TLorentzVector( (*pchlvKcm)+(*pi0lvKcm)  );
    ((TH1F*)gDirectory->FindObject(h_KPcm->GetName()))->Fill( KlvKcm->P() );
    //    ((TH1F*)gDirectory->FindObject(h_mofK->GetName()))->Fill( (*KlvKcm).M() );

    TVector3 pch3vKcm = pchlvKcm->Vect();
    TVector3 pi03vKcm = pi0lvKcm->Vect();
    //    ((TH1F*)gDirectory->FindObject(h_angle_pch_pi0_Kcm->GetName()))->Fill( pch3vKcm.Angle(pi03vKcm) );

    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"ph_pi_dist_LKr");
    cutcounter++;



    //    cout<<"20"<<endl;

    //------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------
    //                                        general cuts
    //------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------


    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~

    double mpi0 = pi0lv->M();
    if(mpi0<0.125 || mpi0>0.145)    return 0;
    //    if(mpi0<0.115 || mpi0>0.155)    {cleanup();return 0;}
    //
    dirs[cutcounter-1]->cd();

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

    ((TH1F*)gDirectory->FindObject(h_pchtrktime->GetName()))        ->Fill( pchtrktime );
    ((TH1F*)gDirectory->FindObject(h_pchtrkhodtime->GetName()))        ->Fill( pchtrkhodtime );
    ((TH1F*)gDirectory->FindObject(h_pchcltime->GetName()))            ->Fill( pchcltime );
    ((TH1F*)gDirectory->FindObject(h_pchtrktimediff->GetName()))    ->Fill( pchtrktime-pchtrkhodtime );
    ((TH1F*)gDirectory->FindObject(h_pchtrkcltimediff->GetName()))    ->Fill( pchtrktime-pchcltime );

    ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill(elp1_uncorr[0]);
    ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill(elp1_uncorr[1]);
    ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill(elv1_uncorr[0]);
    ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill(elv1_uncorr[1]);
    ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill(elp1_uncorr[0], elp1_uncorr[1] );
    ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill(elv1_uncorr[0], elv1_uncorr[1] );

    ((TH1F*)gDirectory->FindObject(h_bx_pch->GetName()))            ->Fill(pchp_uncorr[0]);
    ((TH1F*)gDirectory->FindObject(h_by_pch->GetName()))            ->Fill(pchp_uncorr[1]);
    ((TH1F*)gDirectory->FindObject(h_bdxdz_pch->GetName()))            ->Fill(pchv_uncorr[0]);
    ((TH1F*)gDirectory->FindObject(h_bdydz_pch->GetName()))            ->Fill(pchv_uncorr[1]) ;
    ((TH1F*)gDirectory->FindObject(h_bx_vs_by_pch->GetName()))        ->Fill(pchp_uncorr[0], pchp_uncorr[1]);
    ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_pch->GetName()))->Fill(pchv_uncorr[0], pchv_uncorr[1]);

    ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( el1lv->P() );
    ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( el2lv->P() );
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

    ((TH1F*)gDirectory->FindObject(h_gE->GetName()))->Fill( (*glv).E() );
    ((TH1F*)gDirectory->FindObject(h_gP->GetName()))->Fill( (*glv).P() );
    ((TH1F*)gDirectory->FindObject(h_gPt->GetName()))->Fill( (*glv).Pt() );
    ((TH1F*)gDirectory->FindObject(h_gtrkRadDHC1->GetName()))->Fill( gtrkRadDHC1 );

    ((TH1F*)gDirectory->FindObject(h_pchlvE->GetName()))->Fill( (*pchlv).E() );
    ((TH1F*)gDirectory->FindObject(h_pchlvP->GetName()))->Fill( (*pchlv).P() );
    ((TH1F*)gDirectory->FindObject(h_pchlvPt->GetName()))->Fill( (*pchlv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_angle_e1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_angle_ph_e1e2->GetName()))->Fill( (*eelv).Vect().Angle(glv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_angle_pch_pi0->GetName()))->Fill( (*pi0lv).Vect().Angle(pchlv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_a_e1e2_phe1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()), (*eelv).Vect().Angle(glv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_a_e1e2_pchpi0->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()), (*pi0lv).Vect().Angle(pchlv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_a_phe1e2_pchpi0->GetName()))->Fill( (*eelv).Vect().Angle(glv->Vect()), (*pi0lv).Vect().Angle(pchlv->Vect()) );

    ////    ((TH1F*)gDirectory->FindObject(h_angle_pch_pi0_Kcm->GetName()))->Fill( pch3vKcm.Angle(pi03vKcm) );

    ((TH1F*)gDirectory->FindObject(h_cda_ee->GetName()))->Fill(cda);
    ((TH1F*)gDirectory->FindObject(h_vx->GetName()))->Fill(vertex[0]);
    ((TH1F*)gDirectory->FindObject(h_vy->GetName()))->Fill(vertex[1]);
    ((TH1F*)gDirectory->FindObject(h_vz->GetName()))->Fill(vertex[2]);
    ((TH2F*)gDirectory->FindObject(h_vx_vs_vy->GetName()))->Fill(vertex[0],vertex[1]);

    ((TH1F*)gDirectory->FindObject(h_cda_e1pch->GetName()))->Fill(cda_e1pch);
    ((TH1F*)gDirectory->FindObject(h_vx_e1pch->GetName()))->Fill(vertex_e1pch[0]);
    ((TH1F*)gDirectory->FindObject(h_vy_e1pch->GetName()))->Fill(vertex_e1pch[1]);
    ((TH1F*)gDirectory->FindObject(h_vz_e1pch->GetName()))->Fill(vertex_e1pch[2]);
    ((TH2F*)gDirectory->FindObject(h_vx_vs_vy_e1pch->GetName()))->Fill(vertex_e1pch[0],vertex_e1pch[1]);

    ((TH1F*)gDirectory->FindObject(h_vz_e1pch_diff->GetName()))->Fill( vertex[2] - vertex_e1pch[2] );

    ((TH1F*)gDirectory->FindObject(h_pi0E->GetName()))->Fill(  pi0lv->E() );
    ((TH1F*)gDirectory->FindObject(h_pi0P->GetName()))->Fill(  pi0lv->P() );
    ((TH1F*)gDirectory->FindObject(h_pi0Pt->GetName()))->Fill( pi0lv->Pt() );
    ((TH1F*)gDirectory->FindObject(h_mofpi0->GetName()))->Fill( (*pi0lv).M() );

    ((TH1F*)gDirectory->FindObject(h_pi0Pcm->GetName()))->Fill( pi0lvcm->P() );
    ((TH1F*)gDirectory->FindObject(h_mofpi0diff->GetName()))->Fill( (*pi0lvcm).M()-(*pi0lv).M() );

    //    ((TH1F*)gDirectory->FindObject(h_gpieeE->GetName()))->Fill(  ((*glv)+(*pieelv)).E() );
    //    ((TH1F*)gDirectory->FindObject(h_gpieeP->GetName()))->Fill(  ((*glv)+(*pieelv)).P() );
    //    ((TH1F*)gDirectory->FindObject(h_gpieePt->GetName()))->Fill( ((*glv)+(*pieelv)).Pt());
    //    ((TH1F*)gDirectory->FindObject(h_gpieeM->GetName()))->Fill(  ((*glv)+(*pieelv)).M() );

    ((TH1F*)gDirectory->FindObject(h_KE->GetName()))->Fill(Klv->E());
    ((TH1F*)gDirectory->FindObject(h_KP->GetName()))->Fill(Klv->P());
    ((TH1F*)gDirectory->FindObject(h_KPcm->GetName()))->Fill( KlvKcm->P() );
    ((TH1F*)gDirectory->FindObject(h_KPt->GetName()))->Fill(Klv->Pt());
    ((TH1F*)gDirectory->FindObject(h_mofK->GetName()))->Fill(Klv->M());

    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"pi0_mass_wndw");
    cutcounter++;



    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~

    double pK = Klv->P();
    if(pK<54. || pK>66.)    {cleanup();return 0;}
    //
    dirs[cutcounter-1]->cd();

    //    ((TH1F*)gDirectory->FindObject(h_el1trktime->GetName()))        ->Fill( el1trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trkhodtime->GetName()))        ->Fill( el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trktimediff->GetName()))    ->Fill( el1trktime - el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime1->GetName()))            ->Fill( el1cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e1trkcltimediff->GetName()))    ->Fill( el1trktime- el1cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_el2trktime->GetName()))        ->Fill( el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trkhodtime->GetName()))        ->Fill( el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trktimediff->GetName()))    ->Fill( el2trktime - el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime2->GetName()))            ->Fill( el2cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e2trkcltimediff->GetName()))    ->Fill( el2trktime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_e1e2timediff->GetName()))        ->Fill( el1trktime-el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_e1e2cltimediff->GetName()))    ->Fill( el1cltime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_pchtrktime->GetName()))        ->Fill( pchtrktime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrkhodtime->GetName()))        ->Fill( pchtrkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_pchcltime->GetName()))            ->Fill( pchcltime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrktimediff->GetName()))    ->Fill( pchtrktime-pchtrkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrkcltimediff->GetName()))    ->Fill( pchtrktime-pchcltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill(elp1_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill(elp1_uncorr[1]);
    //    ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill(elv1_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill(elv1_uncorr[1]);
    //    ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill(elp1_uncorr[0], elp1_uncorr[1] );
    //    ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill(elv1_uncorr[0], elv1_uncorr[1] );
    //
    //    ((TH1F*)gDirectory->FindObject(h_bx_pch->GetName()))            ->Fill(pchp_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_by_pch->GetName()))            ->Fill(pchp_uncorr[1]);
    //    ((TH1F*)gDirectory->FindObject(h_bdxdz_pch->GetName()))            ->Fill(pchv_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_bdydz_pch->GetName()))            ->Fill(pchv_uncorr[1]) ;
    //    ((TH1F*)gDirectory->FindObject(h_bx_vs_by_pch->GetName()))        ->Fill(pchp_uncorr[0], pchp_uncorr[1]);
    //    ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_pch->GetName()))->Fill(pchv_uncorr[0], pchv_uncorr[1]);
    //
    //    ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( el1lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( el2lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el1lv->E() );
    //    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el2lv->E() );
    //    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el1lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el2lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el1lv->Pt() );
    //    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el2lv->Pt() );

    //    ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
    ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);

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

    //    ((TH1F*)gDirectory->FindObject(h_angle_pch_pi0_Kcm->GetName()))->Fill( pch3vKcm.Angle(pi03vKcm) );

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

    ((TH1F*)gDirectory->FindObject(h_pi0E->GetName()))->Fill(  pi0lv->E() );
    ((TH1F*)gDirectory->FindObject(h_pi0P->GetName()))->Fill(  pi0lv->P() );
    ((TH1F*)gDirectory->FindObject(h_pi0Pt->GetName()))->Fill( pi0lv->Pt() );
    ((TH1F*)gDirectory->FindObject(h_mofpi0->GetName()))->Fill( (*pi0lv).M() );

    ((TH1F*)gDirectory->FindObject(h_pi0Pcm->GetName()))->Fill( pi0lvcm->P() );
    ((TH1F*)gDirectory->FindObject(h_mofpi0diff->GetName()))->Fill( (*pi0lvcm).M()-(*pi0lv).M() );

    //    ((TH1F*)gDirectory->FindObject(h_gpieeE->GetName()))->Fill(  ((*glv)+(*pieelv)).E() );
    //    ((TH1F*)gDirectory->FindObject(h_gpieeP->GetName()))->Fill(  ((*glv)+(*pieelv)).P() );
    //    ((TH1F*)gDirectory->FindObject(h_gpieePt->GetName()))->Fill( ((*glv)+(*pieelv)).Pt());
    //    ((TH1F*)gDirectory->FindObject(h_gpieeM->GetName()))->Fill(  ((*glv)+(*pieelv)).M() );

    ((TH1F*)gDirectory->FindObject(h_KE->GetName()))->Fill(Klv->E());
    ((TH1F*)gDirectory->FindObject(h_KP->GetName()))->Fill(Klv->P());
    ((TH1F*)gDirectory->FindObject(h_KPcm->GetName()))->Fill( KlvKcm->P() );
    ((TH1F*)gDirectory->FindObject(h_KPt->GetName()))->Fill(Klv->Pt());
    ((TH1F*)gDirectory->FindObject(h_mofK->GetName()))->Fill(Klv->M());

    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"K_mom_wndw");
    cutcounter++;



    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~

    double ptK = Klv->Pt();
    if(ptK>0.02) {cleanup();return 0;}
    //
    dirs[cutcounter-1]->cd();

    //    ((TH1F*)gDirectory->FindObject(h_el1trktime->GetName()))        ->Fill( el1trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trkhodtime->GetName()))        ->Fill( el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el1trktimediff->GetName()))    ->Fill( el1trktime - el1trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime1->GetName()))            ->Fill( el1cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e1trkcltimediff->GetName()))    ->Fill( el1trktime- el1cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_el2trktime->GetName()))        ->Fill( el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trkhodtime->GetName()))        ->Fill( el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_el2trktimediff->GetName()))    ->Fill( el2trktime - el2trkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_ecltime2->GetName()))            ->Fill( el2cltime );
    //    ((TH1F*)gDirectory->FindObject(h_e2trkcltimediff->GetName()))    ->Fill( el2trktime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_e1e2timediff->GetName()))        ->Fill( el1trktime-el2trktime );
    //    ((TH1F*)gDirectory->FindObject(h_e1e2cltimediff->GetName()))    ->Fill( el1cltime-el2cltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_pchtrktime->GetName()))        ->Fill( pchtrktime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrkhodtime->GetName()))        ->Fill( pchtrkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_pchcltime->GetName()))            ->Fill( pchcltime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrktimediff->GetName()))    ->Fill( pchtrktime-pchtrkhodtime );
    //    ((TH1F*)gDirectory->FindObject(h_pchtrkcltimediff->GetName()))    ->Fill( pchtrktime-pchcltime );
    //
    //    ((TH1F*)gDirectory->FindObject(h_bx->GetName()))                ->Fill(elp1_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_by->GetName()))                ->Fill(elp1_uncorr[1]);
    //    ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))                ->Fill(elv1_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))                ->Fill(elv1_uncorr[1]);
    //    ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))            ->Fill(elp1_uncorr[0], elp1_uncorr[1] );
    //    ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))    ->Fill(elv1_uncorr[0], elv1_uncorr[1] );
    //
    //    ((TH1F*)gDirectory->FindObject(h_bx_pch->GetName()))            ->Fill(pchp_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_by_pch->GetName()))            ->Fill(pchp_uncorr[1]);
    //    ((TH1F*)gDirectory->FindObject(h_bdxdz_pch->GetName()))            ->Fill(pchv_uncorr[0]);
    //    ((TH1F*)gDirectory->FindObject(h_bdydz_pch->GetName()))            ->Fill(pchv_uncorr[1]) ;
    //    ((TH1F*)gDirectory->FindObject(h_bx_vs_by_pch->GetName()))        ->Fill(pchp_uncorr[0], pchp_uncorr[1]);
    //    ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_pch->GetName()))->Fill(pchv_uncorr[0], pchv_uncorr[1]);
    //
    //    ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( el1lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( el2lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el1lv->E() );
    //    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el2lv->E() );
    //    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el1lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el2lv->P() );
    //    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el1lv->Pt() );
    //    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el2lv->Pt() );

    //    ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
    ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);

    //    ((TH1F*)gDirectory->FindObject(h_gE->GetName()))->Fill( (*glv).E() );
    //    ((TH1F*)gDirectory->FindObject(h_gP->GetName()))->Fill( (*glv).P() );
    //    ((TH1F*)gDirectory->FindObject(h_gPt->GetName()))->Fill( (*glv).Pt() );
    //    ((TH1F*)gDirectory->FindObject(h_gtrkRadDHC1->GetName()))->Fill( gtrkRadDHC1 );
    //
    //    ((TH1F*)gDirectory->FindObject(h_pchlvE->GetName()))->Fill( (*pchlv).E() );
    //    ((TH1F*)gDirectory->FindObject(h_pchlvP->GetName()))->Fill( (*pchlv).P() );
    //    ((TH1F*)gDirectory->FindObject(h_pchlvPt->GetName()))->Fill( (*pchlv).Pt() );


    //    ((TH1F*)gDirectory->FindObject(h_angle_e1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_angle_ph_e1e2->GetName()))->Fill( (*eelv).Vect().Angle(glv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_angle_pch_pi0->GetName()))->Fill( (*pi0lv).Vect().Angle(pchlv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_a_e1e2_phe1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()), (*eelv).Vect().Angle(glv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_a_e1e2_pchpi0->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()), (*pi0lv).Vect().Angle(pchlv->Vect()) );
    //    ((TH1F*)gDirectory->FindObject(h_a_phe1e2_pchpi0->GetName()))->Fill( (*eelv).Vect().Angle(glv->Vect()), (*pi0lv).Vect().Angle(pchlv->Vect()) );

    //    ((TH1F*)gDirectory->FindObject(h_angle_pch_pi0_Kcm->GetName()))->Fill( pch3vKcm.Angle(pi03vKcm) );

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

    ((TH1F*)gDirectory->FindObject(h_pi0E->GetName()))->Fill(  pi0lv->E() );
    ((TH1F*)gDirectory->FindObject(h_pi0P->GetName()))->Fill(  pi0lv->P() );
    ((TH1F*)gDirectory->FindObject(h_pi0Pt->GetName()))->Fill( pi0lv->Pt() );
    ((TH1F*)gDirectory->FindObject(h_mofpi0->GetName()))->Fill( (*pi0lv).M() );

    ((TH1F*)gDirectory->FindObject(h_pi0Pcm->GetName()))->Fill( pi0lvcm->P() );
    ((TH1F*)gDirectory->FindObject(h_mofpi0diff->GetName()))->Fill( (*pi0lvcm).M()-(*pi0lv).M() );

    ((TH1F*)gDirectory->FindObject(h_gpieeE->GetName()))->Fill(  ((*glv)+(*pieelv)).E() );
    ((TH1F*)gDirectory->FindObject(h_gpieeP->GetName()))->Fill(  ((*glv)+(*pieelv)).P() );
    ((TH1F*)gDirectory->FindObject(h_gpieePt->GetName()))->Fill( ((*glv)+(*pieelv)).Pt());
    ((TH1F*)gDirectory->FindObject(h_gpieeM->GetName()))->Fill(  ((*glv)+(*pieelv)).M() );

    ((TH1F*)gDirectory->FindObject(h_KE->GetName()))->Fill(Klv->E());
    ((TH1F*)gDirectory->FindObject(h_KP->GetName()))->Fill(Klv->P());
    ((TH1F*)gDirectory->FindObject(h_KPcm->GetName()))->Fill( KlvKcm->P() );
    ((TH1F*)gDirectory->FindObject(h_KPt->GetName()))->Fill(Klv->Pt());
    ((TH1F*)gDirectory->FindObject(h_mofK->GetName()))->Fill(Klv->M());

    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"K_tmom_wndw");
    cutcounter++;




    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~

    double mK = Klv->M();
    if(mK<0.475 || mK>0.510) {cleanup(); return 0;}
    //
    dirs[cutcounter-1]->cd();

    ((TH1I*)gDirectory->FindObject(h_evttypebytrks->GetName()))->Fill(evttypebytrks);

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

    ((TH1F*)gDirectory->FindObject(h_pchtrktime->GetName()))        ->Fill( pchtrktime );
    ((TH1F*)gDirectory->FindObject(h_pchtrkhodtime->GetName()))        ->Fill( pchtrkhodtime );
    ((TH1F*)gDirectory->FindObject(h_pchcltime->GetName()))            ->Fill( pchcltime );
    ((TH1F*)gDirectory->FindObject(h_pchtrktimediff->GetName()))    ->Fill( pchtrktime-pchtrkhodtime );
    ((TH1F*)gDirectory->FindObject(h_pchtrkcltimediff->GetName()))    ->Fill( pchtrktime-pchcltime );

    for (i=0; i<ntrack; i++) {
        p_track        = sevt->track[i].p;
        icl_track    = sevt->track[i].iClus;
        e_track      = sevt->cluster[icl_track].energy;
        eop_track      = e_track / p_track;
        ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_track);
        ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_track, p_track);
    }

    ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( el1lv->P() );
    ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( el2lv->P() );
    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el1lv->E() );
    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el2lv->E() );
    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el1lv->P() );
    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el2lv->P() );
    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el1lv->Pt() );
    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el2lv->Pt() );

    ((TH1F*)gDirectory->FindObject(h_gE->GetName()))->Fill( (*glv).E() );
    ((TH1F*)gDirectory->FindObject(h_gP->GetName()))->Fill( (*glv).P() );
    ((TH1F*)gDirectory->FindObject(h_gPt->GetName()))->Fill( (*glv).Pt() );
    ((TH1F*)gDirectory->FindObject(h_gtrkRadDHC1->GetName()))->Fill( gtrkRadDHC1 );

    ((TH1F*)gDirectory->FindObject(h_pchlvE->GetName()))->Fill( (*pchlv).E() );
    ((TH1F*)gDirectory->FindObject(h_pchlvP->GetName()))->Fill( (*pchlv).P() );
    ((TH1F*)gDirectory->FindObject(h_pchlvPt->GetName()))->Fill( (*pchlv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
    ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_pi0E->GetName()))->Fill(  pi0lv->E() );
    ((TH1F*)gDirectory->FindObject(h_pi0P->GetName()))->Fill(  pi0lv->P() );
    ((TH1F*)gDirectory->FindObject(h_pi0Pt->GetName()))->Fill( pi0lv->Pt() );
    ((TH1F*)gDirectory->FindObject(h_pi0Pcm->GetName()))->Fill( pi0lvcm->P() );
    ((TH1F*)gDirectory->FindObject(h_mofpi0->GetName()))->Fill( (*pi0lv).M() );
    ((TH1F*)gDirectory->FindObject(h_mofpi0diff->GetName()))->Fill( (*pi0lvcm).M()-(*pi0lv).M() );

    ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);
    ((TH2F*)gDirectory->FindObject(h_pi0_dalitz->GetName()))->Fill(msq12, msq23);

    ((TH1F*)gDirectory->FindObject(h_pieeE->GetName()))->Fill(pieelv->E());
    ((TH1F*)gDirectory->FindObject(h_pieeP->GetName()))->Fill(pieelv->P());
    ((TH1F*)gDirectory->FindObject(h_pieePt->GetName()))->Fill(pieelv->Pt());
    ((TH1F*)gDirectory->FindObject(h_pieeM->GetName()))->Fill(pieelv->M());

    ((TH1F*)gDirectory->FindObject(h_gpieeE->GetName()))->Fill(  ((*glv)+(*pieelv)).E() );
    ((TH1F*)gDirectory->FindObject(h_gpieeP->GetName()))->Fill(  ((*glv)+(*pieelv)).P() );
    ((TH1F*)gDirectory->FindObject(h_gpieePt->GetName()))->Fill( ((*glv)+(*pieelv)).Pt());
    ((TH1F*)gDirectory->FindObject(h_gpieeM->GetName()))->Fill(  ((*glv)+(*pieelv)).M() );

    ((TH1F*)gDirectory->FindObject(h_KE->GetName()))->Fill(Klv->E());
    ((TH1F*)gDirectory->FindObject(h_KP->GetName()))->Fill(Klv->P());
    ((TH1F*)gDirectory->FindObject(h_KPcm->GetName()))->Fill( KlvKcm->P() );
    ((TH1F*)gDirectory->FindObject(h_KPt->GetName()))->Fill(Klv->Pt());
    ((TH1F*)gDirectory->FindObject(h_mofK->GetName()))->Fill(Klv->M());


    ((TH1F*)gDirectory->FindObject(h_eldistDCH1plane->GetName()))->Fill( eldistDCH1plane );
    //    ((TH2F*)gDirectory->FindObject(h_eldistDCH1_ee_vz->GetName()))->Fill( eldistDCH1plane, vertex[2] );
    //    ((TH2F*)gDirectory->FindObject(h_eldistDCH1_pch_vz->GetName()))->Fill( eldistDCH1plane, vertex_e1pch[2] );
    ((TH1F*)gDirectory->FindObject(h_eldistLKrplane->GetName()))->Fill( eldistLKrplane );
    ((TH1F*)gDirectory->FindObject(h_gel1distLKrplane->GetName()))->Fill( gel1distLKrplane );
    ((TH1F*)gDirectory->FindObject(h_gel2distLKrplane->GetName()))->Fill( gel2distLKrplane );
    ((TH1F*)gDirectory->FindObject(h_pche1distDHC1plane->GetName()))->Fill( pche1distDHC1plane );
    ((TH1F*)gDirectory->FindObject(h_pche2distDHC1plane->GetName()))->Fill( pche2distDHC1plane );
    ((TH1F*)gDirectory->FindObject(h_pche1distLKrplane->GetName()))->Fill( pche1distLKrplane );
    ((TH1F*)gDirectory->FindObject(h_pche2distLKrplane->GetName()))->Fill( pche2distLKrplane );
    ((TH1F*)gDirectory->FindObject(h_gpchdistLKrplane->GetName()))->Fill( gpchdistLKrplane );

    ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill(elp1_uncorr[0]);
    ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill(elp1_uncorr[1]);
    ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill(elv1_uncorr[0]);
    ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill(elv1_uncorr[1]);
    ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill(elp1_uncorr[0], elp1_uncorr[1] );
    ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill(elv1_uncorr[0], elv1_uncorr[1] );

    ((TH1F*)gDirectory->FindObject(h_bx_pch->GetName()))            ->Fill(pchp_uncorr[0]);
    ((TH1F*)gDirectory->FindObject(h_by_pch->GetName()))            ->Fill(pchp_uncorr[1]);
    ((TH1F*)gDirectory->FindObject(h_bdxdz_pch->GetName()))            ->Fill(pchv_uncorr[0]);
    ((TH1F*)gDirectory->FindObject(h_bdydz_pch->GetName()))            ->Fill(pchv_uncorr[1]) ;
    ((TH1F*)gDirectory->FindObject(h_bx_vs_by_pch->GetName()))        ->Fill(pchp_uncorr[0], pchp_uncorr[1]);
    ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_pch->GetName()))->Fill(pchv_uncorr[0], pchv_uncorr[1]);

    ((TH1F*)gDirectory->FindObject(h_cda_ee->GetName()))->Fill(cda);
    ((TH1F*)gDirectory->FindObject(h_vx->GetName()))->Fill(vertex[0]);
    ((TH1F*)gDirectory->FindObject(h_vy->GetName()))->Fill(vertex[1]);
    ((TH1F*)gDirectory->FindObject(h_vz->GetName()))->Fill(vertex[2]);
    ((TH2F*)gDirectory->FindObject(h_vx_vs_vy->GetName()))->Fill(vertex[0],vertex[1]);

    ((TH1F*)gDirectory->FindObject(h_cda_e1pch->GetName()))->Fill(cda_e1pch);
    ((TH1F*)gDirectory->FindObject(h_vx_e1pch->GetName()))->Fill(vertex_e1pch[0]);
    ((TH1F*)gDirectory->FindObject(h_vy_e1pch->GetName()))->Fill(vertex_e1pch[1]);
    ((TH1F*)gDirectory->FindObject(h_vz_e1pch->GetName()))->Fill(vertex_e1pch[2]);
    ((TH2F*)gDirectory->FindObject(h_vx_vs_vy_e1pch->GetName()))->Fill(vertex_e1pch[0],vertex_e1pch[1]);

    ((TH1F*)gDirectory->FindObject(h_vz_e1pch_diff->GetName()))->Fill( vertex[2] - vertex_e1pch[2] );

    ((TH1F*)gDirectory->FindObject(h_angle_e1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_angle_ph_e1e2->GetName()))->Fill( (*eelv).Vect().Angle(glv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_angle_pch_pi0->GetName()))->Fill( (*pi0lv).Vect().Angle(pchlv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_a_e1e2_phe1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()), (*eelv).Vect().Angle(glv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_a_e1e2_pchpi0->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()), (*pi0lv).Vect().Angle(pchlv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_a_phe1e2_pchpi0->GetName()))->Fill( (*eelv).Vect().Angle(glv->Vect()), (*pi0lv).Vect().Angle(pchlv->Vect()) );

    //    ((TH1F*)gDirectory->FindObject(h_angle_pch_pi0_Kcm->GetName()))->Fill( pch3vKcm.Angle(pi03vKcm) );


    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"K_mass_wndw");
    cutcounter++;




    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~

    if((*pchlv).P()<10.) {cleanup(); return 0;}
    //
    dirs[cutcounter-1]->cd();

    ((TH1I*)gDirectory->FindObject(h_evttypebytrks->GetName()))->Fill(evttypebytrks);

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

    ((TH1F*)gDirectory->FindObject(h_pchtrktime->GetName()))        ->Fill( pchtrktime );
    ((TH1F*)gDirectory->FindObject(h_pchtrkhodtime->GetName()))        ->Fill( pchtrkhodtime );
    ((TH1F*)gDirectory->FindObject(h_pchcltime->GetName()))            ->Fill( pchcltime );
    ((TH1F*)gDirectory->FindObject(h_pchtrktimediff->GetName()))    ->Fill( pchtrktime-pchtrkhodtime );
    ((TH1F*)gDirectory->FindObject(h_pchtrkcltimediff->GetName()))    ->Fill( pchtrktime-pchcltime );

    for (i=0; i<ntrack; i++) {
        p_track        = sevt->track[i].p;
        icl_track    = sevt->track[i].iClus;
        e_track      = sevt->cluster[icl_track].energy;
        eop_track      = e_track / p_track;
        ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_track);
        ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_track, p_track);
    }

    ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( el1lv->P() );
    ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( el2lv->P() );
    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el1lv->E() );
    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el2lv->E() );
    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el1lv->P() );
    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el2lv->P() );
    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el1lv->Pt() );
    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el2lv->Pt() );

    ((TH1F*)gDirectory->FindObject(h_gE->GetName()))->Fill( (*glv).E() );
    ((TH1F*)gDirectory->FindObject(h_gP->GetName()))->Fill( (*glv).P() );
    ((TH1F*)gDirectory->FindObject(h_gPt->GetName()))->Fill( (*glv).Pt() );
    ((TH1F*)gDirectory->FindObject(h_gtrkRadDHC1->GetName()))->Fill( gtrkRadDHC1 );

    ((TH1F*)gDirectory->FindObject(h_pchlvE->GetName()))->Fill( (*pchlv).E() );
    ((TH1F*)gDirectory->FindObject(h_pchlvP->GetName()))->Fill( (*pchlv).P() );
    ((TH1F*)gDirectory->FindObject(h_pchlvPt->GetName()))->Fill( (*pchlv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
    ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_pi0E->GetName()))->Fill(  pi0lv->E() );
    ((TH1F*)gDirectory->FindObject(h_pi0P->GetName()))->Fill(  pi0lv->P() );
    ((TH1F*)gDirectory->FindObject(h_pi0Pt->GetName()))->Fill( pi0lv->Pt() );
    ((TH1F*)gDirectory->FindObject(h_pi0Pcm->GetName()))->Fill( pi0lvcm->P() );
    ((TH1F*)gDirectory->FindObject(h_mofpi0->GetName()))->Fill( (*pi0lv).M() );
    ((TH1F*)gDirectory->FindObject(h_mofpi0diff->GetName()))->Fill( (*pi0lvcm).M()-(*pi0lv).M() );

    ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);
    ((TH2F*)gDirectory->FindObject(h_pi0_dalitz->GetName()))->Fill(msq12, msq23);

    ((TH1F*)gDirectory->FindObject(h_pieeE->GetName()))->Fill(pieelv->E());
    ((TH1F*)gDirectory->FindObject(h_pieeP->GetName()))->Fill(pieelv->P());
    ((TH1F*)gDirectory->FindObject(h_pieePt->GetName()))->Fill(pieelv->Pt());
    ((TH1F*)gDirectory->FindObject(h_pieeM->GetName()))->Fill(pieelv->M());

    ((TH1F*)gDirectory->FindObject(h_gpieeE->GetName()))->Fill(  ((*glv)+(*pieelv)).E() );
    ((TH1F*)gDirectory->FindObject(h_gpieeP->GetName()))->Fill(  ((*glv)+(*pieelv)).P() );
    ((TH1F*)gDirectory->FindObject(h_gpieePt->GetName()))->Fill( ((*glv)+(*pieelv)).Pt());
    ((TH1F*)gDirectory->FindObject(h_gpieeM->GetName()))->Fill(  ((*glv)+(*pieelv)).M() );

    ((TH1F*)gDirectory->FindObject(h_KE->GetName()))->Fill(Klv->E());
    ((TH1F*)gDirectory->FindObject(h_KP->GetName()))->Fill(Klv->P());
    ((TH1F*)gDirectory->FindObject(h_KPcm->GetName()))->Fill( KlvKcm->P() );
    ((TH1F*)gDirectory->FindObject(h_KPt->GetName()))->Fill(Klv->Pt());
    ((TH1F*)gDirectory->FindObject(h_mofK->GetName()))->Fill(Klv->M());


    ((TH1F*)gDirectory->FindObject(h_eldistDCH1plane->GetName()))->Fill( eldistDCH1plane );
    //    ((TH2F*)gDirectory->FindObject(h_eldistDCH1_ee_vz->GetName()))->Fill( eldistDCH1plane, vertex[2] );
    //    ((TH2F*)gDirectory->FindObject(h_eldistDCH1_pch_vz->GetName()))->Fill( eldistDCH1plane, vertex_e1pch[2] );
    ((TH1F*)gDirectory->FindObject(h_eldistLKrplane->GetName()))->Fill( eldistLKrplane );
    ((TH1F*)gDirectory->FindObject(h_gel1distLKrplane->GetName()))->Fill( gel1distLKrplane );
    ((TH1F*)gDirectory->FindObject(h_gel2distLKrplane->GetName()))->Fill( gel2distLKrplane );
    ((TH1F*)gDirectory->FindObject(h_pche1distDHC1plane->GetName()))->Fill( pche1distDHC1plane );
    ((TH1F*)gDirectory->FindObject(h_pche2distDHC1plane->GetName()))->Fill( pche2distDHC1plane );
    ((TH1F*)gDirectory->FindObject(h_pche1distLKrplane->GetName()))->Fill( pche1distLKrplane );
    ((TH1F*)gDirectory->FindObject(h_pche2distLKrplane->GetName()))->Fill( pche2distLKrplane );
    ((TH1F*)gDirectory->FindObject(h_gpchdistLKrplane->GetName()))->Fill( gpchdistLKrplane );

    ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill(elp1_uncorr[0]);
    ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill(elp1_uncorr[1]);
    ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill(elv1_uncorr[0]);
    ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill(elv1_uncorr[1]);
    ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill(elp1_uncorr[0], elp1_uncorr[1] );
    ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill(elv1_uncorr[0], elv1_uncorr[1] );

    ((TH1F*)gDirectory->FindObject(h_bx_pch->GetName()))            ->Fill(pchp_uncorr[0]);
    ((TH1F*)gDirectory->FindObject(h_by_pch->GetName()))            ->Fill(pchp_uncorr[1]);
    ((TH1F*)gDirectory->FindObject(h_bdxdz_pch->GetName()))            ->Fill(pchv_uncorr[0]);
    ((TH1F*)gDirectory->FindObject(h_bdydz_pch->GetName()))            ->Fill(pchv_uncorr[1]) ;
    ((TH1F*)gDirectory->FindObject(h_bx_vs_by_pch->GetName()))        ->Fill(pchp_uncorr[0], pchp_uncorr[1]);
    ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_pch->GetName()))->Fill(pchv_uncorr[0], pchv_uncorr[1]);

    ((TH1F*)gDirectory->FindObject(h_cda_ee->GetName()))->Fill(cda);
    ((TH1F*)gDirectory->FindObject(h_vx->GetName()))->Fill(vertex[0]);
    ((TH1F*)gDirectory->FindObject(h_vy->GetName()))->Fill(vertex[1]);
    ((TH1F*)gDirectory->FindObject(h_vz->GetName()))->Fill(vertex[2]);
    ((TH2F*)gDirectory->FindObject(h_vx_vs_vy->GetName()))->Fill(vertex[0],vertex[1]);

    ((TH1F*)gDirectory->FindObject(h_cda_e1pch->GetName()))->Fill(cda_e1pch);
    ((TH1F*)gDirectory->FindObject(h_vx_e1pch->GetName()))->Fill(vertex_e1pch[0]);
    ((TH1F*)gDirectory->FindObject(h_vy_e1pch->GetName()))->Fill(vertex_e1pch[1]);
    ((TH1F*)gDirectory->FindObject(h_vz_e1pch->GetName()))->Fill(vertex_e1pch[2]);
    ((TH2F*)gDirectory->FindObject(h_vx_vs_vy_e1pch->GetName()))->Fill(vertex_e1pch[0],vertex_e1pch[1]);

    ((TH1F*)gDirectory->FindObject(h_vz_e1pch_diff->GetName()))->Fill( vertex[2] - vertex_e1pch[2] );

    ((TH1F*)gDirectory->FindObject(h_angle_e1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_angle_ph_e1e2->GetName()))->Fill( (*eelv).Vect().Angle(glv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_angle_pch_pi0->GetName()))->Fill( (*pi0lv).Vect().Angle(pchlv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_a_e1e2_phe1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()), (*eelv).Vect().Angle(glv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_a_e1e2_pchpi0->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()), (*pi0lv).Vect().Angle(pchlv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_a_phe1e2_pchpi0->GetName()))->Fill( (*eelv).Vect().Angle(glv->Vect()), (*pi0lv).Vect().Angle(pchlv->Vect()) );

    //    ((TH1F*)gDirectory->FindObject(h_angle_pch_pi0_Kcm->GetName()))->Fill( pch3vKcm.Angle(pi03vKcm) );


    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"pi_P");
    cutcounter++;




    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~

    p_track        = sevt->track[pchtrkindx].p;
    icl_track    = sevt->track[pchtrkindx].iClus;
    e_track      = sevt->cluster[icl_track].energy;
    eop_track      = e_track / p_track;

    if(eop_track<0.1) {cleanup(); return 0;}
    //
    dirs[cutcounter-1]->cd();

    ((TH1I*)gDirectory->FindObject(h_evttypebytrks->GetName()))->Fill(evttypebytrks);

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

    ((TH1F*)gDirectory->FindObject(h_pchtrktime->GetName()))        ->Fill( pchtrktime );
    ((TH1F*)gDirectory->FindObject(h_pchtrkhodtime->GetName()))        ->Fill( pchtrkhodtime );
    ((TH1F*)gDirectory->FindObject(h_pchcltime->GetName()))            ->Fill( pchcltime );
    ((TH1F*)gDirectory->FindObject(h_pchtrktimediff->GetName()))    ->Fill( pchtrktime-pchtrkhodtime );
    ((TH1F*)gDirectory->FindObject(h_pchtrkcltimediff->GetName()))    ->Fill( pchtrktime-pchcltime );

    for (i=0; i<ntrack; i++) {
        p_track        = sevt->track[i].p;
        icl_track    = sevt->track[i].iClus;
        e_track      = sevt->cluster[icl_track].energy;
        eop_track      = e_track / p_track;
        ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_track);
        ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_track, p_track);
    }

    ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( el1lv->P() );
    ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( el2lv->P() );
    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el1lv->E() );
    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el2lv->E() );
    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el1lv->P() );
    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el2lv->P() );
    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el1lv->Pt() );
    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el2lv->Pt() );

    ((TH1F*)gDirectory->FindObject(h_gE->GetName()))->Fill( (*glv).E() );
    ((TH1F*)gDirectory->FindObject(h_gP->GetName()))->Fill( (*glv).P() );
    ((TH1F*)gDirectory->FindObject(h_gPt->GetName()))->Fill( (*glv).Pt() );
    ((TH1F*)gDirectory->FindObject(h_gtrkRadDHC1->GetName()))->Fill( gtrkRadDHC1 );

    ((TH1F*)gDirectory->FindObject(h_pchlvE->GetName()))->Fill( (*pchlv).E() );
    ((TH1F*)gDirectory->FindObject(h_pchlvP->GetName()))->Fill( (*pchlv).P() );
    ((TH1F*)gDirectory->FindObject(h_pchlvPt->GetName()))->Fill( (*pchlv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
    ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_pi0E->GetName()))->Fill(  pi0lv->E() );
    ((TH1F*)gDirectory->FindObject(h_pi0P->GetName()))->Fill(  pi0lv->P() );
    ((TH1F*)gDirectory->FindObject(h_pi0Pt->GetName()))->Fill( pi0lv->Pt() );
    ((TH1F*)gDirectory->FindObject(h_pi0Pcm->GetName()))->Fill( pi0lvcm->P() );
    ((TH1F*)gDirectory->FindObject(h_mofpi0->GetName()))->Fill( (*pi0lv).M() );
    ((TH1F*)gDirectory->FindObject(h_mofpi0diff->GetName()))->Fill( (*pi0lvcm).M()-(*pi0lv).M() );

    ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);
    ((TH2F*)gDirectory->FindObject(h_pi0_dalitz->GetName()))->Fill(msq12, msq23);

    ((TH1F*)gDirectory->FindObject(h_pieeE->GetName()))->Fill(pieelv->E());
    ((TH1F*)gDirectory->FindObject(h_pieeP->GetName()))->Fill(pieelv->P());
    ((TH1F*)gDirectory->FindObject(h_pieePt->GetName()))->Fill(pieelv->Pt());
    ((TH1F*)gDirectory->FindObject(h_pieeM->GetName()))->Fill(pieelv->M());

    ((TH1F*)gDirectory->FindObject(h_gpieeE->GetName()))->Fill(  ((*glv)+(*pieelv)).E() );
    ((TH1F*)gDirectory->FindObject(h_gpieeP->GetName()))->Fill(  ((*glv)+(*pieelv)).P() );
    ((TH1F*)gDirectory->FindObject(h_gpieePt->GetName()))->Fill( ((*glv)+(*pieelv)).Pt());
    ((TH1F*)gDirectory->FindObject(h_gpieeM->GetName()))->Fill(  ((*glv)+(*pieelv)).M() );

    ((TH1F*)gDirectory->FindObject(h_KE->GetName()))->Fill(Klv->E());
    ((TH1F*)gDirectory->FindObject(h_KP->GetName()))->Fill(Klv->P());
    ((TH1F*)gDirectory->FindObject(h_KPcm->GetName()))->Fill( KlvKcm->P() );
    ((TH1F*)gDirectory->FindObject(h_KPt->GetName()))->Fill(Klv->Pt());
    ((TH1F*)gDirectory->FindObject(h_mofK->GetName()))->Fill(Klv->M());


    ((TH1F*)gDirectory->FindObject(h_eldistDCH1plane->GetName()))->Fill( eldistDCH1plane );
    //    ((TH2F*)gDirectory->FindObject(h_eldistDCH1_ee_vz->GetName()))->Fill( eldistDCH1plane, vertex[2] );
    //    ((TH2F*)gDirectory->FindObject(h_eldistDCH1_pch_vz->GetName()))->Fill( eldistDCH1plane, vertex_e1pch[2] );
    ((TH1F*)gDirectory->FindObject(h_eldistLKrplane->GetName()))->Fill( eldistLKrplane );
    ((TH1F*)gDirectory->FindObject(h_gel1distLKrplane->GetName()))->Fill( gel1distLKrplane );
    ((TH1F*)gDirectory->FindObject(h_gel2distLKrplane->GetName()))->Fill( gel2distLKrplane );
    ((TH1F*)gDirectory->FindObject(h_pche1distDHC1plane->GetName()))->Fill( pche1distDHC1plane );
    ((TH1F*)gDirectory->FindObject(h_pche2distDHC1plane->GetName()))->Fill( pche2distDHC1plane );
    ((TH1F*)gDirectory->FindObject(h_pche1distLKrplane->GetName()))->Fill( pche1distLKrplane );
    ((TH1F*)gDirectory->FindObject(h_pche2distLKrplane->GetName()))->Fill( pche2distLKrplane );
    ((TH1F*)gDirectory->FindObject(h_gpchdistLKrplane->GetName()))->Fill( gpchdistLKrplane );

    ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill(elp1_uncorr[0]);
    ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill(elp1_uncorr[1]);
    ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill(elv1_uncorr[0]);
    ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill(elv1_uncorr[1]);
    ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill(elp1_uncorr[0], elp1_uncorr[1] );
    ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill(elv1_uncorr[0], elv1_uncorr[1] );

    ((TH1F*)gDirectory->FindObject(h_bx_pch->GetName()))            ->Fill(pchp_uncorr[0]);
    ((TH1F*)gDirectory->FindObject(h_by_pch->GetName()))            ->Fill(pchp_uncorr[1]);
    ((TH1F*)gDirectory->FindObject(h_bdxdz_pch->GetName()))            ->Fill(pchv_uncorr[0]);
    ((TH1F*)gDirectory->FindObject(h_bdydz_pch->GetName()))            ->Fill(pchv_uncorr[1]) ;
    ((TH1F*)gDirectory->FindObject(h_bx_vs_by_pch->GetName()))        ->Fill(pchp_uncorr[0], pchp_uncorr[1]);
    ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_pch->GetName()))->Fill(pchv_uncorr[0], pchv_uncorr[1]);

    ((TH1F*)gDirectory->FindObject(h_cda_ee->GetName()))->Fill(cda);
    ((TH1F*)gDirectory->FindObject(h_vx->GetName()))->Fill(vertex[0]);
    ((TH1F*)gDirectory->FindObject(h_vy->GetName()))->Fill(vertex[1]);
    ((TH1F*)gDirectory->FindObject(h_vz->GetName()))->Fill(vertex[2]);
    ((TH2F*)gDirectory->FindObject(h_vx_vs_vy->GetName()))->Fill(vertex[0],vertex[1]);

    ((TH1F*)gDirectory->FindObject(h_cda_e1pch->GetName()))->Fill(cda_e1pch);
    ((TH1F*)gDirectory->FindObject(h_vx_e1pch->GetName()))->Fill(vertex_e1pch[0]);
    ((TH1F*)gDirectory->FindObject(h_vy_e1pch->GetName()))->Fill(vertex_e1pch[1]);
    ((TH1F*)gDirectory->FindObject(h_vz_e1pch->GetName()))->Fill(vertex_e1pch[2]);
    ((TH2F*)gDirectory->FindObject(h_vx_vs_vy_e1pch->GetName()))->Fill(vertex_e1pch[0],vertex_e1pch[1]);

    ((TH1F*)gDirectory->FindObject(h_vz_e1pch_diff->GetName()))->Fill( vertex[2] - vertex_e1pch[2] );

    ((TH1F*)gDirectory->FindObject(h_angle_e1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_angle_ph_e1e2->GetName()))->Fill( (*eelv).Vect().Angle(glv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_angle_pch_pi0->GetName()))->Fill( (*pi0lv).Vect().Angle(pchlv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_a_e1e2_phe1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()), (*eelv).Vect().Angle(glv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_a_e1e2_pchpi0->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()), (*pi0lv).Vect().Angle(pchlv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_a_phe1e2_pchpi0->GetName()))->Fill( (*eelv).Vect().Angle(glv->Vect()), (*pi0lv).Vect().Angle(pchlv->Vect()) );

    //    ((TH1F*)gDirectory->FindObject(h_angle_pch_pi0_Kcm->GetName()))->Fill( pch3vKcm.Angle(pi03vKcm) );


    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"pi_eop");
    cutcounter++;




    // ~~~ CUT ~~~ CUT ~~~ CUT ~~~

    if((*pi0lv).Vect().Angle(pchlv->Vect())<0.012) {cleanup(); return 0;}
    //
    dirs[cutcounter-1]->cd();

    ((TH1I*)gDirectory->FindObject(h_evttypebytrks->GetName()))->Fill(evttypebytrks);

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

    ((TH1F*)gDirectory->FindObject(h_pchtrktime->GetName()))        ->Fill( pchtrktime );
    ((TH1F*)gDirectory->FindObject(h_pchtrkhodtime->GetName()))        ->Fill( pchtrkhodtime );
    ((TH1F*)gDirectory->FindObject(h_pchcltime->GetName()))            ->Fill( pchcltime );
    ((TH1F*)gDirectory->FindObject(h_pchtrktimediff->GetName()))    ->Fill( pchtrktime-pchtrkhodtime );
    ((TH1F*)gDirectory->FindObject(h_pchtrkcltimediff->GetName()))    ->Fill( pchtrktime-pchcltime );

    for (i=0; i<ntrack; i++) {
        p_track        = sevt->track[i].p;
        icl_track    = sevt->track[i].iClus;
        e_track      = sevt->cluster[icl_track].energy;
        eop_track      = e_track / p_track;
        ((TH1F*)gDirectory->FindObject(h_eop->GetName()))->Fill(eop_track);
        ((TH2F*)gDirectory->FindObject(h_eop_vs_p->GetName()))->Fill(eop_track, p_track);
    }

    ((TH1F*)gDirectory->FindObject(h_e1p->GetName()))->Fill( el1lv->P() );
    ((TH1F*)gDirectory->FindObject(h_e2p->GetName()))->Fill( el2lv->P() );
    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el1lv->E() );
    ((TH1F*)gDirectory->FindObject(h_eE->GetName()))->Fill( el2lv->E() );
    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el1lv->P() );
    ((TH1F*)gDirectory->FindObject(h_eP->GetName()))->Fill( el2lv->P() );
    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el1lv->Pt() );
    ((TH1F*)gDirectory->FindObject(h_ePt->GetName()))->Fill( el2lv->Pt() );

    ((TH1F*)gDirectory->FindObject(h_gE->GetName()))->Fill( (*glv).E() );
    ((TH1F*)gDirectory->FindObject(h_gP->GetName()))->Fill( (*glv).P() );
    ((TH1F*)gDirectory->FindObject(h_gPt->GetName()))->Fill( (*glv).Pt() );
    ((TH1F*)gDirectory->FindObject(h_gtrkRadDHC1->GetName()))->Fill( gtrkRadDHC1 );

    ((TH1F*)gDirectory->FindObject(h_pchlvE->GetName()))->Fill( (*pchlv).E() );
    ((TH1F*)gDirectory->FindObject(h_pchlvP->GetName()))->Fill( (*pchlv).P() );
    ((TH1F*)gDirectory->FindObject(h_pchlvPt->GetName()))->Fill( (*pchlv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_elsE->GetName()))->Fill( (*eelv).E() );
    ((TH1F*)gDirectory->FindObject(h_elsP->GetName()))->Fill( (*eelv).P() );
    ((TH1F*)gDirectory->FindObject(h_elsPt->GetName()))->Fill( (*eelv).Pt() );

    ((TH1F*)gDirectory->FindObject(h_pi0E->GetName()))->Fill(  pi0lv->E() );
    ((TH1F*)gDirectory->FindObject(h_pi0P->GetName()))->Fill(  pi0lv->P() );
    ((TH1F*)gDirectory->FindObject(h_pi0Pt->GetName()))->Fill( pi0lv->Pt() );
    ((TH1F*)gDirectory->FindObject(h_pi0Pcm->GetName()))->Fill( pi0lvcm->P() );
    ((TH1F*)gDirectory->FindObject(h_mofpi0->GetName()))->Fill( (*pi0lv).M() );
    ((TH1F*)gDirectory->FindObject(h_mofpi0diff->GetName()))->Fill( (*pi0lvcm).M()-(*pi0lv).M() );

    ((TH1F*)gDirectory->FindObject(h_mofpi0els->GetName()))->Fill(eelvmass);
    ((TH2F*)gDirectory->FindObject(h_pi0_dalitz->GetName()))->Fill(msq12, msq23);

    ((TH1F*)gDirectory->FindObject(h_pieeE->GetName()))->Fill(pieelv->E());
    ((TH1F*)gDirectory->FindObject(h_pieeP->GetName()))->Fill(pieelv->P());
    ((TH1F*)gDirectory->FindObject(h_pieePt->GetName()))->Fill(pieelv->Pt());
    ((TH1F*)gDirectory->FindObject(h_pieeM->GetName()))->Fill(pieelv->M());

    ((TH1F*)gDirectory->FindObject(h_gpieeE->GetName()))->Fill(  ((*glv)+(*pieelv)).E() );
    ((TH1F*)gDirectory->FindObject(h_gpieeP->GetName()))->Fill(  ((*glv)+(*pieelv)).P() );
    ((TH1F*)gDirectory->FindObject(h_gpieePt->GetName()))->Fill( ((*glv)+(*pieelv)).Pt());
    ((TH1F*)gDirectory->FindObject(h_gpieeM->GetName()))->Fill(  ((*glv)+(*pieelv)).M() );

    ((TH1F*)gDirectory->FindObject(h_KE->GetName()))->Fill(Klv->E());
    ((TH1F*)gDirectory->FindObject(h_KP->GetName()))->Fill(Klv->P());
    ((TH1F*)gDirectory->FindObject(h_KPcm->GetName()))->Fill( KlvKcm->P() );
    ((TH1F*)gDirectory->FindObject(h_KPt->GetName()))->Fill(Klv->Pt());
    ((TH1F*)gDirectory->FindObject(h_mofK->GetName()))->Fill(Klv->M());


    ((TH1F*)gDirectory->FindObject(h_eldistDCH1plane->GetName()))->Fill( eldistDCH1plane );
    //    ((TH2F*)gDirectory->FindObject(h_eldistDCH1_ee_vz->GetName()))->Fill( eldistDCH1plane, vertex[2] );
    //    ((TH2F*)gDirectory->FindObject(h_eldistDCH1_pch_vz->GetName()))->Fill( eldistDCH1plane, vertex_e1pch[2] );
    ((TH1F*)gDirectory->FindObject(h_eldistLKrplane->GetName()))->Fill( eldistLKrplane );
    ((TH1F*)gDirectory->FindObject(h_gel1distLKrplane->GetName()))->Fill( gel1distLKrplane );
    ((TH1F*)gDirectory->FindObject(h_gel2distLKrplane->GetName()))->Fill( gel2distLKrplane );
    ((TH1F*)gDirectory->FindObject(h_pche1distDHC1plane->GetName()))->Fill( pche1distDHC1plane );
    ((TH1F*)gDirectory->FindObject(h_pche2distDHC1plane->GetName()))->Fill( pche2distDHC1plane );
    ((TH1F*)gDirectory->FindObject(h_pche1distLKrplane->GetName()))->Fill( pche1distLKrplane );
    ((TH1F*)gDirectory->FindObject(h_pche2distLKrplane->GetName()))->Fill( pche2distLKrplane );
    ((TH1F*)gDirectory->FindObject(h_gpchdistLKrplane->GetName()))->Fill( gpchdistLKrplane );

    ((TH1F*)gDirectory->FindObject(h_bx->GetName()))            ->Fill(elp1_uncorr[0]);
    ((TH1F*)gDirectory->FindObject(h_by->GetName()))            ->Fill(elp1_uncorr[1]);
    ((TH1F*)gDirectory->FindObject(h_bdxdz->GetName()))            ->Fill(elv1_uncorr[0]);
    ((TH1F*)gDirectory->FindObject(h_bdydz->GetName()))            ->Fill(elv1_uncorr[1]);
    ((TH2F*)gDirectory->FindObject(h_bx_vs_by->GetName()))        ->Fill(elp1_uncorr[0], elp1_uncorr[1] );
    ((TH2F*)gDirectory->FindObject(h_bdxdz_vs_bdydz->GetName()))->Fill(elv1_uncorr[0], elv1_uncorr[1] );

    ((TH1F*)gDirectory->FindObject(h_bx_pch->GetName()))            ->Fill(pchp_uncorr[0]);
    ((TH1F*)gDirectory->FindObject(h_by_pch->GetName()))            ->Fill(pchp_uncorr[1]);
    ((TH1F*)gDirectory->FindObject(h_bdxdz_pch->GetName()))            ->Fill(pchv_uncorr[0]);
    ((TH1F*)gDirectory->FindObject(h_bdydz_pch->GetName()))            ->Fill(pchv_uncorr[1]) ;
    ((TH1F*)gDirectory->FindObject(h_bx_vs_by_pch->GetName()))        ->Fill(pchp_uncorr[0], pchp_uncorr[1]);
    ((TH1F*)gDirectory->FindObject(h_bdxdz_vs_bdydz_pch->GetName()))->Fill(pchv_uncorr[0], pchv_uncorr[1]);

    ((TH1F*)gDirectory->FindObject(h_cda_ee->GetName()))->Fill(cda);
    ((TH1F*)gDirectory->FindObject(h_vx->GetName()))->Fill(vertex[0]);
    ((TH1F*)gDirectory->FindObject(h_vy->GetName()))->Fill(vertex[1]);
    ((TH1F*)gDirectory->FindObject(h_vz->GetName()))->Fill(vertex[2]);
    ((TH2F*)gDirectory->FindObject(h_vx_vs_vy->GetName()))->Fill(vertex[0],vertex[1]);

    ((TH1F*)gDirectory->FindObject(h_cda_e1pch->GetName()))->Fill(cda_e1pch);
    ((TH1F*)gDirectory->FindObject(h_vx_e1pch->GetName()))->Fill(vertex_e1pch[0]);
    ((TH1F*)gDirectory->FindObject(h_vy_e1pch->GetName()))->Fill(vertex_e1pch[1]);
    ((TH1F*)gDirectory->FindObject(h_vz_e1pch->GetName()))->Fill(vertex_e1pch[2]);
    ((TH2F*)gDirectory->FindObject(h_vx_vs_vy_e1pch->GetName()))->Fill(vertex_e1pch[0],vertex_e1pch[1]);

    ((TH1F*)gDirectory->FindObject(h_vz_e1pch_diff->GetName()))->Fill( vertex[2] - vertex_e1pch[2] );

    ((TH1F*)gDirectory->FindObject(h_angle_e1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_angle_ph_e1e2->GetName()))->Fill( (*eelv).Vect().Angle(glv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_angle_pch_pi0->GetName()))->Fill( (*pi0lv).Vect().Angle(pchlv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_a_e1e2_phe1e2->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()), (*eelv).Vect().Angle(glv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_a_e1e2_pchpi0->GetName()))->Fill( (*el1lv).Vect().Angle(el2lv->Vect()), (*pi0lv).Vect().Angle(pchlv->Vect()) );
    ((TH1F*)gDirectory->FindObject(h_a_phe1e2_pchpi0->GetName()))->Fill( (*eelv).Vect().Angle(glv->Vect()), (*pi0lv).Vect().Angle(pchlv->Vect()) );

    //    ((TH1F*)gDirectory->FindObject(h_angle_pch_pi0_Kcm->GetName()))->Fill( pch3vKcm.Angle(pi03vKcm) );

    file1->cd();
    h_cutflow->Fill(cutcounter);
    h_cutflow->GetXaxis()->SetBinLabel(cutcounter, cuts[cutcounter-1].c_str());//"pi_eop");
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
    if(pchlv)             delete pchlv;     pchlv =0;
    if(pieelv)             delete pieelv;    pieelv =0;
    if(glv)             delete glv;       glv =0;
    if(Klv)             delete Klv;       Klv =0;
    if(el1lvcm)         delete el1lvcm;   el1lvcm =0;
    if(el2lvcm)         delete el2lvcm;   el2lvcm =0;
    if(glvcm)             delete glvcm;     glvcm =0;
    if(pi0lvcm)         delete pi0lvcm;   pi0lvcm =0;
    if(pi0lv)             delete pi0lv;     pi0lv =0;
    if(pi0lvKcm)         delete pi0lvKcm;  pi0lvKcm =0;
    if(pchlvKcm)         delete pchlvKcm;  pchlvKcm =0;
    if(KlvKcm)             delete KlvKcm;    KlvKcm =0;

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
