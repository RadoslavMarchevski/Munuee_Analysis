/********************************************************/
/* COmPACT user routine: user_init()                    */
/*                                                      */
/* User routine called upon program startup to allow    */
/* initialization of the user files, variables etc.     */
/*                                          RWM 20/6/97 */
/********************************************************/

#include "cmpio.h"
//#include "user.h"
#include "user_NEW.h"
#include "reader.h"

#include "TFile.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TLorentzVector.h"

#include <stdio.h>
#include <string>
#include <vector>
#include <math.h>
#include "Hist_dir.h"

superMcEvent* mcevent;
superCmpEvent* sevt;



TH1I* HistnoBursts;
TH1I* HistnoGoodBursts;

TH1I* RejectedBursts;

float CPDpos_leftDownCorner[256][2];
float CELLpos_leftDownCorner[256][64][2];
float CPDineff[256][50];     // 0-9: total, 10-19: Ptrack 15-20GeV, 20-29: Ptrack 20-25GeV, 30-39: Ptrack 25-35GeV, 40-49: Ptrack 35-65GeV
float CELLineff[256][64][3]; // 0: E/p>0.6, 1: E/p 0.6-0.95, 2: E/p>1.1
int   CPDindex, CELLindex;
float CELLlength = 1.975;
float CPDlength = 8 * CELLlength;


int runNo;
int Ktype;
int PrevrunNo = 0;
int magnetCurrent;
int AchromatCurrent;
int noBursts;
int noBurstsPerRun;
int noSelectedEvents;
int burstCounter;
int IS_DATA;          // data type == 1
int IS_MC;            // data type == 2


// E/p corrections for each cell
FILE *EopCorrfile;
char  EopCorrfilename[200];
char rootfilename[200];
float EopCorr[100][256][64]; // ab v65: corrections for different (sub-) periods (up to 100)
int   periodFlag;            // ab v65: period to be defined in superBurst.c


int  Particle_type[20];
int  Npart;
float Particle_production_zvtx[20];
float Particle_decay_zvtx[20];
float  DKaon[3];
TLorentzVector True_Momentum[20];
//TLorentzVector mc_Three_Track_Momentum;
//TLorentzVector mc_Two_Track_Momentum;

Hist_dir* Initial_dir;
Hist_dir* K3pi_selection;
Hist_dir* dir1;
Hist_dir* dir2;
Hist_dir* dir3;
Hist_dir* dir4;
Hist_dir* dir5;
Hist_dir* dir6;
Hist_dir* dir7;
Hist_dir* dir8;
Hist_dir* dir9;
Hist_dir* dir10;
Hist_dir* dir11;
Hist_dir* dir12;
Hist_dir* dir13;
Hist_dir* dir14;
Hist_dir* dir15;
Hist_dir* dir16;
Hist_dir* dir17;
Hist_dir* dir18;
Hist_dir* dir19;
Hist_dir* MC_reweight;

// missing compact declarations


int user_init() {
/* WARNING: do not alter things before this line */
/*---------- Add user C code here ----------*/

    int   i, j, k;
    int   l, m, n;
    //Choise of directory name and type:
    //0 - Initial
    //1 - K3pi selection
    //2 - Kmunuee selection
    Initial_dir =  new  Hist_dir("No cuts",0);
    K3pi_selection =  new  Hist_dir("K3pi selection",1);
    dir1 =  new  Hist_dir("Munuee index matching",2);
    dir2 =  new  Hist_dir("WS no Pt and mu P cut",3);
    dir3 =  new  Hist_dir("Munuee Geometry Cut",2);
    dir4 =  new  Hist_dir("Munuee Momentum Cut",2);
    dir5 =  new  Hist_dir("WS Pt > 10 MeV",3);
    dir6 =  new  Hist_dir("WS Pt > 20 MeV",3);
    dir7 =  new  Hist_dir("Munuee e+e- Invariant Mass Cut",2);
    dir8 =  new  Hist_dir("Munuee Timing Cut",2);
    dir9 =  new  Hist_dir("WS Pt > 30 MeV",3);
    dir10=  new  Hist_dir("Munuee lda3 and missing mass Cut",2);
    dir11=  new  Hist_dir("Munuee z vtx and muon status Cut",2);
    dir12=  new  Hist_dir("Munuee accidentals",2);
    //dir13=  new  Hist_dir("Munuee in time",2);

    dir13=  new  Hist_dir("No Pt No mu P cuts",3);
    dir15=  new  Hist_dir("Pt 10 MeV",3);
    dir16=  new  Hist_dir("Pt 20 MeV",3);
    dir17=  new  Hist_dir("Pt 30 MeV",3);
    dir18=  new  Hist_dir("Kcharge = +-3",3);
    dir19=  new  Hist_dir("WS Accidentals",3);
    MC_reweight =  new  Hist_dir("MC after EoP and lda3 reweighting",3);

    HistnoBursts                 = new TH1I("noBursts", "number of bursts read",2000,15000,17000);
    HistnoGoodBursts             = new TH1I("noGoodBursts", "number of good bursts read",2000,15000,17000);


/*----------- End of user C code -----------*/
    return 0;
}
