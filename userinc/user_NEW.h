extern char gString[50];

extern float CPDpos_leftDownCorner[256][2];
extern float CELLpos_leftDownCorner[256][64][2];
extern int   CPDindex, CELLindex;
extern float CPDlength, CELLlength;

extern int runNo;
extern int magnetCurrent;
extern int noBursts;
extern int burstCounter;
extern int IS_DATA, IS_MC;



// missing compact declarations
extern "C" int accep_(int* run, float* x, float* y);
extern "C" int closap_double_(double p1[3], double p2[3], double v1[3], double v2[3], double *dmin, double vtx[3]);
extern "C" void blue_tack_(int *nchar, float *tmom, float Vxyz[3], float vpnt[2], float vdir[2]);

//extern int LKr_acc(int nrun, double pos_x, double pos_y, float par);
extern int LKr_acc(int nrun, float pos_x, float pos_y, float par);        // offizielle Routine in lkraccep_2007.c benutzt ab v26

extern void GetCpdCellIndex(float pos_x, float pos_y, int *cpd_index, int *cell_index);

// E/p corrections for each cell
extern FILE *EopCorrfile;
extern char  EopCorrfilename[200];
extern float EopCorr[100][256][64]; // ab v65: corrections for different (sub-) periods (up to 100)
extern int   periodFlag;            // ab v65: period to be defined in superBurst.c


#define SQR(x) ((x)*(x))

/*************************************************************************/
/*                        root include files                             */
/*************************************************************************/
#include "TROOT.h"
#include "TMath.h"
#include "TH3.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TList.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TDirectory.h"
#include <stdio.h>
#include <string>
#include <vector>
#include <math.h>
<<<<<<< HEAD
//#include "MC_Charged_Particle.h"

using namespace std;
extern superMcEvent* mcevent;
extern superCmpEvent* sevt;
extern int Particle_type[20];
extern int Npart;
extern float Particle_production_zvtx[20];
extern float Particle_decay_zvtx[20];
extern float  DKaon[3];
extern TLorentzVector True_Momentum[20];
//extern TLorentzVector mc_Three_Track_Momentum;
//extern TLorentzVector mc_Two_Track_Momentum;
=======
using namespace std;
extern superMcEvent* mcevent;
//extern superCmpEvent* sevt;
>>>>>>> 434503a23b9e52f2ae2a1b9612608c6ffdc3ceeb
