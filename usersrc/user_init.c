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




TH1F* HistnoBursts;
TH1F* HistnoGoodBursts;

TH1I* RejectedBursts;

float CPDpos_leftDownCorner[256][2];
float CELLpos_leftDownCorner[256][64][2];
float CPDineff[256][50];     // 0-9: total, 10-19: Ptrack 15-20GeV, 20-29: Ptrack 20-25GeV, 30-39: Ptrack 25-35GeV, 40-49: Ptrack 35-65GeV
float CELLineff[256][64][3]; // 0: E/p>0.6, 1: E/p 0.6-0.95, 2: E/p>1.1
int   CPDindex, CELLindex;
float CELLlength = 1.975;
float CPDlength = 8 * CELLlength;


int runNo;
int magnetCurrent;
int noBursts;
int burstCounter;
int IS_DATA;          // data type == 1
int IS_MC;            // data type == 2


// E/p corrections for each cell
FILE *EopCorrfile;
char  EopCorrfilename[200];
char rootfilename[200];
float EopCorr[100][256][64]; // ab v65: corrections for different (sub-) periods (up to 100)
int   periodFlag;            // ab v65: period to be defined in superBurst.c

Hist_dir* Initial_dir;
Hist_dir* K3pi_selection;
Hist_dir* dir1;
Hist_dir* dir2;

int user_init() {
/* WARNING: do not alter things before this line */
/*---------- Add user C code here ----------*/

    int   i, j, k;
    int   l, m, n;
    //Choise of directory name and type:
    //0 - Initial
    //1 - K3pi selection
    //2 - Kmunuee selection
    Initial_dir =  new  Hist_dir("Initial_histograms",0);
    K3pi_selection =  new  Hist_dir("K3pi selection",1);
    dir1 =  new  Hist_dir("Signal selection_Cut_1",2);
    dir2 =  new  Hist_dir("Signal selection_Cut_2",2);
    //dir2 =  new  Hist_dir("Test2");



    //###############################################################################################
    //
    //      Histogramme
    //
    //###############################################################################################





/*----------- End of user C code -----------*/
    return 0;
}
