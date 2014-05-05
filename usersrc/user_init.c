/********************************************************/
/* COmPACT user routine: user_init()                    */
/*                                                      */
/* User routine called upon program startup to allow    */
/* initialization of the user files, variables etc.     */
/*                                          RWM 20/6/97 */
/********************************************************/

#include "cmpio.h"
#include "user.h"
#include "reader.h"

superMcEvent* mcevent;

char rootfilename[200];
char rootFileName[200];

char histoname[200];
char histoname2[200];

TFile *file1;

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
float EopCorr[100][256][64]; // ab v65: corrections for different (sub-) periods (up to 100)
int   periodFlag;            // ab v65: period to be defined in superBurst.c


int user_init() {
/* WARNING: do not alter things before this line */
/*---------- Add user C code here ----------*/

  int   i, j, k;
  int   l, m, n;

  
  fprt=fopen("compact.txt","w");

  if (strcmp (rootfilename, "") == 0)
    {
      sprintf (rootFileName, "run2007.root");
    }
  else
    {
      sprintf (rootFileName, "%s.root", rootfilename);
    }
  printf ("### saving histograms to file %s\n", rootFileName);

  file1 = new TFile(rootFileName, "RECREATE");
  //  file1 = new TFile("run2007.root", "RECREATE");
  
  
  //###############################################################################################
  //
  //      Histogramme
  //
  //###############################################################################################


 HistnoBursts = new TH1F("HistonoBursts","Number of read bursts",5,0,5);


  noBursts = 0;
  burstCounter = 0;

  
  for (i=0; i<16; i++)
    for (j=0; j<16; j++)
      {
	k = i*16 + j;	
	CPDpos_leftDownCorner[k][0] = (-1)*(7-i)*CPDlength;  // LKr-Bezugssystem ist linksh�ndig
	CPDpos_leftDownCorner[k][1] = (7-j)*CPDlength;
	//printf ("CPD %d: position left down corner = %.2f, \t%.2f\n", k, CPDpos_leftDownCorner[k][0], CPDpos_leftDownCorner[k][1]);	      

	for (m=0; m<8; m++)
	  for (n=0; n<8; n++)
	    {
	      l = m*8 + n;
	      CELLpos_leftDownCorner[k][l][0] = CPDpos_leftDownCorner[k][0] - (7-m)*CELLlength;
	      CELLpos_leftDownCorner[k][l][1] = CPDpos_leftDownCorner[k][1] + (7-n)*CELLlength;
	      //printf ("CELL %d in CPD %d: position left down corner = %.2f, \t%.2f\n", 
	      //      l, k, CELLpos_leftDownCorner[k][l][0], CELLpos_leftDownCorner[k][l][1]);
	    }
      }
  ////////////////////////////////////
  //  E/p corrections for each cell //
  ////////////////////////////////////
  int    cpd, cell;

  // Period 5, Runs 20410-20415 (= period 51)
  sprintf (EopCorrfilename, "files/eopCorrfile_p5_20410-20415_v63.dat");
  EopCorrfile = fopen (EopCorrfilename, "r");
  printf ("### Period 51 = Runs 20410-20415: Reading E/p corrections for each cell from file %s\n", EopCorrfilename);  
  
  for (i=0; i<256; i++)
    for (j=0; j<64; j++)
      {	    
	fscanf (EopCorrfile, "%i %i %f\n", &cpd, &cell, &EopCorr[51][i][j]);	    
	//printf ("cpd%i\tcell%i\t%f\n", cpd, cell, EopCorr[51][i][j]);
      }

  // Period 5, Runs 20416-20485 (= period 52)
  sprintf (EopCorrfilename, "files/eopCorrfile_p5_20416-20485_v63.dat");
  EopCorrfile = fopen (EopCorrfilename, "r");
  printf ("### Period 52 = Runs 20416-20485: Reading E/p corrections for each cell from file %s\n", EopCorrfilename);  
  
  for (i=0; i<256; i++)
    for (j=0; j<64; j++)
      {	    
	fscanf (EopCorrfile, "%i %i %f\n", &cpd, &cell, &EopCorr[52][i][j]);	    
	//printf ("cpd%i\tcell%i\t%f\n", cpd, cell, EopCorr[52][i][j]);
      }


/*----------- End of user C code -----------*/
  return 0;
}
