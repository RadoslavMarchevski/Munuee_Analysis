/**************************************************************/
/* COmPACT user routine: user_superMcEvent(superMcEvent *evt) */
/*                                                            */
/* User routine called everytime an event `*evt' is           */
/* loaded. A return value of greater than zero denotes        */
/* an error condition has occured.                            */
/*                                   BH 13/2/98   RWM 11/7/97 */
/**************************************************************/

#include <TMath.h>
#include <math.h>
#include "cmpio.h"
#include "user.h"
#include "reader.h"
#include <constants.h>

// Datei iostream einbinden
#include <iostream>

// using Anweisung
using std::cout;
using std::endl;
int user_superMcEvent(superBurst *sbur,superMcEvent *evt) {
  /* WARNING: do not alter things before this line */
  /*---------- Add user C code here ----------*/
  static int nuserevt=0;
  if(nuserevt<20)
    {
      printSuperMcEvent(evt,fprt);
      nuserevt++;
    }



  // int *Ptype;
  // float *pvtx,*dvtx;
  // float *vtxdiff;

//  if (IS_MC)
//    {
//      k_dvtx[0]=evt->decay.dvertex[0];
//      k_dvtx[1]=evt->decay.dvertex[1];
//      k_dvtx[2]=evt->decay.dvertex[2];
//      Npart = evt->Npart;
//
//      for (int i=0;i<Npart;i++){
//	pType[i]=evt->part[i].type;
//	pvtx[i]=evt->part[i].pvertex[2];
//	dvtx[i]=evt->part[i].dvertex[2];
//	true_V[i] = new TLorentzVector;
//	true_V[i]->SetPxPyPzE(evt->part[i].p[1],evt->part[i].p[2],evt->part[i].p[3],evt->part[i].p[0]);
//      }
//      Muee_three_track =  (*true_V[1]) + (*true_V[2]) + (*true_V[3])  ;
//      mc_ee = (*true_V[2]) + (*true_V[3]) ;
//      Nu_true = (*true_V[0]) - (*true_V[1]) - (*true_V[2]) - (*true_V[3]) ;
//    }
//
  user_superCmpEvent(sbur,&(evt->scmpevt));


  /*----------- End of user C code -----------*/
  return 0;
}
