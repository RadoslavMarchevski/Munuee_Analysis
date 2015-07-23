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
<<<<<<< HEAD
#include "user_NEW.h"
#include "reader.h"
#include <constants.h>
#include "MC_Charged_Particle.h"
#include "Charged_Particle.h"
=======
#include "user.h"
#include "reader.h"
#include <constants.h>

// Datei iostream einbinden
>>>>>>> 434503a23b9e52f2ae2a1b9612608c6ffdc3ceeb
#include <iostream>

// using Anweisung
using std::cout;
using std::endl;
<<<<<<< HEAD
int user_superMcEvent(superBurst *sbur, superMcEvent *evt) {
    /* WARNING: do not alter things before this line */
    /*---------- Add user C code here ----------*/
    static int nuserevt=0;
    //if(nuserevt<20)
    //{
    //printSuperMcEvent(evt,fprt);
    //  nuserevt++;
    //}
    if (IS_MC)
    {

        Npart = evt->Npart;

        //int  Particle_type[20];
        //float Particle_production_zvtx[20];
        //float Particle_decay_zvtx[20];
        //float  DKaon[3];

        DKaon[0] = evt->decay.dvertex[0];
        DKaon[1] = evt->decay.dvertex[1];
        DKaon[2] = evt->decay.dvertex[2];
        //
        for (int i=0;i<Npart;i++){
            Particle_type[i]           =evt->part[i].type;
            Particle_production_zvtx[i]=evt->part[i].pvertex[2];
            Particle_decay_zvtx[i]     =evt->part[i].dvertex[2];
            True_Momentum[i].SetPxPyPzE(evt->part[i].p[1],evt->part[i].p[2],evt->part[i].p[3],evt->part[i].p[0]);
        }
        //mc_Three_Track_Momentum = True_Momentum[1] + True_Momentum[2] + True_Momentum[3]  ;
        //mc_Two_Track_Momentum   = True_Momentum[2] + True_Momentum[3] ;
        //    //Nu_true = (*true_V[0]) - (*true_V[1]) - (*true_V[2]) - (*true_V[3]) ;
    }

    user_superCmpEvent(sbur,&(evt->scmpevt));


    /*----------- End of user C code -----------*/
    return 0;
=======
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
>>>>>>> 434503a23b9e52f2ae2a1b9612608c6ffdc3ceeb
}
