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
#include "user_NEW.h"
#include "reader.h"
#include <constants.h>
#include "MC_Charged_Particle.h"
#include "Charged_Particle.h"
#include <iostream>

// using Anweisung
using std::cout;
using std::endl;
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
        Ktype = evt->decay.Ktype;

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
            True_Momentum[i].SetPxPyPzE(evt->part[i].p[1],evt->part[i].p[2],evt->part[i].p[3], evt->part[i].p[0]);
            //cout << "i = " << i << "Type = "<< Particle_type[i] << " Mass = " << TMath::Sqrt(evt->part[i].p[0]*evt->part[i].p[0] - evt->part[i].p[1]*evt->part[i].p[1] - evt->part[i].p[2]*evt->part[i].p[2] - evt->part[i].p[3]*evt->part[i].p[3])<< " P = " << True_Momentum[i].P() << endl;
        }
        //cout << "---------------------------------------------------" << endl;
        //mc_Three_Track_Momentum = True_Momentum[1] + True_Momentum[2] + True_Momentum[3]  ;
        //mc_Two_Track_Momentum   = True_Momentum[2] + True_Momentum[3] ;
        //    //Nu_true = (*true_V[0]) - (*true_V[1]) - (*true_V[2]) - (*true_V[3]) ;
    }

    user_superCmpEvent(sbur,&(evt->scmpevt));


    /*----------- End of user C code -----------*/
    return 0;
}
