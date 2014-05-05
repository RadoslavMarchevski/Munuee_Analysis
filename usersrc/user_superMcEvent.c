/**************************************************************/
/* COmPACT user routine: user_superMcEvent(superMcEvent *evt) */
/*                                                            */
/* User routine called everytime an event `*evt' is           */
/* loaded. A return value of greater than zero denotes        */
/* an error condition has occured.                            */
/*                                   BH 13/2/98   RWM 11/7/97 */
/**************************************************************/

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

   if (IS_MC)
     {


	   mcevent = evt;
       // MC true values for acceptance calculation  (ab v551)
//       zVertexTrue = evt->decay.dvertex[2];
//       pTrackTrue = sqrt( SQR(evt->part[1].p[1]) + SQR(evt->part[1].p[2]) + SQR(evt->part[1].p[3]) );

     }
   
  // printf ("### is MC event - calling user_superCmpEvent ###\n");
  user_superCmpEvent(sbur,&(evt->scmpevt));

  
  /*----------- End of user C code -----------*/
  return 0;
}
