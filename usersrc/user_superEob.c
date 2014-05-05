/******************************************************************************/
/* COmPACT user routine: user_superEob(superBurst *bur, superEndofBurst *eob) */
/*                                                                            */
/* User routine called at the end of every SuperCOmPACt burst and is          */
/*  passed the burst header and burst trailer.                                */
/*                                                                            */
/*                                                                BH 1/3/98   */
/******************************************************************************/

#include "cmpio.h"
#include "user.h"
#include "reader.h"

int user_superEob(superBurst *bur,superEndofBurst *eob) {
/* WARNING: do not alter things before this line */
/*---------- Add user C code here ----------*/
  printSuperEndofBurst(eob,fprt);
  //HistnoBursts->Fill(0);

  
  if (burstCounter == 200)
    {
      printf ("### writing histograms to root file\n");
      if (noBursts > burstCounter)
	file1->Write("",TObject::kOverwrite);
      else 
	file1->Write();
      burstCounter = 0;
    }
  
/*----------- End of user C code -----------*/
  return 0;
}
