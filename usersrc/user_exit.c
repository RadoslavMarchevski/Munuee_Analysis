/********************************************************/
/* COmPACT user routine: user_exit()                    */
/*                                                      */
/* User routine called once all data has been read. This*/
/* allows output files to be closed and any other user  */
/* resources to be tidied up.                           */
/*                                          RWM 20/6/97 */
/********************************************************/

#include <math.h>
#include "cmpio.h"
#include "user.h"
#include "reader.h"
#include "vscompact.h"

int user_exit() {
/* WARNING: do not alter things before this line */
/*---------- Add user C code here ----------*/
  

  file1->Write("",TObject::kOverwrite);
  //file1->Write();
  file1->Close();

  
  fclose(fprt);

/*----------- End of user C code -----------*/
  return 0;
}
