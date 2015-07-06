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
#include "user_NEW.h"
//#include "user.h"
#include "reader.h"
#include "vscompact.h"
#include "Hist_dir.h"

int user_exit() {

/* WARNING: do not alter things before this line */
/*---------- Add user C code here ----------*/

    TFile* file = new TFile("run2003_2004.root", "RECREATE");
    Initial_dir->AddToFile(file);
    K3pi_selection->AddToFile(file);
    dir1->AddToFile(file);
    dir2->AddToFile(file);
    printf ("### saving histograms to file run2003_2004.root \n");
    //delete dir1;
    //delete dir2;
    file->Write();
    file->Close();

/*----------- End of user C code -----------*/
    return 0;
}
