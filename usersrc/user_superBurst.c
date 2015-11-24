/***********************************************************/
/* COmPACT user routine: user_superBurst(superBurst *sbur) */
/*                                                         */
/* User routine called everytime a SuperCOmPACT burst      */
/* `*bur' is loaded.                                       */
/* A return value of greater than zero denotes an error    */
/* condition has occured.                                  */
/*                                             BH 2/3/98   */
/***********************************************************/

#include "cmpio.h"
//#include "user.h"
#include "user_NEW.h"
#include "reader.h"

int user_superBurst(superBurst *sbur) {
/* WARNING: do not alter things before this line */
/*---------- Add user C code here ----------*/
    sbur->BadB.Skip = 0; /* see user_superBurst.example.c to learn to use it */

    noBursts++;            // total number of bursts
    burstCounter++;        // counter set to 0 in superEob if certain number of bursts reached to write the root file
    magnetCurrent =  sbur->MNP33current;
    runNo = sbur->nrun;

    if (sbur->BadB.Phys!=0)    // exclude bad physics bursts
        sbur->BadB.Skip = 1;
    if (sbur->BadB.Dch!=0)    // exclude bad DCH bursts
        sbur->BadB.Skip = 1;

    if(sbur->BadB.Skip>0) return 0;
    //if (sbur->time != 1190177600)
    //sbur->BadB.Skip = 1;

    // Ab v65: Determine in which period we are and set flag
    //

    // Determine Data type
    IS_DATA = false;
    IS_MC = false;
    if (sbur->brtype == 1)
        IS_DATA = true;
    else if (sbur->brtype == 2)
        IS_MC = true;
    //printf ("IS_DATA = %i \t IS_MC = %i\n", IS_DATA, IS_MC);
    if (IS_DATA){
        if( sbur->BadB.Lkr == 1||
            sbur->BadB.Mbx == 1||
            sbur->BadB.Muv == 1||
            sbur->BadB.HodC== 1||
            sbur->BadB.Dch == 1){
            sbur->BadB.Skip=1;
        }
        if(sbur->nrun > 15900 && sbur->nrun < 16121 ){
            sbur->BadB.Skip=1;
        }

    }

    //  printf ("burst number = %d  -  run %d  -  BadB.Phys = %i\n", noBursts, runNo, sbur->BadB.Phys);
    printf ("burst number = %d  -  run %d  -  period %d  -  burst time %d  -  BadB.Skip = %i  -  magnet current = %.1f  -  dataType = %i\n",
            noBursts, runNo, periodFlag, sbur->time, sbur->BadB.Skip, (float)(magnetCurrent/1000.), sbur->brtype);
    //printf ("magnet current = %d\n", magnetCurrent);

    // HistnoBursts->Fill(0);
    //if (sbur->BadB.Skip == 0)
    //HistnoGoodBursts->Fill(0);

    //  printf ("n2track_total = %i \t n2track_total_cut05 = %i \t nKe3_total = %i \t nKe3_total_psum60_120 = %i\n",
    //  n2track_total, n2track_total_cut05, nKe3_total, nKe3_total_psum60_120);

/*----------- End of user C code -----------*/
    return 0;
}
