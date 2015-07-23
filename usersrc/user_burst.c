/********************************************************/
/* COmPACT user routine: user_burst(Burst *bur)         */
/*                                                      */
/* User routine called everytime a burst `*bur' is      */
/* loaded. A return value of greater than zero denotes  */
/* an error condition has occured.                      */
/*                                          RWM 20/6/97 */
/********************************************************/

#include "cmpio.h"
<<<<<<< HEAD
#include "user_NEW.h"
=======
#include "user.h"
>>>>>>> 434503a23b9e52f2ae2a1b9612608c6ffdc3ceeb
#include "reader.h"
int user_burst(Burst *bur) {
/* WARNING: do not alter things before this line */
/*---------- Add user C code here ----------*/
  static int Nburst=0;
  static int RunCurrent=0;
<<<<<<< HEAD
  // if(!Nburst) PrintCmpGeom(Geom, fprt);
=======
  if(!Nburst) PrintCmpGeom(Geom, fprt);
>>>>>>> 434503a23b9e52f2ae2a1b9612608c6ffdc3ceeb
  Nburst++;
/*  if(Nburst++ < 5) printBurst(bur,fprt);*/
/******************************************************************/
/* The following lines show how to call the printing routine      */
/* for the list of LKR dead cells, for the current run            */
/******************************************************************/
  if(bur->nrun != RunCurrent)
    {
      RunCurrent=bur->nrun;
<<<<<<< HEAD
      //  CmpLkrDeadPrint(fprt);
=======
      CmpLkrDeadPrint(fprt);
>>>>>>> 434503a23b9e52f2ae2a1b9612608c6ffdc3ceeb
/******************************************************************/
/* The following lines show how to call the printing routine      */
/* for the list of timing offsets, for the current run            */
/******************************************************************/
<<<<<<< HEAD
      //CmpTimeOffsetPrint(bur,fprt);
=======
      CmpTimeOffsetPrint(bur,fprt);
>>>>>>> 434503a23b9e52f2ae2a1b9612608c6ffdc3ceeb
    }
/* To print out some of the database data */
/*
  printf(" From Burst nBur=%d nBurBad=%d \n",rdb->nBur, rdb->nBurBad);
  printf(" From Burst BURST=%4d %4d %4d %4d \n",
<<<<<<< HEAD
         bdb->nFiltCh, bdb->nGoodCh, bdb->nFiltNe,  bdb->nGoodNe);
=======
         bdb->nFiltCh, bdb->nGoodCh, bdb->nFiltNe,  bdb->nGoodNe);   
>>>>>>> 434503a23b9e52f2ae2a1b9612608c6ffdc3ceeb
*/
/*----------- End of user C code -----------*/
  return 0;
}
