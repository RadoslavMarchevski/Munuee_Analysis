#include "compactF77.h"
#include "CompactGeomF77.h"
      FUNCTION lda3(ISE,QSE,itr)
C
C  B.Bloch may 2005
C  return a linear combination of 3 variables (E/p, radial rms and radial
C  normalized distance track to cluster ) for a given track itr in Supercompact
C  usage : out_lda = lda3(ISE,QSE,itr)
C  where itr is the track number ( fortran convention, 1.....n)
C  if no LKr cluster is associated to the track, it returns 0.
C  user has to decide on the cut value applied to the lda output,according
C  to the performances he wants to obtain. Have a look at my web page for
C  plots and numbers !
C  http://aleph.web.cern.ch/aleph/~bbloch/NA48/ke4/
C  for an electron efficiency over 95%, I suggest a value> 0.95
C  the pion effciency on tracks with E/P>0.9 is then a few percents
C
      IMPLICIT NONE      
      INTEGER ISB(2),ISE(2)    /* fake array sizes since C structs */
      REAL*4  QSB(2),QSE(2)    /* fake array sizes since C structs */
#include "CompactCommons.h"
#include "rlibFunctionsF77.inc"
C
      REAL*4 lda3,SIGMOIDL
      external SIGMOIDL
      integer itr,iclu
      REAL*4 geo_lkrZ,Zdch,pit,eovp,rmsr,xlkr,ylkr,distx,disty,disn
      REAL*4 RIN1,out1

      geo_lkrZ = GEOM_LKR_Z
      Zdch     = GEOM_DCH_PLANE_Z(4,1)


      lda3 = 0.
      if (itr.le.0) go to 100
      iclu = STRACK_ICLUS(itr)+1
      if(iclu.le.0) go to 100
      pit = p_corr_ab(STRACK_P(itr),STRACK_Q(itr))
      eovp = SCLUSTER_ENERGY(iclu)/pit
      rmsr   = sqrt(SCLUSTER_RMSX(iclu)*SCLUSTER_RMSX(iclu) +
     &              SCLUSTER_RMSY(iclu)*SCLUSTER_RMSY(iclu))
      xlkr   = STRACK_X(itr) + STRACK_DXDZ(itr)*(geo_lkrZ-Zdch)
      ylkr   = STRACK_Y(itr) + STRACK_DYDZ(itr)*(geo_lkrZ-Zdch)
      distx  = xlkr - SCLUSTER_X(iclu)
      disty  = ylkr - SCLUSTER_Y(iclu)
      disn = sqrt(pit*(distx*distx+disty*disty))

C
      RIN1 = 39.49729
     > -82.20312 * eovp
     > +32.76240 * rmsr
     > +1.597024 * disn
C
      OUT1 = SIGMOIDL(RIN1)
C
      lda3 = 1.000964e+00
     > +(-9.863838e-01) * OUT1
C
 100  return
      END

C******************************************************************
      REAL FUNCTION SIGMOIDL(X)
      IF(X.GT.37.) THEN
         SIGMOIDL = 1.
      ELSE IF(X.LT.-37.) THEN
         SIGMOIDL = 0.
      ELSE
         SIGMOIDL = 1./(1.+EXP(-X))
      ENDIF
      END
