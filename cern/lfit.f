*
* $Id: lfit.F,v 1.1.1.1 1996/04/01 15:02:28 mclareni Exp $
*
* $Log: lfit.F,v $
* Revision 1.1.1.1  1996/04/01 15:02:28  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_FORTRAN)||!defined(CERNLIB_CDC)
      SUBROUTINE LFIT(X,Y,L,KEY,A,B,E)
C
C     TO FIT A STRAIGHT LINE    Y=A*X+B    TO L POINTS WITH ERROR E
C     SEE MENZEL , FORMULAS OF PHYSICS P.116
C     POINTS WITH Y=0 ARE IGNOERD IF KEY=0
C     L IS NO. OF POINTS
C
      REAL*4 X,Y,A,B,E
      DIMENSION X(1),Y(1)
C
C     CALCULATE SUMS
      IF(L-2) 25,1,1
    1 COUNT=0.0
      SUMX=0.0
      SUMY=0.0
      SUMXY=0.0
      SUMXX=0.0
      SUMYY=0.0
      DO 10 J=1,L
      IF(Y(J).EQ.0..AND.KEY.EQ.0) GO TO 10
      SUMX=SUMX+X(J)
      SUMY=SUMY+Y(J)
      COUNT=COUNT+1.0
   10 CONTINUE
      IF(COUNT.LE.1.) GO TO 25
      YMED=SUMY/COUNT
      XMED=SUMX/COUNT
      DO 20 J=1,L
      IF(Y(J).EQ.0..AND.KEY.EQ.0) GO TO 20
      SCARTX=X(J)-XMED
      SCARTY=Y(J)-YMED
      SUMXY=SUMXY+SCARTX   *SCARTY
      SUMXX=SUMXX+SCARTX   *SCARTX
      SUMYY=SUMYY+SCARTY   *SCARTY
   20 CONTINUE
C
C     FIT PARAMETERS
      IF(SUMXX.EQ.0.) GO TO 25
      A=SUMXY/SUMXX
      B=YMED-A*XMED
      IF(COUNT.LT.3.) GO TO 101
      E=(SUMYY-SUMXY*A          )/(COUNT-2.0)
      GO TO 100
C
C     ISUFFICIENT POINTS
   25 A=0.0
      B=0.0
  101 E=0.0
  100 RETURN
      END
#endif
