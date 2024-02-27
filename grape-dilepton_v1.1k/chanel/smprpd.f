************************************************************************
      SUBROUTINE SMPRPD(APROP, AMOMQ, AMASSQ, AMAG)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 APROP
      REAL*8     AMOMQ, AMASSQ, AMAG
*
*    Calculate denominator of propagator.
*
*     APROP  : in/out : product of denominator of propagators.
*     AMOMP  : input  : square of mementum.
*     AMASSQ : input  : square of mass.
*     AMAG   : input  : mass * width.
*-----------------------------------------------------------------------
      IF (AMOMQ .GT. 0) THEN
        APROP = - APROP*DCMPLX(AMOMQ - AMASSQ, AMAG)
      ELSE
        APROP = - APROP*DCMPLX(AMOMQ - AMASSQ, 0.0D0)
      ENDIF
 
*     CALL CTIME('SMPRPD')
      RETURN
      END
