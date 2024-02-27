************************************************************************
*980110
************************************************************************
      SUBROUTINE SMINTV(LP, AM, PI, EP, EW, VM, IGAUG)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER LP
      DIMENSION PI(4), EP(4,LP), EW(LP)
      common /chcntl/ jwidth, jtgamm
*
*     print *,'smintv ... called'
      if(jtgamm .eq. 0 .or. vm .gt. 0.0d0 .or. am .gt. 0.0d0 ) then
          call sminv0(lp,am,pi,ep,ew,vm,igaug)
*         print *,'inv0'
      else
          call smintp(lp,am,pi,ep,ew,vm,igaug)
*         print *,'intp'
      endif

      return
      end
