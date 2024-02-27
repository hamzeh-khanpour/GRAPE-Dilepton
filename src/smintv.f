************************************************************************
*980110
************************************************************************
      SUBROUTINE SMINTV(LP, AM, PI, EP, EW, VM, IGAUG)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER LP
      DIMENSION PI(4), EP(4,LP), EW(LP)
      common /chcntl/ jwidth, jtgamm
*(TA)
      external Alpha_Run
      common /smcnst/ TApi,TApi2,TArad,TAgevpb,TAalpha,TAalphas,alpha0
      common/CMN_ISRPOL/isrpol

*
*     print *,'smintv ... called'
      if(jtgamm .eq. 0 .or. vm .gt. 0.0d0 .or. am .gt. 0.0d0 ) then
          call sminv0(lp,am,pi,ep,ew,vm,igaug)
*         print *,'inv0'
      else
          call smintp(lp,am,pi,ep,ew,vm,igaug)
*         print *,'intp'
      endif

*(TA)
      if (isrpol .GT. 0) then
        if (AM .EQ. 0D0) then   ! This vector boson is photon.
           Alpha_tmp = Alpha_Run(VM,alpha0)  ! TA(24/09/2000)
           do 1000 i = 1, LP
            EW(i) = EW(i) *Alpha_tmp
 1000      continue
        endif
      endif

      return
      end
