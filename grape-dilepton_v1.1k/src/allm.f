**************************************************************
*   << Calculation of W1,W2 using ALLM parameterization >>   *
*                   Written by Tetsuo ABE                    *
*                       on 12/03/2000                        *
**************************************************************
***************************************************
* Input:  Q2      in GeV^2        (real*8)
*         W       in GeV          (real*8)
*         Year                    (integer*4)
* Output: W1,W2                   (real*8)
***************************************************
      subroutine ALLM( Q2, W, Year    ! Input
     &                ,W_1, W_2       ! Output
     &               )
      implicit NONE
*-------- Argument --------
      double precision  Q2, W
      integer           Year
      double precision  W_1,W_2
*--------------------------
*-------- Function --------
      double precision  R_LT
      external          R_LT
*--------------------------
*------ Local variables ------
      integer  icount
      save     icount
      double precision  numerator, denominator, one_x, W2,nu
*(ALLM parameter values)
      real  m0m0,mPmP,mRmR, Q0Q0,Lambda2
     &     ,aP1,aP2,aP3, bP1,bP2,bP3, cP1,cP2,cP3
     &     ,aR1,aR2,aR3, bR1,bR2,bR3, cR1,cR2,cR3
      save  m0m0,mPmP,mRmR, Q0Q0,Lambda2
     &     ,aP1,aP2,aP3, bP1,bP2,bP3, cP1,cP2,cP3
     &     ,aR1,aR2,aR3, bR1,bR2,bR3, cR1,cR2,cR3
*(ALLM parameters)
      double precision  t, x,xP,xR, aP,bP,cP, aR,bR,cR, F2P,F2R,F2
     &     ,ln_Q0Q0_Lambda2
      save  ln_Q0Q0_Lambda2
*(Mathmatical constant)
      double precision  pi
      save              pi
*(Physics constants)
      double precision  Mp,Mp2,DMp, R
      save              Mp,Mp2,DMp, R
*-----------------------------
*---------- DATA ----------
      data  icount/0/
*--------------------------
      if     ( Q2 .LT. 0 ) then
        write(6,*) '!!!Error in ALLM_TA!!!'
        write(6,*) '  ---> Q2(=', real(Q2),') should be positive.'
        write(6,*) '  ---> Good-bye!'
        STOP
      elseif ( Q2 .GT. 5000 ) then
C        write(6,*) '!!!Warning in ALLM_TA!!!'
C        write(6,*) '  ---> Q2(=', real(Q2),') is outside '
C        write(6,*) '       the fitting region.'
      endif
      W2 = W*W
      if     ( W2 .LE. Mp2 ) then
        write(6,*) '!!!Error in ALLM_TA!!!'
        write(6,*) '  ---> W2(=', real(W2),') should be >Mp^2GeV^2.'
        write(6,*) '  ---> Good-bye!'
        STOP
      elseif ( W2 .LT. 2.9999 ) then
        write(6,*) '!!!Warning in ALLM_TA!!!'
        write(6,*) '  ---> W2(=', real(W2),') is ouside '
        write(6,*) '       the fitting region.'
        write(6,*) '  ---> Good-bye!'
      endif
      if (icount .EQ. 0) then
C        write(6,*) '!!!!!!Message from ALLM_TA!!!!!!'
C (Ref: DESY-report97-251, hep-ph/9712415)
        if     (Year .EQ. 1991) then
          m0m0 = 0.30508
          mPmP = 10.676
          mRmR = 0.20623
          Lambda2 = 0.06527
          Q0Q0    = 0.27799 + Lambda2
          cP1 = 0.26550
          cP2 = 0.04856
          cP3 = 1.04682
          aP1 = -0.04503
          aP2 = -0.36407
          aP3 = 8.17091
          bP1 = 0.49222**2
          bP2 = 0.52116**2
          bP3 = 3.5515
          cR1 = 0.67639
          cR2 = 0.49027
          cR3 = 2.66275
          aR1 = 0.60408
          aR2 = 0.17353
          aR3 = 1.61812
          bR1 = 1.26066**2
          bR2 = 1.83624**2
          bR3 = 0.81141
C          write(6,*) '  ---> ALLM91 parameter values are used.'
        elseif (Year .EQ. 1997) then
          m0m0 = 0.31985
          mPmP = 49.457
          mRmR = 0.15052
          Lambda2 = 0.06527
          Q0Q0    = 0.46017 + Lambda2
          cP1 = 0.28067
          cP2 = 0.22291
          cP3 = 2.1979
          aP1 = -0.0808
          aP2 = -0.44812
          aP3 = 1.1709
          bP1 = 0.60243**2
          bP2 = 1.3754**2
          bP3 = 1.8439
          cR1 = 0.80107
          cR2 = 0.97307
          cR3 = 3.4942
          aR1 = 0.58400
          aR2 = 0.37888
          aR3 = 2.6063
          bR1 = 0.10711**2
          bR2 = 1.9386**2
          bR3 = 0.49338
C          write(6,*) '  ---> ALLM97 parameter values are used.'
        else
          write(6,*) '!!!Error in ALLM_TA!!!'
          write(6,*) '  ---> You specified unknown Year(=',Year,').'
          write(6,*) '  ---> Good-bye!'
          STOP
        endif
        pi = acos(-1D0)
        ln_Q0Q0_Lambda2 = log( Q0Q0/Lambda2 )
        Mp    = 0.93827231D0
        Mp2   = 0.93827231D0**2
        DMp   = 0.93827231D0*2D0
        R = 0.18
C        write(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        icount = icount + 1
      endif
      if (Q2 .GT. 0) then
        numerator = log( (Q2+Q0Q0)/Lambda2 )
        if (numerator .GT. ln_Q0Q0_Lambda2) then
          t = log( numerator/ln_Q0Q0_Lambda2 )
        else
          t = 0D0
        endif
      else
         t = 0D0
      endif
C      denominator = 1D0 + (W2-Mp2)/Q2
C       x = 1D0/denominator
      denominator = Q2 + (W2-Mp2)
       x = Q2/denominator
      if     ( x.LT.0 .OR. x.GE.1 ) then
        write(6,*) '!!!Error in ALLM_TA!!!'
        write(6,*) '  ---> x(=', real(x),') should be in [0,1).'
        write(6,*) '  ---> Good-bye!'
        STOP
C      elseif ( x.GT.0.85 ) then
C        write(6,*) '!!!Warning in ALLM_TA!!!'
C        write(6,*) '  ---> x(=', real(x),') is outside '
C        write(6,*) '       the fitting region.'
      endif
      R = R_LT(x,Q2)
C      denominator = 1D0 + (W2-Mp2)/(Q2+mPmP)
C       xP = 1D0/denominator
C      denominator = 1D0 + (W2-Mp2)/(Q2+mRmR)
C       xR = 1D0/denominator
      denominator = (Q2+mPmP) + (W2-Mp2)
       xP = (Q2+mPmP)/denominator
      denominator = (Q2+mRmR) + (W2-Mp2)
       xR = (Q2+mRmR)/denominator
      aP = aP1 + (aP1-aP2)*(1D0/(1D0+t**aP3)-1D0)
      bP = bP1 + bP2*t**bP3
      cP = cP1 + (cP1-cP2)*(1D0/(1D0+t**cP3)-1D0)
      aR = aR1 + aR2*t**aR3
      bR = bR1 + bR2*t**bR3
      cR = cR1 + cR2*t**cR3
      one_x = 1D0-x
      F2P = cP *xP**aP *(one_x)**bP
      F2R = cR *xR**aR *(one_x)**bR
      F2 = Q2 /(Q2+m0m0) *(F2P+F2R)
      nu = (W2+Q2-Mp2) /DMp
      W_2 = F2 /nu  !*Mp
C      W_1 = W_2 *(Q2+nu*nu) /Q2 /(1D0+R)
      W_1 = 1D0 /(Q2+m0m0) *(F2P+F2R)   ! = F2/Q2
     &      /nu *Mp
     &      *(Q2+nu*nu) /(1D0+R)
      return
      end
