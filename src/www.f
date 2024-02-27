C      #######################################################
C      #              HADRON TENSOR (version2.0)             #
C      #                  written by T.Abe                   #
C      #                    (12/03/2000)                     #
C      #######################################################
      double precision function  WWW(mu,nu, Q2,W, p,q, IHD, Mn)
C   mu,nu : Lorentz subscripts
C   Q2    : Momentum transfer squared to the hadronic system in GeV^2 (>0)
C   W     : Mass of the final hadronic system in GeV
C   p(4)  : 4-momentum of incoming nucleon
C   q(4)  : 4-momentum of virtual photon coming into the hadron system
C   IHD   : =+-1(HYDROGEN), =0(DEUTERIUM)
C   Mn    : Mass of incoming nucleon
      implicit NONE
      include './inc/graepia.h'
*---------- Argument ----------
      integer            mu,nu, IHD
      double precision   Q2,W, p(4),q(4), Mn
*------------------------------
*----------- COMMON -----------
      double precision  Dgmunu(4,4)
      common /TA_GMN/   Dgmunu
*------------------------------
*--- Functions defined below ---
      double precision   gmunu
       external          gmunu
*-------------------------------
*------ Local variables ------
      double precision   V
      integer  icnt, i,j
      save     icnt
      double precision  Q2_prev,W_prev,W1,W2,p_q
      save              Q2_prev,W_prev,W1,W2,p_q
*-----------------------------
*---------- DATA ----------
      data icnt/0/, Q2_prev/-1D0/,W_prev/-1D0/
*--------------------------
      if (icnt .EQ. 0) then
        do 1000 i = 1, 4
        do 1000 j = 1, 4
          Dgmunu(i,j) = gmunu(i,j)
 1000   continue
        icnt = icnt + 1
      endif
      if (Q2.EQ.Q2_prev .AND. W.EQ.W_prev)  GOTO 2000
      Q2_prev = Q2
      W_prev  = W
      if     (Istrf .EQ. 1) then
        if (IHD .NE. 1) then
          write(6,*) '!!!Error in WWW!!!'
          write(6,*) '  ---> Invalid IHD =', IHD, '  for Istrf=1'
          write(6,*) '  ---> Good-bye!'
          STOP
        endif
        if (W .LT. 2) then
          call BRASSE(Q2,W, W1,W2)
        else
          call ALLM(Q2,W,1997, W1,W2)
        endif
      else
        write(6,*) '!!!Error in WWW!!!'
        write(6,*) '  ---> Invalid Istrf =', Istrf
        write(6,*) '  ---> Good-bye!'
        STOP
      endif
      p_q = p(4)*q(4) - p(1)*q(1) - p(2)*q(2) - p(3)*q(3)
 2000 continue
      WWW = W1*( -Dgmunu(mu,nu) - q(mu)*q(nu)/Q2 )
     &    + W2/(Mn*Mn) *( p(mu) + p_q/Q2*q(mu) )
     &                 *( p(nu) + p_q/Q2*q(nu) )
      return
      end
C -----------------------------------------------------------------
      double precision function  gmunu(mu,nu)
      implicit NONE
      integer  mu,nu
      if ((mu.EQ.4).AND.(nu.EQ.4)) then
          gmunu =  1.D0
       elseif (mu.EQ.nu) then
          gmunu = -1.D0
       else
          gmunu =  0.D0
      endif
      return
      end
