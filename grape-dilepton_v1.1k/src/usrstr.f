************************************************
*          USER Event Storing Routine          *
*              prepared by T.Abe               *
*             on May. 04 in 1999               *
************************************************
      subroutine USRSTR(Ievt, Ngen)
      implicit NONE
*------------ Arguments ------------
      integer  Ievt, Ngen
* Ngen : # of events to be generated
* Ievt : Counter --- < 1      ===> Initialization   phase
*                    1 - Ngen ===> Event generation phase
*                    > Ngen   ===> Termination      phase
*-----------------------------------
*---------- PYTHIA common ----------
      integer           N,NPAD, K(4000,5)
      double precision  P(4000,5), V(4000,5)
       common /PYJETS/ N,NPAD,K,P,V              !!! Event Record !!!
      integer           MINT(400)
      double precision  VINT(400)
       common /PYINT1/ MINT,VINT
* (See PYTHIA manual for details.)
*-----------------------------------
*--------- Local variables ---------
      integer     LUN1,    LUN2,    LUN3
       parameter (LUN1=41, LUN2=42, LUN3=43)
* You can use the above logical unit numbers.
      integer  ISUB
*-----------------------------------
*-------------- DATA --------------
*     data
*----------------------------------
      ISUB = MINT(1)   !<--- process ID in PYTHIA
******** Initialization of USER Event Storing *******
      if (Ievt .LT. 1) then
      endif
*****************************************************
************* <<< USER Event Storing >>> ************
      if ((Ievt .GE. 1).and.(Ievt .LE. Ngen)) then
      endif
*****************************************************
********* Termination of USER Event Storing *********
      if (Ievt .GT. Ngen) then
      endif
*****************************************************
      return
      end
