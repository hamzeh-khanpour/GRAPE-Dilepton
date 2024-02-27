***************************************************
*             Interface to SOPHIA_2.0             *
*              Written by Tetsuo ABE              *
*                 on 12/03/2000                   *
***************************************************
* Ref: Comput. Phys. Commun., 124 (2000), 290-314
***************************************************
      subroutine  Run_SOPHIA(W, PX, Pin, amp)
      implicit NONE
*--------------- Argument ---------------
      double precision  W      ! CM energy of gamma-p system
     &                 ,PX(4)  ! 4-momentum of gamma-p system in LAB frame
     &                 ,Pin(4) ! 4-momentum of incomming proton in LAB
     &                 ,amp    ! Mass of proton
*----------------------------------------
*--------------- PYTHIA stuff ---------------
      integer           Npy, NPAD, Kpy(4000,5)
      double precision  Ppy(4000,5), Vpy(4000,5)
       common /PYJETS/Npy, NPAD, Kpy, Ppy, Vpy
      integer  Npy_X      ! Position of X in /PYJETS/
      integer           MSTU(200),MSTJ(200)
      double precision  PARU(200),PARJ(200)
       COMMON/PYDAT1/MSTU,PARU,MSTJ,PARJ
*--------------------------------------------
*--------------- SOPHIA stuff ---------------
      double precision P
      integer  LLIST, NP, Ideb
       COMMON /S_PLIST/ P(2000,5), LLIST(2000), NP, Ideb
*--------------------------------------------
*-------------- KINEM_4f stuff --------------
      double precision
     &      ttt1,uuu1,ttt2,uuu2,tttt,uuuu,sss2,sss1,q56,q12,q22
     &     ,d1,g1,t1,u1,t2,u2
      common/atest/
     &      ttt1,uuu1,ttt2,uuu2,tttt,uuuu,sss2,sss1,q56,q12,q22
     &     ,d1,g1,t1,u1,t2,u2
*--------------------------------------------
*---------------- Function ----------------
       integer   KF_from_SOPHIA
       external  KF_from_SOPHIA
*------------------------------------------
*------------ Local variables ------------
       integer  L0,Imode, i,j, icount
       save L0,icount
*-----------------------------------------
*---------- DATA ----------
       data icount/0/
*--------------------------
      if (icount .EQ. 0) then
         L0 = 13   ! 13:proton, 14:neutron
         call INITIAL(L0)
         icount = icount + 1
      endif
*(Generation of exclusive hadronic final state)
      call EventGen(L0,W,PX,Pin,amp,ttt1, Imode)
      if (NP .EQ. 0) then
        write(6,*) '!!!Warning in Run_SOPHIA!!!'
        write(6,*) '  ---> # of generated particles =', NP
        write(6,*) '  ---> Something is wrong!'
        RETURN
      endif
      Npy_X = 0
      do i=1,Npy   ! finding X in /PYJETS/
        if (
     +       ( Ppy(i,5) .GT. 1 )   
     + .and. ( Kpy(i,1) .EQ. 1 )   
     + .and. (      ( abs(Kpy(Kpy(i,3), 2)) .EQ. 2112 )     
     +         .or. ( abs(Kpy(Kpy(i,3), 2)) .EQ. 2212 )  )  
     +     ) then
          Npy_X = i
          GOTO 100
        endif
      enddo
 100  continue
      if (Npy_X .EQ. 0) then
        write(6,*) '!!!Warning in Run_SOPHIA!!!'
        write(6,*) '  ---> Cannot find particle:X in /PYJETS/'
        write(6,*) '  ---> Something is wrong!'
        RETURN
      endif
      Kpy(Npy_X, 1) = 11      
      Kpy(Npy_X, 4) = Npy+1   
      do i=1,NP   !!!!!!!! << Loop for decay products from X >> !!!!!!!!
        Npy = Npy + 1    
        do j=1,5
          Ppy(Npy,j) = P(i,j)   
        enddo
        Kpy(Npy,1) = 1       
        Kpy(Npy,2) = KF_from_SOPHIA(LLIST(i))    ! KF code
        Kpy(Npy,3) = Npy_X   ! Mother
      enddo   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Kpy(Npy_X, 5) = Npy     
      do i = 1, min( max(MSTU(70),1), 10 )
        if ( MSTU(70+i) .EQ. (Npy-NP) ) then
          MSTU(70+i) = Npy      
          MSTU(70)   = i        
          GOTO 200
        endif
      enddo
 200  continue
      return
      end
