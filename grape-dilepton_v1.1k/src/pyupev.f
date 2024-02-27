***************************************************
*          User-defined external process          *
*                 written by T.Abe                *
*                on Aug. 06 in 1998               *
***************************************************
      subroutine PYUPev(ISUB,SIGEV)
      implicit NONE
* ---------- Argument ----------
      integer*4           ISUB
      double precision  SIGEV
* ------------------------------
      include './inc/py_common.h'
      include './inc/graepia.h'
      integer           NUP,KUP(20,7),NFUP,IFUP(10,2)     
      double precision  PUP(20,5),Q2UP(0:10)              
       COMMON /PYUPPR/ NUP,KUP,NFUP,IFUP,PUP,Q2UP
                                                   
      double precision  UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GL
       common /GEP_PDF/ UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GL
      double precision    Q2_QCD
       common /QCD_SCALE/ Q2_QCD
* ---------------- GRACE stuff ----------------
      integer     nextrn
       parameter (nextrn=6)
      double precision  vec(4,nextrn)
       common /sp4vec/  vec
      integer     mxextn
       parameter (mxextn = 10)
      double precision  amass1(mxextn), amass2(mxextn)
       common /kmmass/  amass1,         amass2
      integer           kcharg(mxextn), kfcode(mxextn)
       common /kminfo/  kcharg,         kfcode
      double precision  S(3),W(3),FACT
       COMMON /KINEM1/  S,   W,   FACT
      double precision  func
       external         func
      integer          jproc      
       COMMON /amjprc/ jproc
* ---------------------------------------------
* ------------ BASES common on its result ------------
      integer           ITG,ITF
      real*4            STIME
      double precision  AVGI,SD,CHI2A
      COMMON /BSRSLT/AVGI,SD,CHI2A,STIME,ITG,ITF
* ----------------------------------------------------
* ---------------- Cross-sections ----------------
C      double precision  x_sec(500), x_sec_err(500)
C       common /X_SEC/   x_sec,      x_sec_err
* ------------------------------------------------
* ----------- Kinematical variables -----------
      double precision  P1_lab,P2_lab, E1_lab,E2_lab  
     &                 ,Pcms_lab(3),Ecms_lab(3)
     &                 ,GAMMAcms_lab(3),BETGAMcms_lab(3)
     &                 ,vec_isr(4)
       common /GEP_LAB/ P1_lab,P2_lab, E1_lab,E2_lab
     &                 ,Pcms_lab,Ecms_lab
     &                 ,GAMMAcms_lab,   BETGAMcms_lab
     &                 ,vec_isr
* ---------------------------------------------
* ---------- Functions ----------
      integer    KFencode
      double precision  PYR
       external  KFencode, PYR
* -------------------------------
* ---------- Local variables ----------
      integer*4   Status_Code(10)
      integer*4            i,ii, my_code
      double precision   Pz_lab,E_lab
     &                  ,vec_cms(4,nextrn)
     &                  ,prob_PDF(12)
     &                  ,sum
C      real*4   choice
      double precision  choice
* ---------------- DATA ----------------
      data  Status_Code/2,2,1,1,1,1,1,1,1,1/
      data  prob_PDF/12*0D0/
* --------------------------------------
      call spring( func, mxtry )      !!! Event generation by SPRING !!!
* ------ Recovering PE in LAB. ------
      if (Frame_amp .EQ. 2) then
         do ii=1,NEXTRN
           do i=1,4
             vec(i,ii) = PE_str(i,ii)
           enddo
         enddo
      endif
* ------ Boost : Lab. frame ---> CMS of incoming beams ------
      do NUP=1,NEXTRN
        Pz_lab = vec(3,NUP)
        E_lab  = vec(4,NUP)
          vec_cms(1,NUP) = vec(1,NUP)
          vec_cms(2,NUP) = vec(2,NUP)
          vec_cms(3,NUP) = GAMMAcms_lab(1)  * Pz_lab
     &                   - BETGAMcms_lab(1) * E_lab
          vec_cms(4,NUP) = GAMMAcms_lab(1)  * E_lab
     &                   - BETGAMcms_lab(1) * Pz_lab
      enddo
* ------------ Filling event infomation in /PYUPPR/ ------------
      NUP = 0
      if     ((process.EQ.1).or.(process.EQ.2)) then   
         NUP = 6
         do i = 1, NUP
           KUP(i,4) = 0
           KUP(i,5) = 0
           KUP(i,6) = 0
           KUP(i,7) = 0
         enddo
       elseif (process.EQ.3) then      
         NUP = NUP + 1   
         KUP(NUP,4) = 0
         KUP(NUP,5) = 0
         KUP(NUP,6) = 3
         KUP(NUP,7) = 3
          NUP = NUP + 1   
          KUP(NUP,4) = 0
          KUP(NUP,5) = 0
          KUP(NUP,6) = 0
          KUP(NUP,7) = 0
           NUP = NUP + 1   
           KUP(NUP,4) = 1
           KUP(NUP,5) = 1
           KUP(NUP,6) = 0
           KUP(NUP,7) = 0
            NUP = NUP + 1   
            KUP(NUP,4) = 0
            KUP(NUP,5) = 0
            KUP(NUP,6) = 0
            KUP(NUP,7) = 0
             NUP = NUP + 1   
             KUP(NUP,4) = 0
             KUP(NUP,5) = 0
             KUP(NUP,6) = 0
             KUP(NUP,7) = 0
              NUP = NUP + 1   
              KUP(NUP,4) = 0
              KUP(NUP,5) = 0
              KUP(NUP,6) = 0
              KUP(NUP,7) = 0
       elseif (process.EQ.4) then      
         NUP = NUP + 1   
         KUP(NUP,4) = 0
         KUP(NUP,5) = 0
         KUP(NUP,6) = 2
         KUP(NUP,7) = 2
          NUP = NUP + 1   
          KUP(NUP,4) = 0
          KUP(NUP,5) = 0
          KUP(NUP,6) = 1
          KUP(NUP,7) = 1
           NUP = NUP + 1   
           KUP(NUP,4) = 4
           KUP(NUP,5) = 4
           KUP(NUP,6) = 0
           KUP(NUP,7) = 0
            NUP = NUP + 1   
            KUP(NUP,4) = 3
            KUP(NUP,5) = 3
            KUP(NUP,6) = 0
            KUP(NUP,7) = 0
       elseif (process.EQ.5) then      
         NUP = NUP + 1   
         KUP(NUP,4) = 0
         KUP(NUP,5) = 0
         KUP(NUP,6) = 3
         KUP(NUP,7) = 3
          NUP = NUP + 1   
          KUP(NUP,4) = 0
          KUP(NUP,5) = 0
          KUP(NUP,6) = 0
          KUP(NUP,7) = 0
           NUP = NUP + 1   
           KUP(NUP,4) = 1
           KUP(NUP,5) = 1
           KUP(NUP,6) = 0
           KUP(NUP,7) = 0
            NUP = NUP + 1   
            KUP(NUP,4) = 0
            KUP(NUP,5) = 0
            KUP(NUP,6) = 0
            KUP(NUP,7) = 0
             NUP = NUP + 1   
             KUP(NUP,4) = 0
             KUP(NUP,5) = 0
             KUP(NUP,6) = 0
             KUP(NUP,7) = 0
       else
         write(6,*) '!!!Error in pyupev.f!!!'
         write(6,*) ' ---> unknown process(process=',process,')'
         write(6,*) ' ---> Good-bye!'
         STOP
      endif
      do i = 1,NUP
        KUP(i,1) = Status_Code(i)
        KUP(i,2) = KFcode(i)
        KUP(i,3) = 0      
        PUP(i,1) = vec_cms(1,i)
        PUP(i,2) = vec_cms(2,i)
        PUP(i,3) = vec_cms(3,i)
        PUP(i,4) = vec_cms(4,i)
        PUP(i,5) = amass1(i)
      enddo
* ---------------- Merge-mode in DIS ----------------
      if ((process.EQ.3).and.(merge.NE.0)) then   
         prob_PDF( 1) = UPV + USEA                
         prob_PDF( 2) = USEA                      
         prob_PDF( 3) = DNV/4.D0 + DSEA/4.D0      
         prob_PDF( 4) = DSEA/4.D0                 
         prob_PDF( 5) = STR/4.D0                  
         prob_PDF( 6) = STR/4.D0                  
         prob_PDF( 7) = CHM                       
         prob_PDF( 8) = CHM                       
         prob_PDF( 9) = BOT/4.D0                  
         prob_PDF(10) = BOT/4.D0                  
         prob_PDF(11) = TOP                       
         prob_PDF(12) = TOP                       
         if (qflv .GE. 2) then
           do i = 1, qflv-1
             prob_PDF(i) = 0.D0
           enddo
         endif
         sum = 0.D0
         do i = qflv, 12
           sum = sum + prob_PDF(i)
           if (i.GE.2)  prob_PDF(i) = ( prob_PDF(i) + prob_PDF(i-1) )
         enddo
         do i = qflv, 12
           prob_PDF(i) = prob_PDF(i) / sum
         enddo
 111     continue
C         call RANLUX(choice,1)
         choice = PYR(0)
         do i = qflv, 12
           if (prob_PDF(i) .GT. choice) then   
             my_code = i
             KUP(1,2) = KFencode(i)      
             KUP(3,2) = KUP(1,2)         
              kfcode(1) = KUP(1,2)        
              kfcode(3) = KUP(1,2)        
             GOTO 222
           endif
         enddo
 222     continue
        if ((qflv.EQ. 1).and.(merge .EQ. 1234)
     +                       .and.(jgra_sel .LE. 2)) then    
          if (  .NOT.( (my_code.EQ.1).or.(my_code.EQ.2)
     +             .or.(my_code.EQ.3).or.(my_code.EQ.4) ) )  GOTO 111
         elseif ((qflv.EQ. 1).and.(merge .EQ. 123456)
     +                       .and.(jgra_sel .LE. 2)) then    
          if (  .NOT.( (my_code.EQ.1).or.(my_code.EQ.2)
     +             .or.(my_code.EQ.3).or.(my_code.EQ.4)
     +             .or.(my_code.EQ.5).or.(my_code.EQ.6) ) )  GOTO 111
         elseif ((qflv.EQ. 1).and.(merge .EQ. 12345678)
     +                       .and.(jgra_sel .LE. 2)) then    
          if (  .NOT.( (my_code.EQ.1).or.(my_code.EQ.2)
     +             .or.(my_code.EQ.3).or.(my_code.EQ.4)
     +             .or.(my_code.EQ.5).or.(my_code.EQ.6)
     +             .or.(my_code.EQ.7).or.(my_code.EQ.8) ) )  GOTO 111
         elseif ((qflv.EQ. 1).and.(merge .EQ. 1234567890)
     +                       .and.(jgra_sel .LE. 2)) then    
          if (  .NOT.( (my_code.EQ.1).or.(my_code.EQ.2)
     +             .or.(my_code.EQ.3).or.(my_code.EQ.4)
     +             .or.(my_code.EQ.5).or.(my_code.EQ.6)
     +             .or.(my_code.EQ.7).or.(my_code.EQ.8)
     +             .or.(my_code.EQ.9).or.(my_code.EQ.10)) )  GOTO 111
         elseif ((qflv.EQ. 1).and.(merge .EQ. 17)
     +                       .and.(jgra_sel .GE. 1)) then    
          if (  .NOT.( (my_code.EQ.1).or.(my_code.EQ.7) ) )  GOTO 111
         elseif ((qflv.EQ. 2).and.(merge .EQ. 28)
     +                       .and.(jgra_sel .GE. 1)) then    
          if (  .NOT.( (my_code.EQ.2).or.(my_code.EQ.8) ) )  GOTO 111
         elseif ((qflv.EQ. 3).and.(merge .EQ. 35)
     +                       .and.(jgra_sel .GE. 1)) then    
          if (  .NOT.( (my_code.EQ.3).or.(my_code.EQ.5) ) )  GOTO 111
         elseif ((qflv.EQ. 3).and.(merge .EQ. 359)
     +                       .and.(jgra_sel .GE. 1)) then    
          if (  .NOT.( (my_code.EQ.3).or.(my_code.EQ.5)
     +                               .or.(my_code.EQ.9) ) )  GOTO 111
         elseif ((qflv.EQ. 4).and.(merge .EQ. 46)
     +                       .and.(jgra_sel .GE. 1)) then    
          if (  .NOT.( (my_code.EQ.4).or.(my_code.EQ.6) ) )  GOTO 111
         elseif ((qflv.EQ. 4).and.(merge .EQ. 460)
     +                       .and.(jgra_sel .GE. 1)) then    
          if (  .NOT.( (my_code.EQ.4).or.(my_code.EQ.6)
     +                               .or.(my_code.EQ.10)) )  GOTO 111
         else
           write(6,*) '!!!Error in PYUPev!!!'
           write(6,*) '  ---> Not supported qflv/merge(='
     &                   ,qflv,'/',merge,' )'
           write(6,*) '  ---> Good-bye!'
           STOP
        endif
      endif   
* ------------ Parton Shower Parameters ------------
         Q2UP(0) = Q2_QCD   
         Q2UP(1) = S(3)
         Q2UP(2) = S(3)
         if (MSTP(71) .NE. 0) then   
C           call RANLUX(choice,1)
           choice = PYR(0)
           if ((process.EQ.1).or.(process.EQ.2)) then   
C             NFUP =  2      
             NFUP =  1      
             if (choice.LT.1./3.) then
                 IFUP(1,1) = 5   
                 IFUP(1,2) = 6   
C                 IFUP(2,1) = 3   
C                 IFUP(2,2) = 4   
              elseif (choice.LT.2./3.) then
                 IFUP(1,1) = 4
                 IFUP(1,2) = 6
C                 IFUP(2,1) = 3
C                 IFUP(2,2) = 5
              else
                 IFUP(1,1) = 4
                 IFUP(1,2) = 5
C                 IFUP(2,1) = 3
C                 IFUP(2,2) = 6
             endif
           elseif (process .EQ. 3) then   
             NFUP =  2      
             Q2UP(1) = Q2_QCD*4D0
             if (choice.LT.1./3.) then
                 IFUP(1,1) = 3   
                 IFUP(1,2) = 4   
                 IFUP(2,1) = 5   
                 IFUP(2,2) = 6   
              elseif (choice.LT.2./3.) then
                 IFUP(1,1) = 3
                 IFUP(1,2) = 5
                 IFUP(2,1) = 4
                 IFUP(2,2) = 6
              else
                 IFUP(1,1) = 3
                 IFUP(1,2) = 6
                 IFUP(2,1) = 4
                 IFUP(2,2) = 5
             endif
           elseif (process .EQ. 4) then   
             NFUP =  1
                 IFUP(1,1) = 3
                 IFUP(1,2) = 4
           elseif (process .EQ. 5) then   
             NFUP =  1
             if (choice.LT.1./2.) then
                 IFUP(1,1) = 3
                 IFUP(1,2) = 4
              else
                 IFUP(1,1) = 3
                 IFUP(1,2) = 5
             endif
           else
             NFUP =  0
           endif   
         endif   
* ---------- Parameters for the simulation of hadronic system ----------
      if     ((process .EQ. 1).or.(process .EQ. 2)) then
         MINT(45) = 1
           MINT(46) = 1
           MINT(47) = 1
         PUP(2,4) = PUP(1,4)   
         MSTP(61) = 0          
       elseif (process .EQ. 3) then
         MINT(45) = 2
         if     (Isr_flag .EQ. 0) then
           MINT(46) = 1
           MINT(47) = 3
         elseif (Isr_flag .EQ. 1) then
           MINT(46) = 3
           MINT(47) = 4
         elseif (Isr_flag .EQ. 2) then
           write(6,*) '!!!Error in pyupevt!!!'
           write(6,*) '  ---> Isr_flag =',Isr_flag
     &                        ,' is not yet supported.'
         else
           write(6,*) '!!!Error in pyupevt!!!'
           write(6,*) '  ---> Isr_flag =',Isr_flag
     &                        ,' is not yet supported.'
         endif
       elseif (process .EQ. 4) then
         MINT(45) = 2
         MINT(46) = 2
         MINT(47) = 4
       elseif (process .EQ. 5) then
         MINT(45) = 2
         MINT(46) = 1
         MINT(47) = 3
       else
         write(6,*) '!!!Error in pyupev.f!!!'
         write(6,*) ' ---> unknown process(=',process,')'
         write(6,*) ' ---> Good-bye!'
         STOP
      endif
* ----------------------------------------------------------------------
* --------------------------------------------------------------
C      SIGEV = 100.D0      
C      SIGEV = x_sec(ISUB)      
      SIGEV = AVGI      
* --------------------------------------------------------------
      return
      end
