      program MAINSP
      implicit real*8(a-h,o-z)
      include 'inclk.h'
      common / knmflg / knmgvs, knmhst, knmsph, knmisr
      external func
      parameter (nextn  =6)
      common /sp4vec/ vec(4,mxextn)
      common /amjprc/ jproc
      double precision  PYmass
       external         PYmass
      include './inc/graepia.h'
      include './inc/py_common.h'
* ---------- Kinematic variables ----------
      double precision  P1_lab,P2_lab, E1_lab,E2_lab
     &                 ,Pcms_lab(3),Ecms_lab(3)
     &                 ,GAMMAcms_lab(3),BETGAMcms_lab(3)
     &                 ,vec_isr(4)
       common /GEP_LAB/ P1_lab,P2_lab, E1_lab,E2_lab
     &                 ,Pcms_lab,Ecms_lab
     &                 ,GAMMAcms_lab,   BETGAMcms_lab
     &                 ,vec_isr
* -----------------------------------------
* ------------ BASES common on its result ------------
      integer           ITG,ITF
      real*4            STIME
      double precision  AVGI,SD,CHI2A
      COMMON /BSRSLT/AVGI,SD,CHI2A,STIME,ITG,ITF
* ----------------------------------------------------
* ------------ PYTHIA stuff ------------
      integer          MDCY(500,3)
       common /PYDAT3/ MDCY
      integer     PYCOMP
       external   PYCOMP
* --------------------------------------
* ------------ USER common ------------
      integer          ISPRING
       common /TA_USR/ ISPRING
* -------------------------------------
* -------- Local variables for quasi-elastic process --------
      double precision  P_tmp(4), P_in(4)
* -----------------------------------------------------------
* ---------- Local variables ----------
      integer   LUN_cards,    LUN_rnd,    LUN_nt,    LUN_asc
      parameter(LUN_cards=30, LUN_rnd=31, LUN_nt=32, LUN_asc=33)
      integer  IKF,NKF,KF
       parameter  (NKF=19)
      integer       KFSD(NKF)
      character*16  NASD(NKF)
* -------------------------------------
* ---------------- DATA ----------------
      data KFSD / 310,  221, 3122, 3222, 3212
     &           ,3112, 3322, 3312, 3334
     &           ,411,   421,  431, 4122
     &           ,13,  15,   211, 321, 130, 111 /
* --------------------------------------
      ISPRING = 1
      call START_GRAPE(2)
      print_flag = 1
* ------ Initialization of FFREAD and reading control_cards ------
      call Read_Cards(LUN_cards,  'grape.cards')
      call Write_cards
* ----------------------------------------------------------------
* ------ Getting process ID: /amjprc/jproc ------
      call Get_Proc
* -----------------------------------------------
* ------ Getting flags for Feynman Graph Selection ------
      call Get_Graph_flag
* -------------------------------------------------------
* --------- Suppressing decays of particles ---------
      do IKF=1,NKF
        KF = KFSD(IKF)
        MDCY(PYCOMP(KF), 1) = 0
        call PYNAME(KF,NASD(IKF))
      enddo
C      write(6,*) ' '
      write(6,*)
     &  'Suppressing decay of the following particles in PYTHIA,'
      write(6,12) (NASD(IKF),IKF=1,NKF)
 12   format(2X,7(1X,A10),/,2X,7(1X,A10),/,2X,7(1X,A10))
C      write(6,*) ' '
* ---------------------------------------------------
*-----------------------------------------------------------------------
      knmgvs = 1
      knmhst = 0
      knmsph = 0
      knmisr = 0
************************************************************************
*               initialization of BASES/SPRING 5.1
************************************************************************
*=======================================================================
*          initialization of bases by calling bsinit
*=======================================================================
*         -------------
           call bsinit
*         -------------
*=======================================================================
*          read the probability information from the file
*=======================================================================
           lun = 23
C           open(lun,file='bases.data',status='old',form='unformatted'
C     &          ,err=80)
C           GOTO 90
C 80        write(6,*) '!!!Error in MainSP!!!'
C           write(6,*) '  ---> Input file(bases.data) does not exis!'
C           write(6,*) '  ---> Good-bye!'
C           STOP
*         -------------
 90        call bsread( lun )
*         -------------
C           close( lun )
*=======================================================================
*      initialization of parameters
*          for kinematics and matrix elements
*      initialization of histograms
*=======================================================================
      lu     = 12
      open(lu,file='spring.result',status='unknown',form='formatted')
*         ------------------
           call userin( lu )   
*         ------------------   
* ---------- Initialization of HBOOK ----------
      call HB_init(LUN_nt, NTPYT_flag, NTVEC_flag)
* ---------------------------------------------
* --- Definition of User-defined External Processes ---
* --- Initialization of PYTHIA Generation Procedure ---
      call PY_Proc_Def      
* -----------------------------------------------------
* ------ Initialization of RANDOM number generation ------
      write(6,*) ' '
      call RND_Init(LUN_rnd, LRND_flag)
* --------------------------------------------------------
* ---------- USER Initialization ----------
         call USRSTR(0, Ngen)
* -----------------------------------------
*-----------------------------------------------------------------------
      call UXdate(Iyear, Imonth, Iday, Ihour, Iminutes)
      write(6,*) ' '
      write(6,*) ' '
      write(6,107)  Ihour, Iminutes, Iday, Imonth, mod(Iyear,100)
*-----------------------------------------------------------------------
*=======================================================================
*     Event Generation
*=======================================================================
      mxevnt = Ngen
      write(6,*) ' '
      write(6,*)  'Number of generated events =', mxevnt
      print_flag = 0      
      do 100 nevnt = 1, mxevnt   
         call PYevnt
         if ((process.EQ.1).or.(process.EQ.2)) then
           do i = 1, 5
             P(5,i) = P(3,i)
             if (i .LE. 4) then
                P(6,i) = P(4,i) - vec_isr(i)
              else
                P(6,i) = P(4,i)
             endif
           enddo
         endif
* ---------- Storing ISR-photon(s) in /PYJETS/ ----------
         if (   ( Isr_flag .EQ. 1 )
     +     .and.( (.not.(MSTP(61).EQ.1).and.(MSTJ(41).EQ.2)) ) 
     +     .and.(vec_isr(4) .GT. E_min_isr)
     +      ) then
           N = N + 1
           do i=1,4
             P(N,i) = vec_isr(i)
           enddo
           P(N,5) = 0.D0    
           K(N,1) = 1       
           K(N,2) = 22      
           K(N,3) = 4       
           do i = 1, min( max(MSTU(70),1), 10 )
             if ( MSTU(70+i) .EQ. (N-1) ) then
               MSTU(70+i) = N        
               MSTU(70)   = i        
               GOTO 111
             endif
           enddo
 111       continue
         endif
* -------------------------------------------------------
* ---------- Hadron generation in Quasi-elastic process ----------
         if (process.EQ.2 .AND. Qela_decay.GT.0) then
           do i=1,N
             if ((abs(K(i,2)).EQ.2212).and.(K(i,1).EQ.1)) then
               do j=1,4
                 P_tmp(j) = P(i,j)   
               enddo
               GOTO 222
             endif
           enddo
 222       continue
           if     (Qela_decay .EQ. 1) then
             do j=1,4
               P_in(j) = P(1,j)
             enddo
             call Run_SOPHIA(P(i,5), P_tmp, P_in, PYmass(2212))
           endif
         endif
* ----------------------------------------------------------------
* ---------- USER Event Storing ----------
         call USRSTR(nevnt, Ngen)
* ----------------------------------------
         if ((LIST_flag).and.(nevnt.LE.Nlist))  call PYlist(1)
         if (Isr_flag.EQ.1) then
            Nisr = 1
          else
            Nisr = 0
         endif
         call FILL_nt(LUN_nt, Ngen, nevnt, merge, nextn, Nisr
     &                                   ,NTPYT_flag,NTVEC_flag)
         if (mod(nevnt,Nmod) .EQ. 0) then
           call UXdate(Iyear, Imonth, Iday, Ihour, Iminutes)
           write(6,117) nevnt,Ihour,Iminutes,Iday,Imonth,mod(Iyear,100)
         endif
  100 continue   
      call UXdate(Iyear, Imonth, Iday, Ihour, Iminutes)
      write(6,127) Ihour, Iminutes, Iday, Imonth, mod(Iyear,100)
* -------- Dumping RANDOM number parameters --------
      write(6,*) ' '
      call RND_Dump(LUN_rnd, LRND_flag)
* --------------------------------------------------
      write(6,*) ' '
      call PYstat(1)
      call spinfo( lu )
      call shplot( lu )
      close( lu )
C      call HB_term
C      close(LUN_asc)
* (Termination in case of nevnt<Ngen)
C      if (nevnt .LT. Ngen) then
C         call FILL_nt(LUN_nt, 0, nevnt, NTPYT_flag, NTVEC_flag)
C      endif
* ----------- USER Termination -----------
         call USRSTR(Ngen+1, Ngen)
* ----------------------------------------
      call Print_x_sec( 1, AVGI, SD )
 107  format('==========> START of SPRING at '
     &        , I2, ':', I2, '(', I2, '/', I2, '/', I2, ')' )
 117  format('   ---> ', I9, ' events have been generated at '
     &        , I2, ':', I2, '(', I2, '/', I2, '/', I2, ').')
 127  format('==========> END   of SPRING at '
     &        , I2, ':', I2, '(', I2, '/', I2, '/', I2, ')' )
      stop
      end
