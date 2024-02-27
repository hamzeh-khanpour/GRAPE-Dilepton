***************************************************************
*    Definition of User-defined External Process in PYTHIA    *
*      and Initialization of PYTHIA Generation Procedure      *
*                      written by T.Abe                       *
*                     on Aug. 06 in 1998                      *
***************************************************************
      subroutine  PY_Proc_Def
      implicit NONE
* --------- Argument ---------
* ----------------------------
      include './inc/py_common.h'
      include './inc/graepia.h'
* ------------- GRACE common variables -------------
      double precision  amw,     amz,     ama,     amg,     amh
     &                 ,amx,     amy, amnu(3), amlp(3), amuq(3)
     &                 ,amdq(3), amp,     amn,    amcp,    amcm
     &                 ,amcz,   amca,    amcg,  amjp1s
      common /smmass/   amw,     amz,     ama,     amg,     amh
     &                 ,amx,     amy,    amnu,    amlp, amuq
     &                 ,amdq,    amp,     amn,    amcp,    amcm
     &                 ,amcz,   amca,    amcg,  amjp1s
      double precision  agw,      agz,      agh,     agx,      agy
     &                 ,aguq(3),  agdq(3),  agcp,    agcm,     agcz
     &                 ,agjp1s
      common /smgmma/   agw,      agz,      agh,     agx,      agy
     &                 ,aguq,     agdq,     agcp,    agcm,     agcz
     &                 ,agjp1s
      integer          jproc
       common /amjprc/ jproc
      integer     mxextn
       parameter (mxextn = 10)
      double precision  amass1(mxextn), amass2(mxextn)
       common /kmmass/  amass1,         amass2
      integer           kcharg(mxextn), kfcode(mxextn)
       common /kminfo/  kcharg,         kfcode
* --------------------------------------------------
* -------------------- Functions --------------------
      integer    PYCOMP, KFencode
       external  PYCOMP, KFencode
* ---------------------------------------------------
* ----------------- Local variables -----------------
      character*2        BEAM, TARGET
      double precision   dummy
      character         title*28
      character*5       e_p, e_X, e_q(12), q_q(12), q_g(12)
      character*9       ll(3)
      character*3       e
      character*2       q(12)
      integer           ISUB
      integer           KF, KC, N_merged, qflv_merged(12), iq
     &                 ,KF_base,KC_base
* ---------------------------------------------------
* ---------------- Cross-sections ----------------
C      double precision  x_sec(500), x_sec_err(500)
C       common /X_SEC/   x_sec,      x_sec_err
* ------------------------------------------------
* ------------ BASES common on its result ------------
      integer           ITG,ITF
      real*4            STIME
      double precision  AVGI,SD,CHI2A
      COMMON /BSRSLT/AVGI,SD,CHI2A,STIME,ITG,ITF
* ----------------------------------------------------
* ---------------- Defining user-defined process ----------------
      include './inc/proc_title.h'
        ISUB = ISUB_start + jproc
        if      ( jproc .LE. num_prc_ela)      then      
           title = e_p//' -> '//e_p//'  '//ll(lpair)
         elseif ( jproc .LE. num_prc_ela_qela) then      
           title = e_p//' -> '//e_X//'  '//ll(lpair)
         elseif ( jproc .LE. num_prc_4f)       then      
           if (merge .EQ. 0) then
             title = e_q(qflv)//' -> '//e_q(qflv)//'  '//ll(lpair)
           else
             call Make_Title_merge(3, lpair, merge, e, title)
           endif   
         elseif ( jproc .LE. num_prc_4f_qq)    then      
           title = q_q(qflv)//' -> '//ll(lpair)
         elseif ( jproc .LE. num_prc)          then      
           title = q_g(qflv)//' -> '//ll(lpair)//'  g'
         else
           write(6,*) '!!!Error in PY_Proc_Def!!!'
           write(6,*) ' ---> unknown process(jproc=',jproc,')'
           write(6,*) ' ---> Good-bye!'
           STOP
        endif
C        x_sec(ISUB)     = AVGI
C        x_sec_err(ISUB) = SD
C        call PYUPin(ISUB,             title,             x_sec(ISUB))
        call PYUPin(ISUB,             title,             AVGI)
        MSUB(ISUB) = 1
      MSEL = 0      
* ------------------ Specification of initial state ------------------
      BEAM   = 'p+'
      P(1,1) = 0.D0                         
      P(1,2) = 0.D0                         
      P(1,3) =  dble(P_p_beam)*1.D-3        
      P(1,5) = amp
      P(1,4) = sqrt(  P(1,1)**2 + P(1,2)**2 + P(1,3)**2 + P(1,5)**2  )
      if      (kfcode(2) .EQ.  11) then
         TARGET = 'e-'
       elseif (kfcode(2) .EQ. -11) then
         TARGET = 'e+'
       else
         write(6,*) '!!!Error in PY_Proc_Def!!!'
         write(6,*) '  ---> KF-code for Particle2 =', kfcode(2)
         write(6,*) '       is not supported.'
         write(6,*) '  ---> Good-bye!'
         STOP
      endif
      P(2,1) = 0.D0                         
      P(2,2) = 0.D0                         
      P(2,3) = -dble(P_e_beam)*1.D-3        
      P(2,5) = amlp(1)
      P(2,4) = sqrt(  P(2,1)**2 + P(2,2)**2 + P(2,3)**2 + P(2,5)**2  )
      PARP(2) = amass1(3)+amass1(4)+amass1(5)+amass1(6)
      call PYinit('FIVE', BEAM, TARGET, dummy)
      call PYfram(1)   
                       
* ------------ Parton Shower Parameters ------------
      MSTP(61) = PS_isr   
      MSTP(71) = PS_fsr   
      MSTJ(41) = PS_bra   
      MSTJ(40) = PS_sup   
      MSTJ(42) = 2   
      MSTJ(43) = 4   
      MSTJ(46) = 3   
                     
      MSTJ(50) = 3   
                     
                     
                     
* ----------------------------------------------------------------
      MSTP(111) = PY_decay   
                             
      MSTP(91) =  Ipri_Pt    
      MSTP(11)  = 0   
      MSTP(131) = 0   
      MSTP(52) = 2   
      MSTP(51) = Ngroup*1000 + Nset   
      MINT(93) = 1000000 + MSTP(51)   
* --------------------------------------------------------------------
* ------------ Re-setting Physical Parameters ------------
C>>> Mass
      KF = 23   ! Z0
      KC = PYCOMP(KF)
       PMAS(KC,1) = amz   ! mass        (in GeV)
       PMAS(KC,2) = agz   ! total width (in GeV)
      KF = 24   ! W+
      KC = PYCOMP(KF)
       PMAS(KC,1) = amw   ! mass        (in GeV)
       PMAS(KC,2) = agw   ! total width (in GeV)
      KF = 25   ! SM Higgs
      KC = PYCOMP(KF)
       PMAS(KC,1) = amh   ! mass        (in GeV)
       PMAS(KC,2) = agh   ! total width (in GeV)
      KF = 11   ! e-
      KC = PYCOMP(KF)
       PMAS(KC,1) = amlp(1)   ! mass        (in GeV)
C       PMAS(KC,2) =    ! total width (in GeV)
      KF = 12   ! nu_e
      KC = PYCOMP(KF)
       PMAS(KC,1) = amnu(1)   ! mass        (in GeV)
C       PMAS(KC,2) =    ! total width (in GeV)
      KF = 13   ! mu-
      KC = PYCOMP(KF)
       PMAS(KC,1) = amlp(2)   ! mass        (in GeV)
C       PMAS(KC,2) =    ! total width (in GeV)
      KF = 14   ! nu_mu
      KC = PYCOMP(KF)
       PMAS(KC,1) = amnu(2)   ! mass        (in GeV)
C       PMAS(KC,2) =    ! total width (in GeV)
      KF = 15   ! tau-
      KC = PYCOMP(KF)
       PMAS(KC,1) = amlp(3)   ! mass        (in GeV)
C       PMAS(KC,2) =    ! total width (in GeV)
      KF = 16   ! nu_tau
      KC = PYCOMP(KF)
       PMAS(KC,1) = amnu(3)   ! mass        (in GeV)
C       PMAS(KC,2) =    ! total width (in GeV)
      KF =  1   ! d
      KC = PYCOMP(KF)
       PMAS(KC,1) = amdq(1)   ! mass        (in GeV)
       PMAS(KC,2) = agdq(1)   ! total width (in GeV)
      KF =  2   ! u
      KC = PYCOMP(KF)
       PMAS(KC,1) = amuq(1)   ! mass        (in GeV)
       PMAS(KC,2) = aguq(1)   ! total width (in GeV)
      KF =  3   ! s
      KC = PYCOMP(KF)
       PMAS(KC,1) = amdq(2)   ! mass        (in GeV)
       PMAS(KC,2) = agdq(2)   ! total width (in GeV)
      KF =  4   ! c
      KC = PYCOMP(KF)
       PMAS(KC,1) = amuq(2)   ! mass        (in GeV)
       PMAS(KC,2) = aguq(2)   ! total width (in GeV)
      KF =  5   ! b
      KC = PYCOMP(KF)
       PMAS(KC,1) = amdq(3)   ! mass        (in GeV)
       PMAS(KC,2) = agdq(3)   ! total width (in GeV)
      KF =  6   ! t
      KC = PYCOMP(KF)
       PMAS(KC,1) = amuq(3)   ! mass        (in GeV)
       PMAS(KC,2) = aguq(3)   ! total width (in GeV)
      KF = 2212 ! p
      KC = PYCOMP(KF)
       PMAS(KC,1) = amp   ! mass        (in GeV)
C       PMAS(KC,2) =    ! total width (in GeV)
      KF = 2112 ! n
      KC = PYCOMP(KF)
       PMAS(KC,1) = amn   ! mass        (in GeV)
C       PMAS(KC,2) =    ! total width (in GeV)
C>>>> For merge-mode, make quark masses to be the same.
      if ((merge .NE. 0).and.(process .GE. 3)) then
        call Get_q_merged(merge, N_merged, qflv_merged)
        if (N_merged .GE. 2) then
          KF_base = KFencode(qflv_merged(1))
          KC_base = PYCOMP(KF_base)
          do iq = 2, N_merged
            KF = KFencode(qflv_merged(iq))
            KC = PYCOMP(KF)
            PMAS(KC,1) =  PMAS(KC_base,1)      ! mass        (in GeV)
          enddo
        endif
      endif
* --------------------------------------------------------
      return
      end
