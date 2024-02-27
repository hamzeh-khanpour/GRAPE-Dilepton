      subroutine HB_init(LUN, NTPYT_flag, NTVEC_flag)
      implicit NONE
* ---------- Argument ----------
      integer  LUN      
      logical  NTPYT_flag, NTVEC_flag
* ------------------------------
* ---------- HBOOK stuff ----------
      integer    NWRD_HLIMIT
       parameter(NWRD_HLIMIT = 500000)
      real*4        P
       common /PAWC/P(NWRD_HLIMIT)
      integer*4      IQUEST(100)
       common /QUEST/IQUEST
* ---------------------------------
* ---------- Ntuple varialbes ----------
      include 'nt_inc.h'
* --------------------------------------
* ------ Local variables ------
      integer  ierr
* -----------------------------
C      include './inc/graepia.h'
* --------- Initialization of HBOOK ----------
      call hlimit(NWRD_HLIMIT)
      IQUEST(10) = IQ10           
      call hbset('BSIZE', LREC, ierr)   
      if (ierr .NE. 0) then
         write(6,*) '!!!ERROR from HBSET in HB_init!!!'
         write(6,*) '  ---> Ierr =', ierr
         write(6,*) '  ---> Good-bye!'
         STOP
      endif
* ------------ Opening a Ntuple_file ------------
      call hropen(LUN , 'grp' , NT_NAME , 'NQ' , LREC , ierr)
      if (ierr .NE. 0) then
         write(6,*) '!!!ERROR from HROPEN in HB_init!!!'
         write(6,*) '  ---> Ierr =', ierr
         write(6,*) '  ---> Good-bye!'
         STOP
      endif
      if (NTPYT_flag) then
        call hbnt(NTID, NT_TITLE, ' ')         
        call hbname(NTID, 'PYTHIA',  Npy,  PYTHIA)
      endif
      call hbnt(NTID+10, 'PROCESS parameters', ' ')
      call hbname(NTID+10, 'INTEGER4', nt_jproc,  NT_PRC_i4)
      call hbname(NTID+10, 'REAL8',    nt_P1_lab, NT_PRC_r8)
      if (NTVEC_flag) then
        call hbnt(NTID+20, 'EVENT variables', ' ')
        call hbname(NTID+20, 'INTEGER4',  nt_Nisr, NT_EVT_i4)
        call hbnt(NTID+30, '4-vectors', ' ')
        call hbname(NTID+30, 'VECTOR',  vec_px, NT_VEC)
        call hbnt(NTID+40, 'KF code', ' ')
        call hbname(NTID+40, 'KF code',  vec_kf, NT_KF)
        call hbnt(NTID+50, 'ISR photons', ' ')
        call hbname(NTID+50, 'ISR', nt_vec_isr(1), NT_ISR)
      endif
      return
      end
* ==========================================================================
      subroutine  FILL_nt(LUN, Ngen, Ievt, merge, Nextn, Nisr
     &                                       ,NTPYT_flag,NTVEC_flag)
      implicit NONE
* ---------- Argument ----------
      integer  LUN, Ngen, Ievt, merge, Nextn, Nisr
      logical  NTPYT_flag, NTVEC_flag
* ------------------------------
* -------- Ntuple variables --------
      include 'nt_inc.h'
* ----------------------------------
* ---------- GRACE stuff ----------
      integer  mxextn
       parameter(mxextn=10)
      double precision  vec(4,mxextn)
       common /sp4vec/  vec
      double precision  amass1(mxextn), amass2(mxextn)
       common /kmmass/  amass1,         amass2
      integer          kcharg(mxextn), kfcode(mxextn)
       common /kminfo/ kcharg,         kfcode
      integer          jproc
       common /amjprc/ jproc
* ---------------------------------
* ------------ BASES common on its result ------------
      integer           ITG,ITF
      real*4            STIME
      double precision  AVGI,SD,CHI2A
      COMMON /BSRSLT/AVGI,SD,CHI2A,STIME,ITG,ITF
* ----------------------------------------------------
* ---------- PYTHIA stuff ----------
      include './inc/py_common.h'
* ----------------------------------
* ---------- Kinematical variables ----------
      double precision  P1_lab,P2_lab, E1_lab,E2_lab  
     &                 ,Pcms_lab(3),Ecms_lab(3)
     &                 ,GAMMAcms_lab(3),BETGAMcms_lab(3)
     &                 ,vec_isr(4)
       common /GEP_LAB/ P1_lab,P2_lab, E1_lab,E2_lab
     &                 ,Pcms_lab,Ecms_lab
     &                 ,GAMMAcms_lab,   BETGAMcms_lab
     &                 ,vec_isr
* -------------------------------------------
* -------- Local variables --------
      integer  i
* ---------------------------------
C      include './inc/graepia.h'
      call hcdir('//grp', ' ')   
******** Initialization of HBOOK *******
      if (Ievt .LE. 1) then
C         call HB_init(LUN, NTPYT_flag, NTVEC_flag)
         nt_jproc = jproc
         nt_merge = merge
         nt_Nextn = Nextn
         nt_Ngen  = Ngen
         nt_P1_lab = P1_lab
         nt_P2_lab = P2_lab
         nt_xsec(1) = AVGI   !!! x-sec
         nt_xsec(2) = SD     !!! error
         call hfnt(NTID+10)
      endif
****************************************
      if ((Ievt .GE. 1).and.(Ievt .LE. Ngen)) then !!!!!!!!!!!!!!
         if (NTPYT_flag) then
            Npy = min(N,N_max)
            do 1000 i=1, Npy          
              px(i) = p(i,1)
              py(i) = p(i,2)
              pz(i) = p(i,3)
              pe(i) = p(i,4)
              pm(i) = p(i,5)
              kf(i)  = K(i,2)
              sta(i) = K(i,1)
              mot(i) = K(i,3)
 1000       continue
            call hfnt(NTID)
         endif 
         if (NTVEC_flag) then
            nt_Nisr = Nisr
            call hfnt(NTID+20)
            do 1100 i=1,Nextn      
              vec_px = vec(1,i)
              vec_py = vec(2,i)
              vec_pz = vec(3,i)
              vec_e  = vec(4,i)
              vec_m  = amass1(i)
              vec_kf = kfcode(i)
              call hfnt(NTID+30)
              call hfnt(NTID+40)
 1100       continue
            do 1200 i=1,Nisr
              nt_vec_isr(1) = vec_isr(1)
              nt_vec_isr(2) = vec_isr(2)
              nt_vec_isr(3) = vec_isr(3)
              nt_vec_isr(4) = vec_isr(4)
              call hfnt(NTID+50)
 1200       continue
         endif 
      endif   
********* Termination of HBOOK *********
      if (Ievt .GE. Ngen) then
        call HB_term(NTPYT_flag, NTVEC_flag)
      endif
****************************************
      return
      end
* ==========================================================================
      subroutine  HB_term(NTPYT_flag, NTVEC_flag)
      implicit NONE
* ---------- Argument ----------
      logical  NTPYT_flag, NTVEC_flag
* ------------------------------
* ---------- Ntuple varialbes ----------
      include 'nt_inc.h'
* --------------------------------------
* ------ Local variables ------
      integer  icycle
* -----------------------------
C      include './inc/graepia.h'
      call hldir('//grp', 'T')   
                                 
      call hcdir('//grp', ' ')   
      if (NTPYT_flag) then
        call hrout(NTID, icycle, ' ')
      endif
      call hrout(NTID+10, icycle, ' ')
      if (NTVEC_flag) then
        call hrout(NTID+20, icycle, ' ')
        call hrout(NTID+30, icycle, ' ')
        call hrout(NTID+40, icycle, ' ')
        call hrout(NTID+50, icycle, ' ')
      endif
      call hrend('grp')               
      return
      end
* ==========================================================================
