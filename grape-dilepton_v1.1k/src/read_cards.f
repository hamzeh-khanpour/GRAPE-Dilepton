***********************************************************
*     Getting Input_parameters from Control_card File     *
*                    written by T.Abe                     *
*                   on Aug. 20 in 1998                    *
***********************************************************
* Modified by T.Abe on May 23, 2002
***********************************************
      subroutine Read_cards(LUN, filename)
      implicit NONE
* -------- Argument --------
      integer        LUN
      character*(*)  filename
* --------------------------
* -------------------- FFREAD stuff --------------------
      integer    NW
       parameter(NW=3000)
      real*4           space(NW)
       common /CFREAD/ space
* ------------------------------------------------------
      include './inc/graepia.h'
      integer         jproc      ! For multi-process
       COMMON /amjprc/jproc
* ------ Local variables ------
      integer   i
* -----------------------------
      write(6,*) ' '
      write(6,*) '======> Start of reading control_cards(Read_Cards)'
      write(6,*) ' '
* =============== Initialization of FFREAD ===============
      call FFinit(NW)   !!! Check that NW is large enough or not !!!
      call FFset('LINP', LUN)    ! LUN is used for input file.
       open(unit=LUN, file=filename, status='OLD', err=80)
      GOTO 90
 80   write(6,*) '!!!Error in Read_Cards!!!'
      write(6,*) '  ---> File with input data cards does not exis!'
      write(6,*) '  ---> Good-bye!'
      STOP
 90   call FFset('SIZE', 8)   ! Up to 8 chars are allowed for one key.
* ========================================================
* =============== Setting DEFAULT values ===============
      KF_Lbeam = -11
      Ebeam_pol(1) = 0.   
      Ebeam_pol(2) = 0.   
      Ebeam_pol(3) = 0.   
      P_e_beam =  27520.
      P_p_beam = 820000.
      process   = 1
      lpair     = 2
      Isr_flag  = 1
      NN_ISR    = -1   ! 27
      ISR_scale = 1
      Factor_ISR= 1.
      qflv      = 1
      merge     = 0
      IQCD_scale = 1
      Ngroup =  5   ! GRV94
      Nset   =  5   !  Leading-Order
      Istrf = 1
      do i = 1, num_gra_flg
        jgra_flag(i) = 0
      enddo
*>>> Electroweak Dilepton Production
      jgra_sel  = 3
      weimj_QED = 1.
      weimj_Z0  = 0.
      weimj_H   = 0.
      am_higgs  = 120.
      ag_higgs  =   0.1
      lee_int_VM = .true.
      helicity_VM(-1) = .true.
      helicity_VM( 0) = .true.
      helicity_VM(+1) = .true.
      lJPprod = .false.
      do i = 1, num_jpgra
        lJPgra(i) = .false.
      enddo
      lYYprod = .false.
      do i = 1, num_yygra
        lYYgra(i) = .false.
      enddo
      mxtry     = 30000
      Ngen      = 100
      PS_isr = 1
      PS_fsr = 1
      PS_bra = 2
      PS_sup = 0
      PY_decay = 1
      Ipri_Pt = 1
      E_min_isr = 0.0001
      x_Range_ME(1)  = -1.
      x_Range_ME(2)  = 1.
      y_Range_ME(1)  = -1.
      y_Range_ME(2)  = 1.
      Q2_Range_ME(1) = -1.
      Q2_Range_ME(2) = 1E20
      W_Range_ME(1)  = -1.
      W_Range_ME(2)  = 1E20
      x_Range_OB(1)  = -1.
      x_Range_OB(2)  = 1.
      y_Range_OB(1)  = -1.
      y_Range_OB(2)  = 1.
      Q2_Range_OB(1) = -1.
      Q2_Range_OB(2) = 1E20
      W_Range_OB(1)  = -1.
      W_Range_OB(2)  = 1E20
      Q2p_cut(1) =  0.
      Q2p_cut(2) =  1.E20
      do i=1,4
        theta_min(i) = 0.
        theta_max(i) = 180.
        E_min(i)     = 0.
        E_max(i)     = 1.E20
        P_min(i)     = 0.
        P_max(i)     = 1.E20
        Pt_min(i)    = 0.
        Pt_max(i)    = 1.E20
      enddo
      lPtMAX_cut = .true.
      PtMAX_cut(1)  =  0.
       PtMAX_cut(2) =  1.E20
      the_PtMAX_cut(1)  =   0.
       the_PtMAX_cut(2) = 180.
      Mass56_cut(1)   = 0.
       Mass56_cut(2)  = 1.E20
      Mass56_cut(3)   = 0.
       Mass56_cut(4)  = 1.E20
      MassELL_cut(1)  = 1.
       MassELL_cut(2) = 1.E20
      MassQLL_cut(1)  = 5.
       MassQLL_cut(2) = 1.E20
      MassQLL_cut(3)  = 5.
       MassQLL_cut(4) = 1.E20
      Ivisi = -1   ! 0
      do i=1,4
        the_visi_min(i) =   0.
        the_visi_max(i) = 180.
        ene_visi_min(i) =   0.
        ene_visi_max(i) = 1.E20
        Pt_visi_min(i)  =   0.
        Pt_visi_max(i)  = 1.E20
      enddo
      lveto = .false.
      Iveto = 2
      do i=1,4
        the_veto_min(i) =   0.
        the_veto_max(i) = 180.
        ene_veto_min(i) =   0.
        ene_veto_max(i) = 1.E20
        Pt_veto_min(i)  =   0.
        Pt_veto_max(i)  = 1.E20
      enddo

      l3D_cut = .false.
      do i = 1, 3
         the_3D_min(i) =   0.
         the_3D_max(i) = 180.
         E_3D_min(i) = 0.
         E_3D_max(i) = 1.E20
         P_3D_min(i) = 0.
         P_3D_max(i) = 1.E20
         Pt_3D_min(i) = 0.
         Pt_3D_max(i) = 1.E20
      enddo
      Mass12_3D(1) = 0.
      Mass12_3D(2) = 1.E20
      Q2_3D_ME(1) = 0.
      Q2_3D_ME(2) = 1.E20
      Q2_3D_OB(1) = 0.
      Q2_3D_OB(2) = 1.E20
      W_3D_ME(1) = 0.
      W_3D_ME(2) = 1.E20
      W_3D_OB(1) = 0.
      W_3D_OB(2) = 1.E20
      X_3D_ME(1) = 0.
      X_3D_ME(2) = 1.
      X_3D_OB(1) = 0.
      X_3D_OB(2) = 1.
      Y_3D_ME(1) = 0.
      Y_3D_ME(2) = 1.
      Y_3D_OB(1) = 0.
      Y_3D_OB(2) = 1.

      l3E_cut = .false.
      do i = 1, 3
         the_3E_min(i) =   0.
         the_3E_max(i) = 180.
         E_3E_min(i) = 0.
         E_3E_max(i) = 1.E20
         P_3E_min(i) = 0.
         P_3E_max(i) = 1.E20
         Pt_3E_min(i) = 0.
         Pt_3E_max(i) = 1.E20
      enddo
      Mass12_3E(1) = 0.
      Mass12_3E(2) = 1.E20
      Q2_3E_ME(1) = 0.
      Q2_3E_ME(2) = 1.E20
      Q2_3E_OB(1) = 0.
      Q2_3E_OB(2) = 1.E20
      W_3E_ME(1) = 0.
      W_3E_ME(2) = 1.E20
      W_3E_OB(1) = 0.
      W_3E_OB(2) = 1.E20
      X_3E_ME(1) = 0.
      X_3E_ME(2) = 1.
      X_3E_OB(1) = 0.
      X_3E_OB(2) = 1.
      Y_3E_ME(1) = 0.
      Y_3E_ME(2) = 1.
      Y_3E_OB(1) = 0.
      Y_3E_OB(2) = 1.

      l3F_cut = .false.
      do i = 1, 3
         the_3F_min(i) =   0.
         the_3F_max(i) = 180.
         E_3F_min(i) = 0.
         E_3F_max(i) = 1.E20
         P_3F_min(i) = 0.
         P_3F_max(i) = 1.E20
         Pt_3F_min(i) = 0.
         Pt_3F_max(i) = 1.E20
      enddo
      Mass12_3F(1) = 0.
      Mass12_3F(2) = 1.E20
      Q2_3F_ME(1) = 0.
      Q2_3F_ME(2) = 1.E20
      Q2_3F_OB(1) = 0.
      Q2_3F_OB(2) = 1.E20
      W_3F_ME(1) = 0.
      W_3F_ME(2) = 1.E20
      W_3F_OB(1) = 0.
      W_3F_OB(2) = 1.E20
      X_3F_ME(1) = 0.
      X_3F_ME(2) = 1.
      X_3F_OB(1) = 0.
      X_3F_OB(2) = 1.
      Y_3F_ME(1) = 0.
      Y_3F_ME(2) = 1.
      Y_3F_OB(1) = 0.
      Y_3F_OB(2) = 1.


      lHOT_flag = .true.
      l2e_visiA_flag = .false.    
      the_2e(1) =   0.
      the_2e(2) = 120.
      the_2e(3) = 180.
      E_min_2e(1) = 0.
      E_min_2e(2) = 0.
      l2e_visiB_flag = .false.    
      the_2eB(1)     =   0.
       the_2eB(2)    =   0.
        the_2eB(3)   = 117.
         the_2eB(4)  = 180.
          the_2eB(5) = 180.
      E_min_2eB(1)    = 0.
       E_min_2eB(2)   = 0.
        E_min_2eB(3)  = 0.
         E_min_2eB(4) = 0.
      l3e_visi_flag = .false.    
      the_3e(1) = 0.
      the_3e(2) = 0.
      the_3e(3) = 180.
      the_3e(4) = 180.
      E_min_3e  = 0.
      lPt_visi_flag = .false.    
      the_min_pt(1) =   0.
      the_min_pt(2) =   0.
      the_min_pt(3) =   0.
      the_max_pt(1) = 180.
      the_max_pt(2) = 180.
      the_max_pt(3) = 180.
      Pt_min_pt(1)  =   0.
      Pt_min_pt(2)  =   0.
      Pt_min_pt(3)  =   0.
      lscattL_flag  = .false.   
      E_scattL(1) = 0.
      E_scattL(2) = 1.E20
        P_prodL(1) = 0.
        P_prodL(2) = 1.E20
      theta_scattL(1) = 0.
      theta_scattL(2) = 180.
        theta_prodL(1) = 0.
        theta_prodL(2) = 180.
      u_cut(1) = -999.
      u_cut(2) =  1.E20
      W_cut(1) =  1.08
      W_cut(2) =  1.E20
* ------ Quasi-elastic ------
      Wmin     =  1.08
      Wmax     =  1.E20
      rk    = 0.05         
      A_NBD = 0.176
      B_NBD = 0.43132
      C_NBD = 0.86224
      P_u      = 0.5       
      P_d      = 0.5
      P_s      = 0.07
      P_ud_bar = 0.02
      Supp_pi0 = 0.4
      A_slope = 0.00175    
      B_slope = 0.00167
      C_slope = 0.353
      red_slope_s = 0.6
* ---------------------------
      isym_34 = 1122   ! 1
      ii34    = 3
      Iscale_xq = 2
      Ireso56  = 1122   ! -1
      thresh56 = 0.9
      Rnn_MJ56 = 2.
      Ireso456    = -1
      Rnn_456     =  2.
      Mas_reso456 = 0.
      Wid_reso456 = 0.
      IcosP3 = 2
      IcosP5 = 2
      Nregion  = 2
      Neps_p = -1
      Rscale_Mn = 1.
      num_it_grid  =  4
      num_it_integ = 10
      acc1_card = 0.2
      acc2_card = 0.01
      num_call = 1000000
      LIST_flag = .true.
      NTPYT_flag = .false.
      NTVEC_flag = .false.
      ASC_flag  = .false.
      Nmod      = 1000
      Nlist     = 10
      ROT_kinem_flag  = .false.
      ROT_pyrand_flag = .true.
      Qela_decay = 1
      Wrt_cards = .false.
      LRND_flag = .false.
      Frame_amp = 1   ! 1:ep, 2:gp
* ======================================================
* ============ Reading user-defined control_cards ============
      call FFkey('KFLBEAM',  KF_Lbeam,  1,  'INTE')
      call FFkey('EBEAM',   P_e_beam, 1, 'REAL')
      call FFkey('PBEAM',   P_p_beam, 1, 'REAL')
      call FFkey('EPOL',  Ebeam_pol,  3, 'REAL')
      call FFkey('PROCESS', process,   1, 'INTE')
      call FFkey('LPAIR',   lpair,     1, 'INTE')
      call FFkey('ISR',     Isr_flag,  1, 'INTE')
      call FFkey('NNISR',     NN_ISR,  1, 'INTE')
      call FFkey('ISRSCALE', ISR_scale,  1, 'INTE')
      call FFkey('FACTISR',  Factor_ISR,  1, 'REAL')
      call FFkey('QFLV',    qflv,      1, 'INTE')
      call FFkey('MERGE',       merge, 1, 'INTE')
      call FFkey('QCDSCALE',  IQCD_scale, 1, 'INTE')
      call FFkey('NGROUP', Ngroup,  1, 'INTE')
      call FFkey('NSET',   Nset,    1, 'INTE')
      call FFkey('STRF',   Istrf,    1, 'INTE')
      call FFkey('GRAFLG',  jgra_flag, num_gra_flg, 'INTE')
*>>> Electroweak Dilepton Production
      call FFkey('GRASEL',   jgra_sel,    1, 'INTE')
      call FFkey('CNTMJWEI', weimj_QED,   1, 'REAL')
      call FFkey('Z0MJWEI',  weimj_Z0,    1, 'REAL')
      call FFkey('H0MJWEI',   weimj_H,     1, 'REAL')
      call FFkey('H0MASS',    am_higgs,    1, 'REAL')
      call FFkey('H0WIDTH',   ag_higgs,    1, 'REAL')
      call FFkey('VMMODEL',   Model_VM,     1,  'INTE')
      call FFkey('VMTSLOPE',  t_slope_VM,   1,  'REAL')
      call FFkey('VMEEINT',   lee_int_VM,   1,  'LOGICAL')
      call FFkey('VMHELI',    helicity_VM(-1),   3,  'LOGICAL')
      call FFkey('JPPROD',   lJPprod,          1,  'LOGICAL')
      call FFkey('JPTYPE',   lJPtype,    num_jptype, 'LOGICAL')
      call FFkey('JPBCORR',    bcorr_jp, num_jptype, 'REAL')
      call FFkey('JPMJWEI',    weimj_jp, num_jptype, 'REAL')
      call FFkey('JPPHASEP',  phaseP_jp, num_jptype, 'REAL')
      call FFkey('JPPHASEG',  phaseG_jp, num_jptype, 'REAL')
      call FFkey('JPGRAPH',   lJPgra,  num_jpgra, 'LOGICAL')
      call FFkey('YYPROD',   lYYprod,          1,  'LOGICAL')
      call FFkey('YYTYPE',   lYYtype,    num_yytype, 'LOGICAL')
      call FFkey('YYBCORR',    bcorr_yy, num_yytype, 'REAL')
      call FFkey('YYMJWEI',    weimj_yy, num_yytype, 'REAL')
      call FFkey('YYPHASEP',  phaseP_yy, num_yytype, 'REAL')
      call FFkey('YYPHASEG',  phaseG_yy, num_yytype, 'REAL')
      call FFkey('YYGRAPH',   lYYgra,  num_yygra, 'LOGICAL')
      call FFkey('MXTRY',       mxtry, 1, 'INTE')
      call FFkey('NGEN',         Ngen, 1, 'INTE')
      call FFkey('PSISR',      PS_isr, 1, 'INTE')
      call FFkey('PSFSR',      PS_fsr, 1, 'INTE')
      call FFkey('PSBRA',      PS_bra, 1, 'INTE')
      call FFkey('PSSUP',      PS_sup, 1, 'INTE')
      call FFkey('PYDECAY',  PY_decay, 1, 'INTE')
      call FFkey('PRIPT',     Ipri_Pt, 1, 'INTE')
      call FFkey('XXRNGME',    x_Range_ME,   2, 'REAL')
      call FFkey('YYRNGME',    y_Range_ME,   2, 'REAL')
      call FFkey('Q2RNGME',   Q2_Range_ME,   2, 'REAL')
      call FFkey('WWRNGME',    W_Range_ME,   2, 'REAL')
      call FFkey('XXRNGOB',    x_Range_OB,   2, 'REAL')
      call FFkey('YYRNGOB',    y_Range_OB,   2, 'REAL')
      call FFkey('Q2RNGOB',   Q2_Range_OB,   2, 'REAL')
      call FFkey('WWRNGOB',    W_Range_OB,   2, 'REAL')
      call FFkey('Q2P',    Q2p_cut,    2, 'REAL')
      call FFkey('THMIN',  theta_min,  4, 'REAL')
      call FFkey('THMAX',  theta_max,  4, 'REAL')
      call FFkey('EMIN',   E_min,      4, 'REAL')
      call FFkey('EMAX',   E_max,      4, 'REAL')
      call FFkey('PMIN',   P_min,      4, 'REAL')
      call FFkey('PMAX',   P_max,      4, 'REAL')
      call FFkey('PTMIN',  Pt_min,     4, 'REAL')
      call FFkey('PTMAX',  Pt_max,     4, 'REAL')
      call FFKey('PTMCTFLG',   lPtMAX_cut,  1, 'LOGICAL')
      call FFkey('PTMXCT',      PtMAX_cut,  2, 'REAL')
      call FFkey('THPTMCT', the_PtMAX_cut,  2, 'REAL')
      call FFkey('MASSLL',  Mass56_cut,  4, 'REAL')
      call FFkey('MASSELL', MassELL_cut, 2, 'REAL')
      call FFkey('MASSQLL', MassQLL_cut, 4, 'REAL')
      call FFkey('IVISI',   Ivisi,  1,   'INTE')
      call FFkey('THEVMIN', the_visi_min, 4, 'REAL')
      call FFkey('THEVMAX', the_visi_max, 4, 'REAL')
      call FFkey('EVMIN',   ene_visi_min, 4, 'REAL')
      call FFkey('EVMAX',   ene_visi_max, 4, 'REAL')
      call FFkey('PTVMIN',   Pt_visi_min, 4, 'REAL')
      call FFkey('PTVMAX',   Pt_visi_max, 4, 'REAL')
      call FFkey('LVETO',   lveto,  1,   'LOGICAL')
      call FFkey('IVETO',   Iveto,  1,   'INTE')
      call FFkey('THEVTMIN', the_veto_min, 4, 'REAL')
      call FFkey('THEVTMAX', the_veto_max, 4, 'REAL')
      call FFkey('EVTMIN',   ene_veto_min, 4, 'REAL')
      call FFkey('EVTMAX',   ene_veto_max, 4, 'REAL')
      call FFkey('PTVTMIN',   Pt_veto_min, 4, 'REAL')
      call FFkey('PTVTMAX',   Pt_veto_max, 4, 'REAL')

      call FFkey('L3DCUT',   l3D_cut,       1, 'LOGICAL')
      call FFkey('TH3DMIN',  the_3D_min,    3, 'REAL')
      call FFkey('TH3DMAX',  the_3D_max,    3, 'REAL')
      call FFkey('E3DMIN',   E_3D_min,      3, 'REAL')
      call FFkey('E3DMAX',   E_3D_max,      3, 'REAL')
      call FFkey('P3DMIN',   P_3D_min,      3, 'REAL')
      call FFkey('P3DMAX',   P_3D_max,      3, 'REAL')
      call FFkey('PT3DMIN',  Pt_3D_min,     3, 'REAL')
      call FFkey('PT3DMAX',  Pt_3D_max,     3, 'REAL')
      call FFkey('M12CUT3D',  Mass12_3D,    2, 'REAL')
      call FFkey('Q2CT3DME', Q2_3D_ME,     2, 'REAL')
      call FFkey('Q2CT3DOB', Q2_3D_OB,     2, 'REAL')
      call FFkey('WCT3DME', W_3D_ME,     2, 'REAL')
      call FFkey('WCT3DOB', W_3D_OB,     2, 'REAL')
      call FFkey('XCT3DME', X_3D_ME,     2, 'REAL')
      call FFkey('XCT3DOB', X_3D_OB,     2, 'REAL')
      call FFkey('YCT3DME', Y_3D_ME,     2, 'REAL')
      call FFkey('YCT3DOB', Y_3D_OB,     2, 'REAL')

      call FFkey('L3ECUT',   l3E_cut,       1, 'LOGICAL')
      call FFkey('TH3EMIN',  the_3E_min,    3, 'REAL')
      call FFkey('TH3EMAX',  the_3E_max,    3, 'REAL')
      call FFkey('E3EMIN',   E_3E_min,      3, 'REAL')
      call FFkey('E3EMAX',   E_3E_max,      3, 'REAL')
      call FFkey('P3EMIN',   P_3E_min,      3, 'REAL')
      call FFkey('P3EMAX',   P_3E_max,      3, 'REAL')
      call FFkey('PT3EMIN',  Pt_3E_min,     3, 'REAL')
      call FFkey('PT3EMAX',  Pt_3E_max,     3, 'REAL')
      call FFkey('M12CUT3E',  Mass12_3E,    2, 'REAL')
      call FFkey('Q2CT3EME', Q2_3E_ME,     2, 'REAL')
      call FFkey('Q2CT3EOB', Q2_3E_OB,     2, 'REAL')
      call FFkey('WCT3EME', W_3E_ME,     2, 'REAL')
      call FFkey('WCT3EOB', W_3E_OB,     2, 'REAL')
      call FFkey('XCT3EME', X_3E_ME,     2, 'REAL')
      call FFkey('XCT3EOB', X_3E_OB,     2, 'REAL')
      call FFkey('YCT3EME', Y_3E_ME,     2, 'REAL')
      call FFkey('YCT3EOB', Y_3E_OB,     2, 'REAL')

      call FFkey('L3FCUT',   l3F_cut,       1, 'LOGICAL')
      call FFkey('TH3FMIN',  the_3F_min,    3, 'REAL')
      call FFkey('TH3FMAX',  the_3F_max,    3, 'REAL')
      call FFkey('E3FMIN',   E_3F_min,      3, 'REAL')
      call FFkey('E3FMAX',   E_3F_max,      3, 'REAL')
      call FFkey('P3FMIN',   P_3F_min,      3, 'REAL')
      call FFkey('P3FMAX',   P_3F_max,      3, 'REAL')
      call FFkey('PT3FMIN',  Pt_3F_min,     3, 'REAL')
      call FFkey('PT3FMAX',  Pt_3F_max,     3, 'REAL')
      call FFkey('M12CUT3F',  Mass12_3F,    2, 'REAL')
      call FFkey('Q2CT3FME', Q2_3F_ME,     2, 'REAL')
      call FFkey('Q2CT3FOB', Q2_3F_OB,     2, 'REAL')
      call FFkey('WCT3FME', W_3F_ME,     2, 'REAL')
      call FFkey('WCT3FOB', W_3F_OB,     2, 'REAL')
      call FFkey('XCT3FME', X_3F_ME,     2, 'REAL')
      call FFkey('XCT3FOB', X_3F_OB,     2, 'REAL')
      call FFkey('YCT3FME', Y_3F_ME,     2, 'REAL')
      call FFkey('YCT3FOB', Y_3F_OB,     2, 'REAL')

      call FFKey('HOTFLG',    lHOT_flag, 1, 'LOGICAL')
      call FFKey('L2EVISIA', l2e_visiA_flag, 1, 'LOGICAL')    
      call FFkey('THE2E',    the_2e,     3, 'REAL')
      call FFkey('EMIN2E',   E_min_2e,   2, 'REAL')
      call FFKey('L2EVISIB', l2e_visiB_flag, 1, 'LOGICAL')    
      call FFkey('THE2EB',    the_2eB,   5, 'REAL')
      call FFkey('EMIN2EB',   E_min_2eB, 4, 'REAL')
      call FFKey('L3EVISI',  l3e_visi_flag,  1, 'LOGICAL')    
      call FFkey('THE3E',    the_3e,     4, 'REAL')
      call FFkey('EMIN3E',   E_min_3e,      1, 'REAL')
      call FFKey('LPTVISI',  lPt_visi_flag,  1, 'LOGICAL')    
      call FFkey('THEMINPT', the_min_pt, 3, 'REAL')
      call FFkey('THEMAXPT', the_max_pt, 3, 'REAL')
      call FFkey('PTMINPT',  Pt_min_pt,  3, 'REAL')
      call FFKey('LSCATL',  lscattL_flag,   1, 'LOGICAL')    
      call FFkey('ESCATL',   E_scattL,     2, 'REAL')
      call FFkey('THESCATL', theta_scattL, 2, 'REAL')
      call FFkey('PPRODL',   P_prodL,     2, 'REAL')
      call FFkey('THEPRODL', theta_prodL, 2, 'REAL')
      call FFkey('UCUT', u_cut, 2, 'REAL')
      call FFkey('MHAD', W_cut, 2, 'REAL')
      call FFkey('RK',           rk,  1, 'REAL')
      call FFkey('ANBD',      A_NBD,  1, 'REAL')
      call FFkey('BNBD',      B_NBD,  1, 'REAL')
      call FFkey('CNBD',      C_NBD,  1, 'REAL')
      call FFkey('PUQ',         P_u,  1, 'REAL')
      call FFkey('PDQ',         P_d,  1, 'REAL')
      call FFkey('PSQ',         P_s,  1, 'REAL')
      call FFkey('PUD',    P_ud_bar,  1, 'REAL')
      call FFkey('SUPPI0', Supp_pi0,  1, 'REAL')
      call FFkey('ASLOPE',  A_slope,  1, 'REAL')
      call FFkey('BSLOPE',  B_slope,  1, 'REAL')
      call FFkey('CSLOPE',  C_slope,  1, 'REAL')
      call FFkey('REDS', red_slope_s, 1, 'REAL')
      call FFkey('ISYM34', isym_34, 1, 'INTE')
      call FFkey('I34',    ii34,    1, 'INTE')
      call FFkey('SCALEXQ', Iscale_xq, 1, 'INTE')
      call FFkey('RESNS56', Ireso56,  1, 'INTE')
      call FFkey('THRES56',  thresh56, 1, 'REAL')
      call FFkey('NCNT',  Rnn_MJ56, 1, 'REAL')
      call FFkey('RESNS456',Ireso456, 1, 'INTE')
      call FFkey('RNN456', Rnn_456, 1, 'REAL')
      call FFkey('MAS456', Mas_reso456, 1, 'REAL')
      call FFkey('WID456', Wid_reso456, 1, 'REAL')
      call FFkey('ICOS3',    IcosP3,  1, 'INTE')
      call FFkey('ICOS5',    IcosP5,  1, 'INTE')
      call FFkey('NREG',    Nregion,  1, 'INTE')
      call FFkey('NEPSP',    Neps_p,  1, 'INTE')
      call FFkey('RSCALEMN', Rscale_Mn,  1, 'REAL')
      call FFkey('ITMX1',  num_it_grid,  1, 'INTE')
      call FFkey('ITMX2',  num_it_integ, 1, 'INTE')
      call FFkey('ACC1',    acc1_card,   1, 'REAL')
      call FFkey('ACC2',    acc2_card,   1, 'REAL')
      call FFkey('NCALL',  num_call,     1, 'INTE')
      call FFKey('PYLIST',  LIST_flag,  1, 'LOGICAL')
      call FFKey('NTPYT',   NTPYT_flag, 1, 'LOGICAL')
      call FFKey('ASCII',   ASC_flag,   1, 'LOGICAL')
      call FFkey('NMOD',    Nmod,       1, 'INTE')
      call FFkey('NLIST',   Nlist,      1, 'INTE')
      call FFKey('NTVEC',   NTVEC_flag, 1, 'LOGICAL')
      call FFkey('EMINISR', E_min_isr,        1, 'REAL')
      call FFKey('ROTKIN',  ROT_kinem_flag,   1, 'LOGICAL')
      call FFKey('ROTPY',   ROT_pyrand_flag,  1, 'LOGICAL')
      call FFkey('QELAX',   Qela_decay,       1, 'INTE')
      call FFkey('WRTCARDS',   Wrt_cards,      1, 'LOGICAL')
      call FFkey('RNDGEN',     LRND_flag,      1, 'LOGICAL')
      call FFkey('FRAMEAMP',   Frame_amp,      1, 'INTE')
* ============================================================
      call FFgo   
      close(LUN)
      write(6,*) ' '
      write(6,*) '========> End of reading control_cards(Read_Cards)'
      write(6,*) ' '
      write(6,*) ' '
      include 'kinem_auto_tune.h'
      include 'warnings_cards.h'
      return
      end
*###################################################################
      subroutine Write_cards
      implicit NONE
* -------- Argument --------
* --------------------------
* -------------------- FFREAD stuff --------------------
* ------------------------------------------------------
      include './inc/graepia.h'
* ------ Local variables ------
      integer   i, N
* -----------------------------
      if (.not.Wrt_cards)  RETURN
* ============ Printing user-defined control_cards ============
      write(6,*) ' '
      write(6,1100)     '   KFLBEAM     ', KF_Lbeam
      write(6,2100)     '   EBEAM       ', P_e_beam
      write(6,2100)     '   PBEAM       ', P_p_beam
      write(6,2300)     '   EPOL        ', (Ebeam_pol(i),i=1,3)
      write(6,1100)     '   PROCESS     ', process
      write(6,1100)     '   LPAIR       ', lpair
      write(6,1100)     '   ISR         ', Isr_flag
      write(6,1100)     '   NNISR       ', NN_ISR
      write(6,1100)     '   ISRSCALE    ', ISR_scale
      write(6,2100)     '   FACTISR     ', Factor_ISR
      write(6,1100)     '   QFLV        ', qflv
      write(6,1100)     '   MERGE       ', merge
      write(6,1100)     '   QCDSCALE    ', IQCD_scale
      write(6,1100)     '   NGROUP      ', Ngroup
      write(6,1100)     '   NSET        ', Nset
      write(6,1100)     '   STRF        ', Istrf
      N = min(num_gra_flg, 5)
      write(6,1500)     '   GRAFLG      ', (jgra_flag(i),i=1,N)
      if (num_gra_flg .GT. 5) then
         N = min(num_gra_flg, 10)
         write(6,1500)  '               ', (jgra_flag(i),i=6,N)
      endif
      if (num_gra_flg .GT.10) then
         N = min(num_gra_flg, 15)
         write(6,1500)  '               ', (jgra_flag(i),i=11,N)
      endif
      if (num_gra_flg .GT.15) then
         N = min(num_gra_flg, 20)
         write(6,1500)  '               ', (jgra_flag(i),i=16,N)
      endif
      if (num_gra_flg .GT.20) then
         N = num_gra_flg   ! min(num_gra_flg, 25)
         write(6,1500)  '               ', (jgra_flag(i),i=21,N)
      endif
      write(6,*)   '======= Electroweak Dilepton Production ======='
      write(6,1100)     '  GRASEL       ', jgra_sel
      write(6,2100)     '  QEDMJWEI     ', weimj_QED
      write(6,2100)     '  Z0MJWEI      ', weimj_Z0
      write(6,2100)     '  HMJWEI       ', weimj_H
      write(6,2100)     '  HMASS        ', am_higgs
      write(6,2100)     '  HWIDTH       ', ag_higgs
      write(6,*) ' '
      write(6,1100)     '  VMMODEL      ', Model_VM
      write(6,2100)     '  VMTSLOPE     ', t_slope_VM
      write(6,100)      '  VMEEINT      ', lee_int_VM
      write(6,130)      '  VMHELI       ', (helicity_VM(i),i=-1,+1)
      write(6,*) ' '
      write(6,100)      '  JPPROD       ', lJPprod
      write(6,*)        '  JPTYPE       ', (lJPtype(i),i=1,num_jptype)
      write(6,*)        '  JPBCORR      ', (bcorr_jp(i),i=1,num_jptype)
      write(6,*)        '  JPMJWEI      ', (weimj_jp(i),i=1,num_jptype)
      write(6,*)        '  JPPHASEP     ', (phaseP_jp(i),i=1,num_jptype)
      write(6,*)        '  JPPHASEG     ', (phaseG_jp(i),i=1,num_jptype)
      write(6,*)        '  JPGRAPH      ', (lJPgra(i),i=1,num_jpgra)
      write(6,*) ' '
      write(6,100)      '  YYPROD       ', lYYprod
      write(6,*)        '  YYTYPE       ', (lYYtype(i),i=1,num_yytype)
      write(6,*)        '  YYBCORR      ', (bcorr_yy(i),i=1,num_yytype)
      write(6,*)        '  YYMJWEI      ', (weimj_yy(i),i=1,num_yytype)
      write(6,*)        '  YYPHASEP     ', (phaseP_yy(i),i=1,num_yytype)
      write(6,*)        '  YYPHASEG     ', (phaseG_yy(i),i=1,num_yytype)
      write(6,*)        '  YYGRAPH      ', (lYYgra(i),i=1,num_yygra)
      write(6,*)   '========== BASES parameters =========='
      write(6,1100)     '  ITMX1        ', num_it_grid
      write(6,1100)     '  ITMX2        ', num_it_integ
      write(6,2100)     '  ACC1         ', acc1_card
      write(6,2100)     '  ACC2         ', acc2_card
      write(6,1100)     '  NCALL        ', num_call
      write(6,1100)     '  SCALEXQ      ', Iscale_xq
      write(6,1100)     '  ISYM34       ', isym_34
      write(6,1100)     '  I34          ', ii34
      write(6,1100)     '  RESNS56      ', Ireso56
      write(6,2100)     '  NNMJ56       ', Rnn_MJ56
      write(6,2100)     '  THRES56      ', thresh56
      write(6,1100)     '  RESNS456     ', Ireso456
      write(6,2100)     '  RNN456       ', Rnn_456
      write(6,2100)     '  MAS456       ', Mas_reso456
      write(6,2100)     '  WID456       ', Wid_reso456
      write(6,1100)     '  ICOS3        ', IcosP3
      write(6,1100)     '  ICOS5        ', IcosP5
      write(6,1100)     '  NREG         ', Nregion
      write(6,1100)     '  NEPSP        ', Neps_p
      write(6,2100)     '  RSCALEMN     ', Rscale_Mn
      write(6,*)   '======================================'
      write(6,*)   '========== SPRING parameters =========='
      write(6,1100)     '  NGEN         ', Ngen
      write(6,1100)     '  NMOD         ', Nmod
      write(6,1100)     '  MXTRY        ', mxtry
      write(6,1100)     '  PSISR        ', PS_isr
      write(6,1100)     '  PSFSR        ', PS_fsr
      write(6,1100)     '  PSBRA        ', PS_bra
      write(6,1100)     '  PSSUP        ', PS_sup
      write(6,1100)     '  PYDECAY      ', PY_decay
      write(6,1100)     '  PRIPT        ', Ipri_Pt
      write(6,*)   '======================================'
      write(6,*)   '========== Output of Generated Events =========='
      write(6,100)      '  PYLIST       ', LIST_flag
      write(6,1100)     '  NLIST        ', Nlist
      write(6,100)      '  NTPYT        ', NTPYT_flag
      write(6,100)      '  NTVEC        ', NTVEC_flag
      write(6,100)      '  ASCII        ', ASC_flag
      write(6,*)   '================================================'
      write(6,*)   '================== Misc. =================='
      write(6,2100)     '  EMINISR      ', E_min_isr
      write(6,100)      '  ROTKIN       ', ROT_kinem_flag
      write(6,100)      '  ROTPY        ', ROT_pyrand_flag
      write(6,1100)     '  QELAX        ', Qela_decay
      write(6,100)      '  WRTCARDS     ', Wrt_cards
      write(6,100)      '  RNDGEN       ', LRND_flag
      write(6,1100)     '  FRAMEAMP     ', Frame_amp
      write(6,*)   '==========================================='
      write(6,*)   '=============== << Cuts >> ==============='
      write(6,2200)     '  XXRNGME      ',  (x_Range_ME(i),i=1,2)
      write(6,2200)     '  XXRNGOB      ',  (x_Range_OB(i),i=1,2)
      write(6,2200)     '  YYRNGME      ',  (y_Range_ME(i),i=1,2)
      write(6,2200)     '  YYRNGOB      ',  (y_Range_OB(i),i=1,2)
      write(6,2200)     '  Q2RNGME      ', (Q2_Range_ME(i),i=1,2)
      write(6,2200)     '  Q2RNGOB      ', (Q2_Range_OB(i),i=1,2)
      write(6,2200)     '  WWRNGME      ',  (W_Range_ME(i),i=1,2)
      write(6,2200)     '  WWRNGOB      ',  (W_Range_OB(i),i=1,2)
      write(6,*) ' ---------- '
      write(6,2200)     '  MHADCUT      ', (W_cut(i),i=1,2)
      write(6,*) ' ---------- '
      write(6,2200)     '  Q2P          ', (Q2p_cut(i),i=1,2)
      write(6,*) ' ---------- '
      write(6,2200)     '  UCUT         ', (u_cut(i),i=1,2)
      write(6,*) ' ---------- '
      write(6,2400)     '  THMIN        ', (theta_min(i),i=1,4)
      write(6,2400)     '  THMAX        ', (theta_max(i),i=1,4)
      write(6,2400)     '  EMIN         ', (E_min(i),    i=1,4)
      write(6,2400)     '  EMAX         ', (E_max(i),    i=1,4)
      write(6,2400)     '  PMIN         ', (P_min(i),    i=1,4)
      write(6,2400)     '  PMAX         ', (P_max(i),    i=1,4)
      write(6,2400)     '  PTMIN        ', (Pt_min(i),   i=1,4)
      write(6,2400)     '  PTMAX        ', (Pt_max(i),   i=1,4)
      write(6,*) ' ---------- '
      write(6,100)      '  PTMCTFLG     ', lPtMAX_cut
      write(6,2200)     '  THPTMCT      ', (the_PtMAX_cut(i),i=1,2)
      write(6,2200)     '  PTMXCT       ', (PtMAX_cut(i),i=1,2)
      write(6,*) ' ---------- '
      write(6,2200)     '  MASSLL       ', (Mass56_cut(i),i=1,2)
      write(6,2200)     '               ', (Mass56_cut(i),i=3,4)
      write(6,2200)     '  MASSELL      ', (MassELL_cut(i),i=1,2)
      write(6,2200)     '  MASSQLL      ', (MassQLL_cut(i),i=1,2)
      write(6,2200)     '               ', (MassQLL_cut(i),i=3,4)
      write(6,*)   ' ------ Hot selection ------ '
      write(6,100)      '  HOTFLG       ', lHOT_flag
      write(6,100)      '   LVETO       ', lveto
      write(6,1100)     '   IVETO       ', Iveto
      write(6,2400)     '    THEVTMIN   ', (the_veto_min(i), i=1,4)
      write(6,2400)     '    THEVTMAX   ', (the_veto_max(i), i=1,4)
      write(6,2400)     '    EVTMIN     ', (ene_veto_min(i), i=1,4)
      write(6,2400)     '    EVTMAX     ', (ene_veto_max(i), i=1,4)
      write(6,2400)     '    PTVTMIN    ', (Pt_veto_min(i),  i=1,4)
      write(6,2400)     '    PTVTMAX    ', (Pt_veto_max(i),  i=1,4)
      write(6,1100)     '  IVISI        ', Ivisi
      write(6,2400)     '    THEVMIN    ', (the_visi_min(i), i=1,4)
      write(6,2400)     '    THEVMAX    ', (the_visi_max(i), i=1,4)
      write(6,2400)     '    EVMIN      ', (ene_visi_min(i), i=1,4)
      write(6,2400)     '    EVMAX      ', (ene_visi_max(i), i=1,4)
      write(6,2400)     '    PTVMIN     ', (Pt_visi_min(i),  i=1,4)
      write(6,2400)     '    PTVMAX     ', (Pt_visi_max(i),  i=1,4)

      write(6,100)      '   L3DCUT      ', l3D_cut
      write(6,2300)     '    TH3DMIN    ', (the_3D_min(i), i=1,3)
      write(6,2300)     '    TH3DMAX    ', (the_3D_max(i), i=1,3)
      write(6,2300)     '    E3DMIN     ', (E_3D_min(i),   i=1,3)
      write(6,2300)     '    E3DMAX     ', (E_3D_max(i),   i=1,3)
      write(6,2300)     '    P3DMIN     ', (P_3D_min(i),   i=1,3)
      write(6,2300)     '    P3DMAX     ', (P_3D_max(i),   i=1,3)
      write(6,2300)     '    PT3DMIN    ', (Pt_3D_min(i),  i=1,3)
      write(6,2300)     '    PT3DMAX    ', (Pt_3D_max(i),  i=1,3)
      write(6,2200)     '    M12CUT3D   ', (Mass12_3D(i),  i=1,2)
      write(6,2200)     '    Q2CT3DME   ', (Q2_3D_ME(i),   i=1,2)
      write(6,2200)     '    Q2CT3DOB   ', (Q2_3D_OB(i),   i=1,2)
      write(6,2200)     '    WCT3DME    ', (W_3D_ME(i),   i=1,2)
      write(6,2200)     '    WCT3DOB    ', (W_3D_OB(i),   i=1,2)
      write(6,2200)     '    XCT3DME    ', (X_3D_ME(i),   i=1,2)
      write(6,2200)     '    XCT3DOB    ', (X_3D_OB(i),   i=1,2)
      write(6,2200)     '    YCT3DME    ', (Y_3D_ME(i),   i=1,2)
      write(6,2200)     '    YCT3DOB    ', (Y_3D_OB(i),   i=1,2)

      write(6,100)      '   L3ECUT      ', l3E_cut
      write(6,2300)     '    TH3EMIN    ', (the_3E_min(i), i=1,3)
      write(6,2300)     '    TH3EMAX    ', (the_3E_max(i), i=1,3)
      write(6,2300)     '    E3EMIN     ', (E_3E_min(i),   i=1,3)
      write(6,2300)     '    E3EMAX     ', (E_3E_max(i),   i=1,3)
      write(6,2300)     '    P3EMIN     ', (P_3E_min(i),   i=1,3)
      write(6,2300)     '    P3EMAX     ', (P_3E_max(i),   i=1,3)
      write(6,2300)     '    PT3EMIN    ', (Pt_3E_min(i),  i=1,3)
      write(6,2300)     '    PT3EMAX    ', (Pt_3E_max(i),  i=1,3)
      write(6,2200)     '    M12CUT3E   ', (Mass12_3E(i),  i=1,2)
      write(6,2200)     '    Q2CT3EME   ', (Q2_3E_ME(i),   i=1,2)
      write(6,2200)     '    Q2CT3EOB   ', (Q2_3E_OB(i),   i=1,2)
      write(6,2200)     '    WCT3EME    ', (W_3E_ME(i),   i=1,2)
      write(6,2200)     '    WCT3EOB    ', (W_3E_OB(i),   i=1,2)
      write(6,2200)     '    XCT3EME    ', (X_3E_ME(i),   i=1,2)
      write(6,2200)     '    XCT3EOB    ', (X_3E_OB(i),   i=1,2)
      write(6,2200)     '    YCT3EME    ', (Y_3E_ME(i),   i=1,2)
      write(6,2200)     '    YCT3EOB    ', (Y_3E_OB(i),   i=1,2)

      write(6,100)      '   L3FCUT      ', l3F_cut
      write(6,2300)     '    TH3FMIN    ', (the_3F_min(i), i=1,3)
      write(6,2300)     '    TH3FMAX    ', (the_3F_max(i), i=1,3)
      write(6,2300)     '    E3FMIN     ', (E_3F_min(i),   i=1,3)
      write(6,2300)     '    E3FMAX     ', (E_3F_max(i),   i=1,3)
      write(6,2300)     '    P3FMIN     ', (P_3F_min(i),   i=1,3)
      write(6,2300)     '    P3FMAX     ', (P_3F_max(i),   i=1,3)
      write(6,2300)     '    PT3FMIN    ', (Pt_3F_min(i),  i=1,3)
      write(6,2300)     '    PT3FMAX    ', (Pt_3F_max(i),  i=1,3)
      write(6,2200)     '    M12CUT3F   ', (Mass12_3F(i),  i=1,2)
      write(6,2200)     '    Q2CT3FME   ', (Q2_3F_ME(i),   i=1,2)
      write(6,2200)     '    Q2CT3FOB   ', (Q2_3F_OB(i),   i=1,2)
      write(6,2200)     '    WCT3FME    ', (W_3F_ME(i),   i=1,2)
      write(6,2200)     '    WCT3FOB    ', (W_3F_OB(i),   i=1,2)
      write(6,2200)     '    XCT3FME    ', (X_3F_ME(i),   i=1,2)
      write(6,2200)     '    XCT3FOB    ', (X_3F_OB(i),   i=1,2)
      write(6,2200)     '    YCT3FME    ', (Y_3F_ME(i),   i=1,2)
      write(6,2200)     '    YCT3FOB    ', (Y_3F_OB(i),   i=1,2)


      write(6,100)      '   L2EVISIA    ', l2e_visiA_flag   
      write(6,2300)     '     THE2E     ', (the_2e(i),i=1,3)
      write(6,2200)     '     EMIN2E    ', (E_min_2e(i),i=1,2)
      write(6,100)      '   L2EVISIB    ', l2e_visiB_flag   
      write(6,2500)     '     THE2EB    ', (the_2eB(i),i=1,5)
      write(6,2400)     '     EMIN2EB   ', (E_min_2eB(i),i=1,4)
      write(6,100)      '   L3EVISI     ', l3e_visi_flag    
      write(6,2400)     '     THE3E     ', (the_3e(i),i=1,4)
      write(6,2100)     '     EMIN3E    ', E_min_3e
      write(6,100)      '   LPTVISI     ', lPt_visi_flag    
      write(6,2300)     '     THEMINPT  ', (the_min_pt(i),i=1,3)
      write(6,2300)     '     THEMAXPT  ', (the_max_pt(i),i=1,3)
      write(6,2300)     '     PTMINPT   ', (Pt_min_pt(i),i=1,3)
      write(6,100)      '   LSCATL      ', lscattL_flag     
      write(6,2200)     '     THESCATL  ', (theta_scattL(i),i=1,2)
      write(6,2200)     '     ESCATL    ', (E_scattL(i),    i=1,2)
      write(6,2200)     '     THEPRODL  ', ( theta_prodL(i),i=1,2)
      write(6,2200)     '     PPRODL    ', ( P_prodL(i),    i=1,2)
      write(6,*)   ' --------------------------- '
      write(6,*)   '=========================================='
      write(6,*)   '========== Parameters in Quasi-elastic =========='
      write(6,2100)     '   RK          ', rk
      write(6,2100)     '   ANBD        ', A_NBD
      write(6,2100)     '   BNBD        ', B_NBD
      write(6,2100)     '   CNBD        ', C_NBD
      write(6,2100)     '   PUQ         ', P_u
      write(6,2100)     '   PDQ         ', P_d
      write(6,2100)     '   PSQ         ', P_s
      write(6,2100)     '   PUD         ', P_ud_bar
      write(6,2100)     '   SUPPI0      ', Supp_pi0
      write(6,2100)     '   ASLOPE      ', A_slope
      write(6,2100)     '   BSLOPE      ', B_slope
      write(6,2100)     '   CSLOPE      ', C_slope
      write(6,2100)     '   REDS        ', red_slope_s
      write(6,*)   '================================================='
 100  format(A15, L10)
 130  format(A15, 3L5)
 1100 format(A15,   I10)
 1500 format(A15, 5(I10,1X))
 2100 format(A15,   G10.4)
 2200 format(A15, 2(G10.4,2X))
 2300 format(A15, 3(G10.4,2X))
 2400 format(A15, 4(G10.4,2X))
 2500 format(A15, 5(G10.4,2X))
* =============================================================
      return
      end
