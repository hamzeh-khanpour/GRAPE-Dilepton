*(graepia.h)
* ------------ # of processes in GRAPE ------------
      integer  num_prc_ela, num_prc_ela_qela
     &        ,num_prc_4f,  num_prc_4f_qq,   num_prc

       parameter ( num_prc_ela      =   3 )   !@  ela
       parameter ( num_prc_ela_qela =   6 )   !@  ela+qela
       parameter ( num_prc_4f       =  42 )   !@  ela+qela+DIS
       parameter ( num_prc_4f_qq    =  78 )   !@  ela+qela+DIS+res
       parameter ( num_prc          = 114 )   !@  ela+qela+DIS+res+dir
* -------------------------------------------------

* -------------- Event generation --------------
      integer      ISUB_start
       parameter ( ISUB_start = 400 )
* ----------------------------------------------

* ------------ GRAPE parameters ------------
      integer      num_gra_flg
       parameter ( num_gra_flg = 10 )
* ------------------------------------------

      integer*4          KF_Lbeam, ipol_Ebeam
      real*4             P_e_beam, P_p_beam,   Ebeam_pol(3)
      logical            Lpol_Ebeam
       common /GEP_BEAM/ KF_Lbeam, ipol_Ebeam
     &                  ,P_e_beam, P_p_beam,   Ebeam_pol
     &                  ,Lpol_Ebeam

      integer*4          process, lpair, qflv
     &                  ,Isr_flag, NN_ISR, ISR_scale
     &                  ,jgra_flag(num_gra_flg)
     &                  ,Ngroup,Nset, Istrf
     &                  ,merge, IQCD_scale
      real*4             Wmin,Wmax, Factor_ISR
       common /GEP_PROC/ process, lpair, qflv
     &                  ,Isr_flag, NN_ISR, ISR_scale, Factor_ISR
     &                  ,jgra_flag
     &                  ,Ngroup,Nset, Istrf
     &                  ,Wmin,Wmax
     &                  ,merge, IQCD_scale

      integer*4        jgra_sel
      real*4           weimj_QED, weimj_Z0, weimj_H, am_higgs,ag_higgs
       common /GEP_EW/ jgra_sel
     &                ,weimj_QED, weimj_Z0, weimj_H, am_higgs,ag_higgs

      integer*4          isym_34, ii34, Iscale_xq, Ireso56, Ireso456
     &                  ,Nregion, IcosP3,IcosP5, Neps_p
      real*4             thresh56, Rnn_MJ56, Mas_reso456,Wid_reso456
     &                  ,Rscale_Mn, Rnn_456
       common /GEP_KINE/ isym_34, ii34, Iscale_xq, Ireso56, Ireso456
     &                  ,Nregion, IcosP3,IcosP5, Neps_p
     &                  ,thresh56, Rnn_MJ56, Mas_reso456,Wid_reso456
     &                  ,Rscale_Mn, Rnn_456

      integer*4        num_it_grid, num_it_integ, num_call, mxtry,Ngen
      real*4           acc1_card, acc2_card
       common /GEP_BS/ num_it_grid, num_it_integ, num_call, mxtry,Ngen
     &                ,acc1_card, acc2_card

      integer*4        PS_isr, PS_fsr, PS_bra, PS_sup, PY_decay
     &                ,Ipri_Pt
      real*4           E_min_isr
       common /GEP_PY/ PS_isr, PS_fsr, PS_bra, PS_sup, PY_decay
     &                ,Ipri_Pt
     &                ,E_min_isr

      logical           lHOT_flag, l2e_visiA_flag, l2e_visiB_flag
     &                 ,l3e_visi_flag,  lPt_visi_flag, lscattL_flag
     &                 ,lPtMAX_cut, lveto, l3D_cut,l3E_cut,l3F_cut
      real*4            
     &                  Q2p_cut(2)
     &                 ,theta_min(4),theta_max(4)
     &                 ,E_min(4),E_max(4), P_min(4),P_max(4)
     &               ,Pt_min(4),Pt_max(4), PtMAX_cut(2),the_PtMAX_cut(2)
     &                 ,E_scattL(2), theta_scattL(2)
     &                   ,P_prodL(2), theta_prodL(2)
     &                 ,Mass56_cut(4), MassELL_cut(2), MassQLL_cut(4)
     &                 ,the_visi_min(4), the_visi_max(4)
     &                 ,ene_visi_min(4), ene_visi_max(4)
     &                 ,Pt_visi_min(4),  Pt_visi_max(4)

     &                 ,the_2e(3)     !@  min, thre, max
     &                 ,E_min_2e(2)   !@  E_minf_2e, E_minb_2e
     &                 ,the_3e(4)     !@  the_fmin_3e,the_fmax_3e
*!@                                      ,the_bmin_3e,the_bmax_3e
     &                 ,E_min_3e
     &                 ,the_min_pt(3), the_max_pt(3), Pt_min_pt(3)

     &                 ,the_2eB(5)
     &                 ,E_min_2eB(4)
     &                 ,u_cut(2), W_cut(2)

     &                 ,the_veto_min(4), the_veto_max(4)
     &                 ,ene_veto_min(4), ene_veto_max(4)
     &                 ,Pt_veto_min(4),  Pt_veto_max(4)

     &      ,x_Range_ME(4), y_Range_ME(4), Q2_Range_ME(4), W_Range_ME(4)
     &      ,x_Range_OB(4), y_Range_OB(4), Q2_Range_OB(4), W_Range_OB(4)

     &                 ,the_3D_min(3), the_3D_max(3)
     &                 ,  E_3D_min(3),   E_3D_max(3)
     &                 ,  P_3D_min(3),   P_3D_max(3)
     &                 , Pt_3D_min(3),  Pt_3D_max(3)
     &                 , Mass12_3D(2),  Q2_3D_ME(2),Q2_3D_OB(2)
     &                 , W_3D_ME(2),W_3D_OB(2)
     &                 , X_3D_ME(2),X_3D_OB(2), Y_3D_ME(2),Y_3D_OB(2)

     &                 ,the_3E_min(3), the_3E_max(3)
     &                 ,  E_3E_min(3),   E_3E_max(3)
     &                 ,  P_3E_min(3),   P_3E_max(3)
     &                 , Pt_3E_min(3),  Pt_3E_max(3)
     &                 , Mass12_3E(2),  Q2_3E_ME(2),Q2_3E_OB(2)
     &                 , W_3E_ME(2),W_3E_OB(2)
     &                 , X_3E_ME(2),X_3E_OB(2), Y_3E_ME(2),Y_3E_OB(2)

     &                 ,the_3F_min(3), the_3F_max(3)
     &                 ,  E_3F_min(3),   E_3F_max(3)
     &                 ,  P_3F_min(3),   P_3F_max(3)
     &                 , Pt_3F_min(3),  Pt_3F_max(3)
     &                 , Mass12_3F(2),  Q2_3F_ME(2),Q2_3F_OB(2)
     &                 , W_3F_ME(2),W_3F_OB(2)
     &                 , X_3F_ME(2),X_3F_OB(2), Y_3F_ME(2),Y_3F_OB(2)

      integer*4         Ivisi, Iveto
       common /GEP_CUT/ 
     &                  Q2p_cut
     &                 ,theta_min,   theta_max
     &                 ,E_min,   E_max,    P_min,   P_max
     &                 ,Pt_min,   Pt_max
     &                 ,lPtMAX_cut,  PtMAX_cut,  the_PtMAX_cut
     &                 ,E_scattL, theta_scattL, P_prodL, theta_prodL
     &                 ,Mass56_cut, MassELL_cut, MassQLL_cut
     &                 ,the_visi_min,  the_visi_max
     &                 ,ene_visi_min,  ene_visi_max
     &                 ,Pt_visi_min,   Pt_visi_max
     &                 ,Ivisi
     &                 ,the_2e        !@  min, thre, max
     &                 ,E_min_2e      !@  E_minf_2e, E_minb_2e
     &                 ,the_3e        !@  the_fmin_3e,the_fmax_3e
*!@                                       ,the_bmin_3e,the_bmax_3e
       common /GEP_CUT/ E_min_3e
     &                 ,the_min_pt, the_max_pt, Pt_min_pt

     &                 ,the_2eB, E_min_2eB

     &                 ,lHOT_flag, l2e_visiA_flag, l2e_visiB_flag
     &                 ,l3e_visi_flag,  lPt_visi_flag, lscattL_flag

     &                 ,u_cut, W_cut

     &                 ,lveto, Iveto
     &                 ,the_veto_min, the_veto_max
     &                 ,ene_veto_min, ene_veto_max
     &                 ,Pt_veto_min,  Pt_veto_max

     &      ,x_Range_ME, y_Range_ME, Q2_Range_ME, W_Range_ME
     &      ,x_Range_OB, y_Range_OB, Q2_Range_OB, W_Range_OB

       common /GEP_CUT/ l3D_cut
     &                 ,the_3D_min, the_3D_max
     &                 ,  E_3D_min,   E_3D_max
     &                 ,  P_3D_min,   P_3D_max
     &                 , Pt_3D_min,  Pt_3D_max
     &                 , Mass12_3D,  Q2_3D_ME, Q2_3D_OB
     &                 , W_3D_ME,W_3D_OB
     &                 , X_3D_ME,X_3D_OB, Y_3D_ME,Y_3D_OB
     &                 ,l3E_cut
     &                 ,the_3E_min, the_3E_max
     &                 ,  E_3E_min,   E_3E_max
     &                 ,  P_3E_min,   P_3E_max
     &                 , Pt_3E_min,  Pt_3E_max
     &                 , Mass12_3E,  Q2_3E_ME, Q2_3E_OB
     &                 , W_3E_ME,W_3E_OB
     &                 , X_3E_ME,X_3E_OB, Y_3E_ME,Y_3E_OB
     &                 ,l3F_cut
     &                 ,the_3F_min, the_3F_max
     &                 ,  E_3F_min,   E_3F_max
     &                 ,  P_3F_min,   P_3F_max
     &                 , Pt_3F_min,  Pt_3F_max
     &                 , Mass12_3F,  Q2_3F_ME, Q2_3F_OB
     &                 , W_3F_ME,W_3F_OB
     &                 , X_3F_ME,X_3F_OB, Y_3F_ME,Y_3F_OB


      logical           LIST_flag, NTPYT_flag, ASC_flag, NTVEC_flag
     &                 ,Wrt_cards, LRND_flag
      integer*4         icount, print_flag, Nmod, Nlist, Frame_amp
     &                 ,Qela_decay
       common /GEP_MISC/icount, print_flag, Nmod, Nlist, Frame_amp
     &                 ,Qela_decay
     &                 ,LIST_flag, NTPYT_flag, ASC_flag, NTVEC_flag
     &                 ,Wrt_cards, LRND_flag

      logical           ROT_kinem_flag, ROT_pyrand_flag
       common /GEP_ROT/ ROT_kinem_flag, ROT_pyrand_flag

      integer           Inn
      double precision  Dnn, Diff_nn
       common /GEP_NN/  Dnn, Diff_nn, Inn

*---------- User-defined functions ----------
      integer    Ibtest
       external  Ibtest
*--------------------------------------------

*---------- For quasi-elastic process ----------
      integer     mult_limit          !@  Max. # of produced particles
       parameter (mult_limit = 400)     !@     (neutral + charged)

      real*4
     &         rk, A_NBD,B_NBD,C_NBD   !@  Input for NBD
     &        ,P_u,P_d,P_s,P_ud_bar, Supp_pi0 !@ Input for flavor selection
     &        ,A_slope, B_slope, C_slope, red_slope_s !@ Input for Pt slope
      common /GEP_QELA/
     &         rk, A_NBD,B_NBD,C_NBD   !@  Input for NBD
     &        ,P_u,P_d,P_s,P_ud_bar, Supp_pi0 !@ Input for flavor selection
     &        ,A_slope, B_slope, C_slope, red_slope_s !@ Input for Pt slope
*-----------------------------------------------

*-------------- Vector Meson production -------------
      integer*4  Model_VM
      real*4     t_slope_VM
      logical    lee_int_VM, helicity_VM(-1:+1)
       common /GEP_VM/  Model_VM, t_slope_VM, lee_int_VM
     &                 ,helicity_VM

*>>> J/psi(psi)
      integer      num_jptype,   num_jpgra
       parameter ( num_jptype=6, num_jpgra=2 )
      logical  lJPprod, lJPtype(num_jptype), lJPgra(num_jpgra)
      real*4   bcorr_jp(num_jptype), weimj_jp(num_jptype)
     &        ,phaseP_jp(num_jptype), phaseG_jp(num_jptype)
       common /GEP_JP/
     &     lJPprod, lJPtype, lJPgra
     &    ,bcorr_jp, weimj_jp, phaseP_jp, phaseG_jp

*>>> Upsilon
      integer      num_yytype,   num_yygra
       parameter ( num_yytype=6, num_yygra=2 )
      logical  lYYprod, lYYtype(num_yytype), lYYgra(num_yygra)
      real*4   bcorr_yy(num_yytype), weimj_yy(num_yytype)
     &        ,phaseP_yy(num_yytype), phaseG_yy(num_yytype)
       common /GEP_YY/
     &     lYYprod, lYYtype, lYYgra
     &    ,bcorr_yy, weimj_yy, phaseP_yy, phaseG_yy
*----------------------------------------------------

*---------- Stuff for Multi-Jacobian method ----------
      integer      max_reso_mj
       parameter ( max_reso_mj=20 )
      integer  num_reso_mj
      double precision
     &         mass_mj(max_reso_mj), mass2_mj(max_reso_mj)
     &        ,width_mj(max_reso_mj)
     &        ,weight_mj(0:max_reso_mj), weisel_mj(0:max_reso_mj)

       common /GEP_MJ/ mass_mj, mass2_mj, width_mj
     &                ,weight_mj, weisel_mj
     &                ,num_reso_mj
*-----------------------------------------------------

*----- Stuff for calculations in various frames -----
       double precision  PE_str(4,12)
        common /GEP_FRM/ PE_str
*----------------------------------------------------
