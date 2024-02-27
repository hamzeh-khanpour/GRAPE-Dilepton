* This kinematics is based on the one by Y.Kurihara.
      SUBROUTINE KINEM_4f(NEXTRN,X_arg,PE,PP,YACOB,NREG,IREG,JUMP)
      IMPLICIT REAL* 8(A-H,O-Z)
      INTEGER NEXTRN
      PARAMETER ( MXDIM = 50 )
      REAL*8  X_arg(MXDIM),X(MXDIM)
      REAL*8  PE(4,NEXTRN), PP(NEXTRN,NEXTRN)
      REAL*8  YACOB
      INTEGER NREG, IREG
      INTEGER JUMP
      INCLUDE 'inclk.h'         
      INCLUDE 'inclp_4f.h'      
      complex*16 aprop
      real*8 Dk(2)
*--------------------------------------------------------------------
      COMMON/KINEM1/S(3),W(3),FACT
      COMMON/CUT001_4f/COSCUT(2,4),ENGYCT(2,4),AMASCT(2,6),ARESNS(2,4)
     .,opncut,swapm2
      common/atest/ttt1,uuu1,ttt2,uuu2,tttt,uuuu,sss2,sss1,q56,q12,q22
     .,d1,g1,t1,u1,t2,u2
      common/kmcntl_4f/iresns(4),icos3,icosq3,icos5,isr,iswap,ident
     &                ,iphi6,ieeee,i34,itag,isym
      common/TAISR/yy2,y2max
*--------------------------------------------------------------------
      COMMON / BPARM1 / XL(MXDIM),XU(MXDIM),NDIM,NWILD,
     &                 IG(MXDIM),NCALL
* -------------------- PDFLIB stuff --------------------
      integer           IFLPRT
       common /W50510/  IFLPRT
      double precision  xmin,xmax, Q2min,Q2max
       common /W50513/  xmin,xmax, Q2min,Q2max
      logical           FIRST
       common /W50516/  FIRST
                                                   
      double precision  UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GL
       common /GEP_PDF/ UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GL
*-------------------------------------------------------
      include './inc/graepia.h'
      integer        jproc      
      common /amjprc/jproc
*-------------------------------------------------------
      logical           lvisi,lvisi3,lvisi4,lvisi5,lvisi6,lvisi_PtMAX(4)
     &                 ,l2e_visi4, l2e_visi5, l2e_visi6, l2e_visiA
     &                                                  ,l2e_visiB
     &                 ,l3e_visi4, l3e_visi5, l3e_visi6, l3e_visi
     &                 ,lPt_visi, lscattL_visi, lprodL_visi(2),lemu_visi

      logical                  !  Cut1-3  Leptons
     &                  l3D_flag(    3,      3  )
     &                 ,l3E_flag(    3,      3  )
     &                 ,l3F_flag(    3,      3  )
      double precision  aM3_12, aM3_23, aM3_31
     &                 ,Q2_1_ME, Q2_2_ME, Q2_3_ME
     &                 ,Q2_1_OB, Q2_2_OB, Q2_3_OB
     &                 ,W_1_ME, W_2_ME, W_3_ME
     &                 ,W_1_OB, W_2_OB, W_3_OB
     &                 ,X_1_ME, X_2_ME, X_3_ME
     &                 ,X_1_OB, X_2_OB, X_3_OB
     &                 ,Y_1_ME, Y_2_ME, Y_3_ME
     &                 ,Y_1_OB, Y_2_OB, Y_3_OB

      integer           i,ii,j, Ivisi1,Ivisi2, num_visi, visi_num(3)
      double precision   momt(4), theta(4), the_d(4)
     &                 ,Et3,   Et4,   Et5,   Et6
     &                 ,Pt3,   Pt4,   Pt5,   Pt6,  Pt(4),Pt_2(4)
     &                 ,Q2_QCD, PtMAX, deg, mass_ee(2)
     &                 ,u_val, W_val,W_val2, Pt_max_tmp(2)
     &                 ,mass_ell,mass_qll(2), PE_g(4), PE_gp(4)
      common /QCD_SCALE/ Q2_QCD
      double precision  x_el_ME, y_el_ME, Q2_el_ME, W_el_ME, S_ME
     &                 ,x_el_OB, y_el_OB, Q2_el_OB, W_el_OB, S_OB
      common /PARM_KINEM4F/
     &                  x_el_ME, y_el_ME, Q2_el_ME, W_el_ME, S_ME
     &                 ,x_el_OB, y_el_OB, Q2_el_OB, W_el_OB, S_OB
*-------------------------------------------------------
      double precision  P1_lab,P2_lab, E1_lab,E2_lab  
     &                 ,Pcms_lab(3),Ecms_lab(3)
     &                 ,GAMMAcms_lab(3),BETGAMcms_lab(3)
     &                 ,vec_isr(4)
       common /GEP_LAB/ P1_lab,P2_lab, E1_lab,E2_lab
     &                 ,Pcms_lab,Ecms_lab
     &                 ,GAMMAcms_lab,   BETGAMcms_lab
     &                 ,vec_isr
      double precision   Pq_lab,Eq_lab, YR,yacob_x, PP_gp_LAB(4)
      double precision  hradiO   
      external          hradiO   
      logical  lee_int
      icount = icount + 1
      if (icount.EQ.1) then
         FIRST = .true.
       else
         FIRST = .false.
      endif
      if ( (lpair .EQ. 1) .AND.      
     &     (         (jgra_sel .GT. 1)   
     &          .OR. (lJPprod .and. lee_int_VM)
     &          .OR. (lYYprod .and. lee_int_VM)
     &     )
     &   ) then
         NREG = Nregion        
         lee_int = .true.
       else
         NREG = 1              
         lee_int = .false.
      endif
      jump   =   0
      YACOB  =   0.D0
      AJACOB=    1.D0
      do i = 1, NDIM
        X(i) = X_arg(i)
      enddo
      if (NDIM .GT. 7) then
        call TAswap(X(4),X(NDIM-1-ipol_Ebeam))
        call TAswap(X(6),X(NDIM-ipol_Ebeam))
      else
        call TAswap(X(4),X(NDIM-ipol_Ebeam))
      endif
********** Determination of the Mass of hadronic system *********
      if (process.EQ.2) then   
        DD = Wmax - Wmin
        if (Rscale_Mn .LE. 1.) then
          AMASS1(3) = dble(Wmin) + DD*X(8)
          AJACOB = AJACOB * 2.D0*AMASS1(3)*DD /(2.D0*PI)
        else
          AMASS1(3) = dble(Wmin) + DD*X(8)**Rscale_Mn
          AJACOB = AJACOB * 2.D0*AMASS1(3)
     &                    * DD*Rscale_Mn*X(8)**(Rscale_Mn-1.)
     &                     /(2.D0*PI)
        endif
        AMASS2(3) = AMASS1(3) * AMASS1(3)
      endif
****************************************************************
      if ( process .EQ. 3 ) then
       if      (Iscale_xq .EQ. 1) then
          YR      = xmax/xmin            
          xq      = xmin * YR**X(8)      
          yacob_x = xq * log(YR)         
        elseif (Iscale_xq .EQ. 2) then
          YR      = xmax-xmin            
          xq      = xmin + YR*X(8)       
          yacob_x = YR                   
        else
          write(6,*) '!!!Error in kinem_4f!!!'
          write(6,*) '  ---> unknown scale for xq(=',Iscale_xq,' )'
          write(6,*) '  ---> Good-bye!'
          STOP
       endif
        Pq_lab = P1_lab * xq
        Eq_lab = sqrt(Pq_lab**2 + AMASS2(1))
       S(2) = AMASS2(1)+AMASS2(2) + 2.D0*(Eq_lab*E2_lab + Pq_lab*P2_lab)
        W(2) = sqrt(  max( S(2), 0.D0 )  )
       Pcms_lab(2) = Pq_lab - P2_lab
       Ecms_lab(2) = Eq_lab + E2_lab
        GAMMAcms_lab(2)  = Ecms_lab(2)/W(2)
        BETGAMcms_lab(2) = Pcms_lab(2)/W(2)
       totmas=amass1(3)+amass1(4)+amass1(5)+amass1(6)
       if (W(2) .EQ. 0) then
	  write(6,*)
     .  ' CM energy is less than the sum of final particle masses.'
          write(6,*) '   (CM energy) =', W(2)
          write(6,*) '   (the sum of final particle masses) =', totmas
          write(6,*) ' Something is wrong in kinem. ---> Good-bye!'
          STOP
        elseif (W(2).LT.totmas) then 
          jump = 1                   
          icount = icount - 1
          GOTO 999
       endif
C<TA>       if (icount.EQ.1) then
C<TA>         write(6,*) '********** Information (in Lab. frame) **********'
C<TA>         write(6,*) '                              (in unit of GeV)'
C<TA>         write(6,*) '  P of electrons    =',real(P2_lab)
C<TA>         write(6,*) '  P of quark        =',real(Pq_lab)
C<TA>         write(6,*) '  Mass of quark     =',real(AMASS1(1))
C<TA>         write(6,*) '  Bjorken-x         =',real(xq)
C<TA>         write(6,*) '  sqrt(S)           =',real(W(2))
C<TA>         write(6,*) '  P of CMS          =',real(Pcms_lab(2))
C<TA>         write(6,*) '  E of CMS          =',real(Ecms_lab(2))
C<TA>         write(6,*) '  gamma of CMS      =',real(GAMMAcms_lab(2))
C<TA>         write(6,*) '  beta*gamma of CMS =',real(BETGAMcms_lab(2))
C<TA>         write(6,*) '*************************************************'
C<TA>       endif
      endif   
      if (ireg.EQ.1) then
       icosq3=2
       call k6nem0(NEXTRN,X,PE,PP,YACOB,NREG,IREG,JUMP)
       if(jump.eq.1) goto 999
       vn10 = TTT1   
       vn20 = TTT2
       vn22 = SSS1
       vn28 = UUU1
       vn30 = Q56
       vn26 = UUU2
       vn52 = UUUU
       vn58 =  - 2.0d0*pp(2,6) + amass2(2) + amass2(6)
       vn56 =  + 2.0d0*pp(3,4) + 2.0d0*pp(3,5) + 2.0d0*pp(4,5)
     &       + amass2(3) + amass2(4) + amass2(5)
       vn50 =  - 2.0d0*pp(1,4) - 2.0d0*pp(1,5) + 2.0d0*pp(4,5)
     &       + amass2(1) + amass2(4) + amass2(5)
       vn48 =  + 2.0d0*pp(4,5) + amass2(4) + amass2(5)
      else if(ireg.eq.2) then
       icosq3=2
       call k6exch(4,6,pe,pp)
       call k6nem0(NEXTRN,X,PE,PP,YACOB,NREG,IREG,JUMP)
       if(jump.eq.1) then
        call k6exch(4,6,pe,pp)
        goto 999
       else
        call k6exch(6,4,pe,pp)
       end if
       vn10 = TTT1   
       vn58 = TTT2
       vn56 = SSS1
       vn50 = UUU1
       vn48 = Q56
       vn26 = UUUU
       vn52 = UUU2
       vn20 =  - 2.0d0*pp(2,4) + amass2(2) + amass2(4)
       vn22 =  + 2.0d0*pp(1,2) - 2.0d0*pp(1,4) - 2.0d0*pp(2,4)
     &       + amass2(1) + amass2(2) + amass2(4)
       vn28 =  - 2.0d0*pp(2,3) - 2.0d0*pp(2,4) + 2.0d0*pp(3,4)
     &       + amass2(2) + amass2(3) + amass2(4)
       vn30 =  + 2.0d0*pp(5,6) + amass2(5) + amass2(6)
      END if
      if (lee_int) then
         Pt4 = sqrt( PE(1,4)*PE(1,4) + PE(2,4)*PE(2,4) )
         Pt6 = sqrt( PE(1,6)*PE(1,6) + PE(2,6)*PE(2,6) )
         if (Pt4 .GT. Pt6)   GOTO 999   ! Selecting 'Pt4<Pt6' region
         aprop=1
         vm = vn20                              
         call smprpd(aprop,vm,ama**2,0.0d0)
         D20=abs(aprop)**2
         aprop=1
         vm = vn58                              
         call smprpd(aprop,vm,ama**2,0.0d0)
         D58=abs(aprop)**2
         aprop=1
         vm = vn30                              
         call smprpd(aprop,vm,ama**2,0.0d0)
         D30=abs(aprop)**2
         aprop=1
         vm = vn48                              
         call smprpd(aprop,vm,ama**2,0.0d0)
         D48=abs(aprop)**2
         D30_VM = 1D0
         D48_VM = 1D0
         if (num_reso_mj .GE. 1) then
           do i = 1, num_reso_mj
             aprop=1
             vm = vn30                              
             call smprpd(aprop,vm,mass2_mj(i),mass_mj(i)*width_mj(i))
             D30_VM = D30_VM *abs(aprop)**2
             aprop=1
             vm = vn48                              
             call smprpd(aprop,vm,mass2_mj(i),mass_mj(i)*width_mj(i))
             D48_VM = D48_VM *abs(aprop)**2
           enddo
         endif
         D26 = 1D0
         D52 = 1D0
         if (.false.
     +       .OR. jgra_sel .EQ. 3
     +       .OR. jgra_sel .EQ. 4
     +       .OR. jgra_sel .EQ. 5
     +       .OR. jgra_sel .EQ.13
     +       .OR. jgra_sel .EQ.14
     +       .OR. (lJPprod.and.lJPgra(2))
     +       .OR. (lYYprod.and.lYYgra(2))
     +      ) then
           aprop=1
           vm = vn26                              
           call smprpd(aprop,vm,amass2(2),0.0d0)
           D26 = D26*abs(aprop)**2
           aprop=1
           vm = vn52                              
           call smprpd(aprop,vm,amass2(2),0.0d0)
           D52 = D52*abs(aprop)**2
         endif
         D28 = 1D0
         D50 = 1D0
         if (process .EQ. 3) then
           if (.false.
     +         .OR. jgra_sel .EQ. 3
     +         .OR. jgra_sel .EQ. 4
     +         .OR. jgra_sel .EQ. 5
     +         .OR. jgra_sel .EQ.13
     +         .OR. jgra_sel .EQ.14
     +         .OR.(jgra_sel .EQ.15 .and. qflv.GT.6)
     +        ) then
             aprop=1
             vm = vn28                              
             call smprpd(aprop,vm,amass2(1),0.0d0)
             D28 = D28*abs(aprop)**2
             aprop=1
             vm = vn50                              
             call smprpd(aprop,vm,amass2(1),0.0d0)
             D50 = D50*abs(aprop)**2
           endif
         endif
         D22 = 1D0
         D56 = 1D0
         if (process .EQ. 3) then
           if (.false.
     +         .OR. jgra_sel .EQ. 3
     +         .OR. jgra_sel .EQ. 4
     +         .OR. jgra_sel .EQ. 5
     +         .OR. jgra_sel .EQ.13
     +         .OR. jgra_sel .EQ.14
     +         .OR.(jgra_sel .EQ.15 .and. qflv.GT.6)
     +        ) then
             aprop=1
             vm = vn22
             call smprpd(aprop,vm,amass2(1),0.0d0)
             D22 = D22*abs(aprop)**2
             aprop=1
            vm = vn56
             call smprpd(aprop,vm,amass2(1),0.0d0)
             D56 = D56*abs(aprop)**2
           endif
         endif
*************************************************************
         Dk(1) = D58 * D48 * D52 * D50 * D56 * D48_VM
         Dk(2) = D20 * D30 * D26 * D28 * D22 * D30_VM
         DDD = Dk(1)+Dk(2)                  
         yacob = yacob * Dk(ireg)/DDD  
      endif
      if ( process .EQ. 3 ) then
        if (IQCD_scale .EQ. 1) then   
           Q2_QCD = abs(vn10)
         else
           write(6,*) '!!!Error in kinem_4f!!!'
           write(6,*) '  ---> unknown choice of QCD scale('
     &                                                ,IQCD_scale,' )'
           write(6,*) '  ---> Good-bye!'
           STOP
        endif
        if (Q2_QCD .LT. Q2min) then     
           GOTO 999
        endif
        call STRUCTM(xq, sqrt(Q2_QCD)      
     +              ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GL)
        if      ((qflv.EQ. 1).and.(merge .EQ. 0)) then
           YACOB = YACOB * (UPV + USEA)/xq * yacob_x      ! u    quark
         elseif ((qflv.EQ. 2).and.(merge .EQ. 0)) then
           YACOB = YACOB *        USEA /xq * yacob_x      ! u^   quark
         elseif ((qflv.EQ. 3).and.(merge .EQ. 0)) then
           YACOB = YACOB * (DNV + DSEA)/xq * yacob_x      ! d    quark
         elseif ((qflv.EQ. 4).and.(merge .EQ. 0)) then
           YACOB = YACOB *        DSEA /xq * yacob_x      ! d^   qurak
         elseif ((qflv.LE. 6).and.(merge .EQ. 0)) then
           YACOB = YACOB *        STR  /xq * yacob_x      ! s,s^ quark
         elseif ((qflv.LE. 8).and.(merge .EQ. 0)) then
           YACOB = YACOB *        CHM  /xq * yacob_x      ! c,c^ quark
         elseif ((qflv.LE.10).and.(merge .EQ. 0)) then
           YACOB = YACOB *        BOT  /xq * yacob_x      ! b,b^ quark
         elseif ((qflv.LE.12).and.(merge .EQ. 0)) then
           YACOB = YACOB *        TOP  /xq * yacob_x      ! t,t^ quark
         elseif ((qflv.EQ. 1).and.(merge .EQ. 1234)
     +                       .and.(jgra_sel .LE. 2)) then    
           YACOB = YACOB * (UPV + 2.D0*USEA + DNV/4.D0 + DSEA/2.D0)
     +               /xq * yacob_x
         elseif ((qflv.EQ. 1).and.(merge .EQ. 123456)
     +                       .and.(jgra_sel .LE. 2)) then    
           YACOB = YACOB * (UPV + 2.D0*USEA + DNV/4.D0 + DSEA/2.D0
     +                       + STR/2.D0)
     +               /xq * yacob_x
         elseif ((qflv.EQ. 1).and.(merge .EQ. 12345678)
     +                       .and.(jgra_sel .LE. 2)) then    
           YACOB = YACOB * (UPV + 2.D0*USEA + DNV/4.D0 + DSEA/2.D0
     +                       + STR/2.D0 + 2.D0*CHM)
     +               /xq * yacob_x
         elseif ((qflv.EQ. 1).and.(merge .EQ. 1234567890)
     +                       .and.(jgra_sel .LE. 2)) then    
           YACOB = YACOB * (UPV + 2.D0*USEA + DNV/4.D0 + DSEA/2.D0
     +                       + STR/2.D0 + 2.D0*CHM + BOT/2.D0)
     +               /xq * yacob_x
         elseif ((qflv.EQ. 1).and.(merge .EQ. 17)   ! u+c
     +                       .and.(jgra_sel .GE. 1)) then    
           YACOB = YACOB * (UPV + USEA + CHM) / xq * yacob_x
         elseif ((qflv.EQ. 2).and.(merge .EQ. 28)   ! u^+c^
     +                       .and.(jgra_sel .GE. 1)) then    
           YACOB = YACOB * (      USEA + CHM) / xq * yacob_x
         elseif ((qflv.EQ. 3).and.(merge .EQ. 35)   ! d+s
     +                       .and.(jgra_sel .GE. 1)) then    
           YACOB = YACOB * (DNV + DSEA + STR) / xq * yacob_x
         elseif ((qflv.EQ. 3).and.(merge .EQ. 359)  ! d+s+b
     +                       .and.(jgra_sel .GE. 1)) then    
           YACOB = YACOB * (DNV + DSEA + STR + BOT) / xq * yacob_x
         elseif ((qflv.EQ. 4).and.(merge .EQ. 46)   ! d^+s^
     +                       .and.(jgra_sel .GE. 1)) then    
           YACOB = YACOB * (DSEA + STR) / xq * yacob_x
         elseif ((qflv.EQ. 4).and.(merge .EQ. 460)  ! d^+s^+b^
     +                       .and.(jgra_sel .GE. 1)) then    
           YACOB = YACOB * (DSEA + STR + BOT) / xq * yacob_x
         else
           write(6,*) '!!!Error in kinem_4f!!!'
           write(6,*) '  ---> unknown qflv/merge(=',qflv,'/',merge,' )'
           write(6,*) '  ---> Good-bye!'
           STOP
        endif
      endif
* ------------ Boost : CMS ---> Lab. frame ------------
       do k=1,NEXTRN
         Pz_cms = PE(3,k)
         E_cms  = PE(4,k)
           PE(3,k) = GAMMAcms_lab(3)  * Pz_cms
     &             + BETGAMcms_lab(3) * E_cms
           PE(4,k) = GAMMAcms_lab(3)  * E_cms
     &             + BETGAMcms_lab(3) * Pz_cms
       enddo
* -----------------------------------------------------
* ------------- Cuts in Lab. frame -------------
       if (lee_int) then
         if (abs(vn20) .LT. abs(vn58)) then
            Iscatt = 4
            Q2_el_ME = max(-vn20,0D0)
         else
            Iscatt = 6
            Q2_el_ME = max(-vn58,0D0)
         endif
       else
            Iscatt = 4
            Q2_el_ME = max(-vn20,0D0)
       endif
*>>> Q2 cut
      Q2_el_OB = 2D0*AMASS2(2) -2D0*( PE(4,Iscatt)*E2_lab
C    &                               -PE(1,Iscatt)*0D0
C    &                               -PE(2,Iscatt)*0D0
     &                               +PE(3,Iscatt)*P2_lab
     &                              )
      Q2_el_OB = max(-Q2_el_OB,0D0)
      if (  .false.
     &   .OR. (Q2_Range_OB(1) .GE. 0)
     &   .OR. (x_Range_OB(1)  .GE. 0)
     &   .OR. (y_Range_OB(1)  .GE. 0)
     &   ) then
         if (Q2_el_OB.LT.Q2_Range_OB(1) .OR. Q2_el_OB.GT.Q2_Range_OB(2))
     &        GOTO 999
      endif
      if (  .false.
     &   .OR. (Q2_Range_ME(1) .GE. 0)
     &   .OR. (x_Range_ME(1)  .GE. 0)
     &   .OR. (y_Range_ME(1)  .GE. 0)
     &   ) then
         if (Q2_el_ME.LT.Q2_Range_ME(1) .OR. Q2_el_ME.GT.Q2_Range_ME(2))
     &        GOTO 999
      endif
*>>> W cut
      do i = 1, 4
        PE_g(i) = PE(i,2) - PE(i,Iscatt)
      enddo
      PE_gp(1) = PE_g(1) + 0D0
      PE_gp(2) = PE_g(2) + 0D0
      PE_gp(3) = PE_g(3) + P1_lab
      PE_gp(4) = PE_g(4) + E1_lab
      W_el_ME = sqrt( max( 0D0,
     &                  PE_gp(4)*PE_gp(4)
     &                 -PE_gp(1)*PE_gp(1)
     &                 -PE_gp(2)*PE_gp(2)
     &                 -PE_gp(3)*PE_gp(3) ) )
C      if (  .false.
C     &   .OR. (W_Range_ME(1) .GE. 0)
C     &   .OR. (x_Range_ME(1) .GE. 0)
C     &   .OR. (y_Range_ME(1) .GE. 0)
C     &   ) then
         if (W_el_ME.LT.W_Range_ME(1) .OR. W_el_ME.GT.W_Range_ME(2))
     &        GOTO 999
C      endif
      PE_g(1) =     0D0 - PE(1,Iscatt)
      PE_g(2) =     0D0 - PE(2,Iscatt)
      PE_g(3) = -P2_lab - PE(3,Iscatt)
      PE_g(4) =  E2_lab - PE(4,Iscatt)
      PE_gp(1) = PE_g(1) + 0D0
      PE_gp(2) = PE_g(2) + 0D0
      PE_gp(3) = PE_g(3) + P1_lab
      PE_gp(4) = PE_g(4) + E1_lab
      W_el_OB = sqrt( max( 0D0,
     &                  PE_gp(4)*PE_gp(4)
     &                 -PE_gp(1)*PE_gp(1)
     &                 -PE_gp(2)*PE_gp(2)
     &                 -PE_gp(3)*PE_gp(3) ) )
C      if (  .false.
C     &   .OR. (W_Range_OB(1) .GE. 0)
C     &   .OR. (x_Range_OB(1) .GE. 0)
C     &   .OR. (y_Range_OB(1) .GE. 0)
C     &   ) then
         if (W_el_OB.LT.W_Range_OB(1) .OR. W_el_OB.GT.W_Range_OB(2))
     &        GOTO 999
C      endif
*>>> x cut
      if (  .false.
     &   .OR. (x_Range_ME(1) .GE. 0)
C    &   .OR. (y_Range_ME(1) .GE. 0)
     &   ) then
         x_el_ME = Q2_el_ME /(Q2_el_ME +W_el_ME*W_el_ME -amp2)
         if (x_el_ME.LT.x_Range_ME(1) .OR. x_el_ME.GT.x_Range_ME(2))
     &        GOTO 999
      endif
      if (  .false.
     &   .OR. (x_Range_OB(1) .GE. 0)
C    &   .OR. (y_Range_OB(1) .GE. 0)
     &   ) then
         x_el_OB = Q2_el_OB /(Q2_el_OB +W_el_OB*W_el_OB -amp2)
         if (x_el_OB.LT.x_Range_OB(1) .OR. x_el_OB.GT.x_Range_OB(2))
     &        GOTO 999
      endif
*>>> y cut
      S_ME = (PE(4,2) + E1_lab)**2
     &      -(PE(1,2) +    0D0)**2
     &      -(PE(2,2) +    0D0)**2
     &      -(PE(3,2) + P1_lab)**2
      if (y_Range_ME(1) .GE. 0) then
         y_el_ME = (Q2_el_ME +W_el_ME*W_el_ME -amp2) /(S_ME-amp2)
         if (y_el_ME.LT.y_Range_ME(1) .OR. y_el_ME.GT.y_Range_ME(2))
     &      GOTO 999
      endif
      if (y_Range_OB(1) .GE. 0) then
         y_el_OB = (Q2_el_OB +W_el_OB*W_el_OB -amp2) /(S(1)-amp2)
         if (y_el_OB.LT.y_Range_OB(1) .OR. y_el_OB.GT.y_Range_OB(2))
     &      GOTO 999
      endif
      if (process .EQ. 3) then   
        if (lpair .EQ. 1) then
           u_val = min( abs(vn28), abs(vn50) )
         else
           u_val = abs(vn28)   
        endif
        if ((u_val .LT. u_cut(1)).or.(u_val .GT. u_cut(2))) GOTO 999
      endif
      if ( (abs(vn10) .LT. Q2p_cut(1)).OR.(abs(vn10) .GT. Q2p_cut(2)) )
     &     GOTO 999
      if (process .EQ. 3) then
        W_val2 = amp2 + amass2(2) + 2.D0*(E1_lab*pe(4,2)-P1_lab*pe(3,2))
     &         + amass2(4) + amass2(5) + amass2(6)   
     &            + (vn48-amass2(4)-amass2(5))
     &            + (vn30-amass2(5)-amass2(6))
     &            + 2.D0*PP(4,6)
     &         -2.D0*( E1_lab*(pe(4,4)+pe(4,5)+pe(4,6))   
     &                    - P1_lab*(pe(3,4)+pe(3,5)+pe(3,6)) 
     &                + (vn20-amass2(2)-amass2(4))/(-2.D0)
     &                + PP(2,5)
     &                + (vn58-amass2(2)-amass2(6))/(-2.D0)
     &               )
        W_val  = sqrt(  max( W_val2, 0.D0 )  )
        if ((W_val .LT. W_cut(1)).or.(W_val .GT. W_cut(2))) GOTO 999
      endif
      do i=1,4                                              ! Energy cut
        ii = i + 2
        if (   (pe(4,ii) .LT. E_min(i))
     +     .or.(pe(4,ii) .GT. E_max(i)) )  GOTO 999
      enddo
      do i=1,4                                              ! Momentum cut
        ii = i + 2
        Pt_2(i) = pe(1,ii)**2 + pe(2,ii)**2
        momt(i) = sqrt( Pt_2(i) + pe(3,ii)**2 )
        if (   (momt(i) .LT. P_min(i))
     +     .or.(momt(i) .GT. P_max(i)) )  GOTO 999
      enddo
      deg = 180.D0 / pi
      do i=1,4
        theta(i) = acos( pe(3,i+2)/momt(i) )
        the_d(i) = theta(i) * deg
        if (   (the_d(i) .LT. theta_min(i))
     +     .or.(the_d(i) .GT. theta_max(i)) )  GOTO 999
      enddo
      do i=1,4                                              ! Pt cut
        Pt(i) = sqrt( Pt_2(i) )
        if (   (Pt(i) .LT. Pt_min(i))
     +     .or.(Pt(i) .GT. Pt_max(i)) )  GOTO 999
      enddo
      if (lPtMAX_cut) then             !!! PtMAX_cut !!!
        do i=2,4
          lvisi_PtMAX(i) = .false.
          if (     (the_d(i) .GT. the_PtMAX_cut(1))
     +        .and.(the_d(i) .LT. the_PtMAX_cut(2)) ) then
            lvisi_PtMAX(i) = .true.
          endif
        enddo
        if     ( (.not.lvisi_PtMAX(2))
     +      .and.(.not.lvisi_PtMAX(3))
     +      .and.(.not.lvisi_PtMAX(4)) ) then      ! no visi.
          GOTO 999
        elseif ( (.not.lvisi_PtMAX(2))
     +      .and.(.not.lvisi_PtMAX(3))
     +      .and.(     lvisi_PtMAX(4)) ) then      ! only 4 visi.
          if (    (Pt(4) .LT. PtMAX_cut(1))
     +        .OR.(Pt(4) .GT. PtMAX_cut(2)) )  GOTO 999
        elseif ( (.not.lvisi_PtMAX(2))
     +      .and.(     lvisi_PtMAX(3))
     +      .and.(.not.lvisi_PtMAX(4)) ) then      ! only 3 visi.
          if (    (Pt(3) .LT. PtMAX_cut(1))
     +        .OR.(Pt(3) .GT. PtMAX_cut(2)) )  GOTO 999
        elseif ( (     lvisi_PtMAX(2))
     +      .and.(.not.lvisi_PtMAX(3))
     +      .and.(.not.lvisi_PtMAX(4)) ) then      ! only 2 visi.
          if (    (Pt(2) .LT. PtMAX_cut(1))
     +        .OR.(Pt(2) .GT. PtMAX_cut(2)) )  GOTO 999
        elseif ( (.not.lvisi_PtMAX(2))
     +      .and.(     lvisi_PtMAX(3))
     +      .and.(     lvisi_PtMAX(4)) ) then      ! only 3,4 visi.
          PtMAX = max( Pt(3), Pt(4) )
          if (    (PtMAX .LT. PtMAX_cut(1))
     +        .OR.(PtMAX .GT. PtMAX_cut(2)) )  GOTO 999
        elseif ( (     lvisi_PtMAX(2))
     +      .and.(.not.lvisi_PtMAX(3))
     +      .and.(     lvisi_PtMAX(4)) ) then      ! only 2,4 visi.
          PtMAX = max( Pt(2),Pt(4) )
          if (    (PtMAX .LT. PtMAX_cut(1))
     +        .OR.(PtMAX .GT. PtMAX_cut(2)) )  GOTO 999
        elseif ( (     lvisi_PtMAX(2))
     +      .and.(     lvisi_PtMAX(3))
     +      .and.(.not.lvisi_PtMAX(4)) ) then      ! only 2,3 visi.
          PtMAX = max( Pt(2),Pt(3) )
          if (    (PtMAX .LT. PtMAX_cut(1))
     +        .OR.(PtMAX .GT. PtMAX_cut(2)) )  GOTO 999
        elseif ( (     lvisi_PtMAX(2))
     +      .and.(     lvisi_PtMAX(3))
     +      .and.(     lvisi_PtMAX(4)) ) then      ! all of 2,3,4 visi.
          PtMAX = max( max(Pt(2),Pt(3)), Pt(4) )
          if (    (PtMAX .LT. PtMAX_cut(1))
     +        .OR.(PtMAX .GT. PtMAX_cut(2)) )  GOTO 999
        else
          write(6,*) '!!!Error in kinem_4f!!!'
          write(6,*) ' ---> PtMAX_cut exception'
          write(6,*) ' ---> Please contact the author.'
          write(6,*) ' ---> Good-bye!'
          STOP
        endif
      endif
* (M_ell cut)
      mass_ell = sqrt( max(SSS2,0D0) )
      if (  (mass_ell .LT. MassELL_cut(1))
     +  .OR.(mass_ell .GT. MassELL_cut(2))  )  GOTO 999
* (M_qll cut)
      if (process .EQ. 3) then                          
        if (lee_int) then                               
           mass_qll(1) = sqrt( max( min(vn22,vn56),0D0 ) )
           mass_qll(2) = sqrt( max( max(vn22,vn56),0D0 ) )
         else
           mass_qll(1) = sqrt(  max( vn22,0D0 )  )
           mass_qll(2) = (MassQLL_cut(3)+MassQLL_cut(4)) /2D0
        endif
        if (    (mass_qll(1) .GT. MassQLL_cut(1))
     +    .AND. (mass_qll(1) .LT. MassQLL_cut(2))
     +    .AND. (mass_qll(2) .GE. MassQLL_cut(3))
     +    .AND. (mass_qll(2) .LE. MassQLL_cut(4))  ) then
           continue
        else
           GOTO 999
        endif
      endif
      if (lpair.EQ.1 .AND. lee_int) then
        mass_ee(1) = sqrt(  max( MIN(vn30,vn48), 0.D0 )  )
        mass_ee(2) = sqrt(  max( MAX(vn30,vn48), 0.D0 )  )
        if (    (mass_ee(1) .GT. Mass56_cut(1))
     +    .AND. (mass_ee(1) .LT. Mass56_cut(2))
     +    .AND. (mass_ee(2) .GT. Mass56_cut(3))
     +    .AND. (mass_ee(2) .LT. Mass56_cut(4))  )  then
           continue
         else
           GOTO 999
        endif
      endif
      if (lHOT_flag) then
*************** Hot selection on ee channel ***************
      if (lveto) then
        lvisi = .false.
        lvisi3 =   (the_d(1) .GT. the_veto_min(1))
     +        .AND.(the_d(1) .LT. the_veto_max(1))
     +        .AND.(pe(4,3)  .GT. ene_veto_min(1))
     +        .AND.(pe(4,3)  .LT. ene_veto_max(1))
     +        .AND.(Pt(1)    .GT. Pt_veto_min(1))
     +        .AND.(Pt(1)    .LT. Pt_veto_max(1))
        lvisi4 =   (the_d(2) .GT. the_veto_min(2))
     +        .AND.(the_d(2) .LT. the_veto_max(2))
     +        .AND.(pe(4,4)  .GT. ene_veto_min(2))
     +        .AND.(pe(4,4)  .LT. ene_veto_max(2))
     +        .AND.(Pt(2)    .GT. Pt_veto_min(2))
     +        .AND.(Pt(2)    .LT. Pt_veto_max(2))
        lvisi5 =   (the_d(3) .GT. the_veto_min(3))
     +        .AND.(the_d(3) .LT. the_veto_max(3))
     +        .AND.(pe(4,5)  .GT. ene_veto_min(3))
     +        .AND.(pe(4,5)  .LT. ene_veto_max(3))
     +        .AND.(Pt(3)    .GT. Pt_veto_min(3))
     +        .AND.(Pt(3)    .LT. Pt_veto_max(3))
        lvisi6 =   (the_d(4) .GT. the_veto_min(4))
     +        .AND.(the_d(4) .LT. the_veto_max(4))
     +        .AND.(pe(4,6)  .GT. ene_veto_min(4))
     +        .AND.(pe(4,6)  .LT. ene_veto_max(4))
     +        .AND.(Pt(4)    .GT. Pt_veto_min(4))
     +        .AND.(Pt(4)    .LT. Pt_veto_max(4))
        if   (Iveto .EQ. 1) then
          if (lvisi3.or.lvisi4.or.lvisi5.or.lvisi6)  lvisi = .true.
        elseif (Iveto .EQ. 2) then
          if (        (lvisi3.AND.lvisi4)
     +            .OR.(lvisi3.AND.lvisi5)
     +            .OR.(lvisi3.AND.lvisi6)
     +            .OR.(lvisi4.AND.lvisi5)
     +            .OR.(lvisi4.AND.lvisi6)
     +            .OR.(lvisi5.AND.lvisi6)
     +        )  lvisi = .true.
        elseif (Iveto .EQ. 3) then
          if (lvisi4.AND.lvisi5.AND.lvisi6)  lvisi = .true.
        elseif (Iveto .EQ. 4) then
          if (lvisi3.AND.lvisi4.AND.lvisi5.AND.lvisi6)  lvisi = .true.
        elseif (Iveto .EQ.56) then
          if (lvisi5.AND.lvisi6)  lvisi = .true.
        else
          write(6,*) '!!!Error in kinem_4f!!!'
          write(6,*) '  ---> unknown IVETO(=',Iveto,' )'
          write(6,*) '  ---> Good-bye!'
          STOP
        endif
        if (lvisi)  GOTO 999
      endif
      if (Ivisi .LT. 0)  GOTO 777
      lvisi = .false.
      if (Ivisi .GT. 0) then
        lvisi3 =   (the_d(1) .GT. the_visi_min(1))
     +        .AND.(the_d(1) .LT. the_visi_max(1))
     +        .AND.(pe(4,3)  .GT. ene_visi_min(1))
     +        .AND.(pe(4,3)  .LT. ene_visi_max(1))
     +        .AND.(Pt(1)    .GT. Pt_visi_min(1))
     +        .AND.(Pt(1)    .LT. Pt_visi_max(1))
        lvisi4 =   (the_d(2) .GT. the_visi_min(2))
     +        .AND.(the_d(2) .LT. the_visi_max(2))
     +        .AND.(pe(4,4)  .GT. ene_visi_min(2))
     +        .AND.(pe(4,4)  .LT. ene_visi_max(2))
     +        .AND.(Pt(2)    .GT. Pt_visi_min(2))
     +        .AND.(Pt(2)    .LT. Pt_visi_max(2))
        lvisi5 =   (the_d(3) .GT. the_visi_min(3))
     +        .AND.(the_d(3) .LT. the_visi_max(3))
     +        .AND.(pe(4,5)  .GT. ene_visi_min(3))
     +        .AND.(pe(4,5)  .LT. ene_visi_max(3))
     +        .AND.(Pt(3)    .GT. Pt_visi_min(3))
     +        .AND.(Pt(3)    .LT. Pt_visi_max(3))
        lvisi6 =   (the_d(4) .GT. the_visi_min(4))
     +        .AND.(the_d(4) .LT. the_visi_max(4))
     +        .AND.(pe(4,6)  .GT. ene_visi_min(4))
     +        .AND.(pe(4,6)  .LT. ene_visi_max(4))
     +        .AND.(Pt(4)    .GT. Pt_visi_min(4))
     +        .AND.(Pt(4)    .LT. Pt_visi_max(4))
        if   (Ivisi .EQ. 1) then
          if (lvisi3.or.lvisi4.or.lvisi5.or.lvisi6)  lvisi = .true.
        elseif (Ivisi .EQ. 2) then
          if (        (lvisi3.AND.lvisi4)
     +            .OR.(lvisi3.AND.lvisi5)
     +            .OR.(lvisi3.AND.lvisi6)
     +            .OR.(lvisi4.AND.lvisi5)
     +            .OR.(lvisi4.AND.lvisi6)
     +            .OR.(lvisi5.AND.lvisi6)
     +        )  lvisi = .true.
        elseif (Ivisi .EQ. 3) then
          if (lvisi4.AND.lvisi5.AND.lvisi6)  lvisi = .true.
        elseif (Ivisi .EQ. 4) then
          if (lvisi3.AND.lvisi4.AND.lvisi5.AND.lvisi6)  lvisi = .true.
        elseif (Ivisi .EQ.56) then
          if (lvisi5.AND.lvisi6)  lvisi = .true.
        else
          write(6,*) '!!!Error in kinem_4f!!!'
          write(6,*) '  ---> unknown IVISI(=',Ivisi,' )'
          write(6,*) '  ---> Good-bye!'
          STOP
        endif
        if (lvisi)  GOTO 777
      endif   ! (Ivisi .GT. 0)
 1000 continue


      aM3_12 = sqrt(max(vn48,0D0))
      aM3_23 = sqrt(max(vn30,0D0))
      vn46 = + 2.0d0*pp(4,6) + amass2(4) + amass2(6)
      aM3_31 = sqrt(max(vn46,0D0))

      Q2_1_ME = abs(vn20)
      vn36 = - 2.0d0*pp(2,5) + amass2(2) + amass2(5)
      Q2_2_ME = abs(vn36)
      Q2_3_ME = abs(vn58)

      I = 4
      Iscatt = I
        Q2_1_OB = AMASS2(2)+AMASS2(I) -2D0*( PE(4,I)*E2_lab
C    &                                      -PE(1,I)*0D0
C    &                                      -PE(2,I)*0D0
     &                                      +PE(3,I)*P2_lab )
        Q2_1_OB = max(-Q2_1_OB,0D0)
        do ii = 1, 4 !!! ME
          PE_g(ii) = PE(ii,2) - PE(ii,Iscatt)
        enddo
        PE_gp(1) = PE_g(1) + 0D0
        PE_gp(2) = PE_g(2) + 0D0
        PE_gp(3) = PE_g(3) + P1_lab
        PE_gp(4) = PE_g(4) + E1_lab
        W_1_ME = sqrt( max( 0D0,
     &                  PE_gp(4)*PE_gp(4)
     &                 -PE_gp(1)*PE_gp(1)
     &                 -PE_gp(2)*PE_gp(2)
     &                 -PE_gp(3)*PE_gp(3) ) )
        do ii = 1, 4 !!! OB
          PE_g(ii) = (PE(ii,2)+vec_isr(ii)) - PE(ii,Iscatt)
        enddo
        PE_gp(1) = PE_g(1) + 0D0
        PE_gp(2) = PE_g(2) + 0D0
        PE_gp(3) = PE_g(3) + P1_lab
        PE_gp(4) = PE_g(4) + E1_lab
        W_1_OB = sqrt( max( 0D0,
     &                  PE_gp(4)*PE_gp(4)
     &                 -PE_gp(1)*PE_gp(1)
     &                 -PE_gp(2)*PE_gp(2)
     &                 -PE_gp(3)*PE_gp(3) ) )
        X_1_ME = Q2_1_ME /(Q2_1_ME +W_1_ME*W_1_ME -amp2)
        X_1_OB = Q2_1_OB /(Q2_1_OB +W_1_OB*W_1_OB -amp2)
        Y_1_ME = (Q2_1_ME +W_1_ME*W_1_ME -amp2) /(S_ME-amp2)
        Y_1_OB = (Q2_1_OB +W_1_OB*W_1_OB -amp2) /(S(1)-amp2)

      I = 5
      Iscatt = I
        Q2_2_OB = AMASS2(2)+AMASS2(I) -2D0*( PE(4,I)*E2_lab
C    &                                      -PE(1,I)*0D0
C    &                                      -PE(2,I)*0D0
     &                                      +PE(3,I)*P2_lab )
        Q2_2_OB = max(-Q2_2_OB,0D0)
        do ii = 1, 4 !!! ME
          PE_g(ii) = PE(ii,2) - PE(ii,Iscatt)
        enddo
        PE_gp(1) = PE_g(1) + 0D0
        PE_gp(2) = PE_g(2) + 0D0
        PE_gp(3) = PE_g(3) + P1_lab
        PE_gp(4) = PE_g(4) + E1_lab
        W_2_ME = sqrt( max( 0D0,
     &                  PE_gp(4)*PE_gp(4)
     &                 -PE_gp(1)*PE_gp(1)
     &                 -PE_gp(2)*PE_gp(2)
     &                 -PE_gp(3)*PE_gp(3) ) )
        do ii = 1, 4 !!! OB
          PE_g(ii) = (PE(ii,2)+vec_isr(ii)) - PE(ii,Iscatt)
        enddo
        PE_gp(1) = PE_g(1) + 0D0
        PE_gp(2) = PE_g(2) + 0D0
        PE_gp(3) = PE_g(3) + P1_lab
        PE_gp(4) = PE_g(4) + E1_lab
        W_2_OB = sqrt( max( 0D0,
     &                  PE_gp(4)*PE_gp(4)
     &                 -PE_gp(1)*PE_gp(1)
     &                 -PE_gp(2)*PE_gp(2)
     &                 -PE_gp(3)*PE_gp(3) ) )
        X_2_ME = Q2_2_ME /(Q2_2_ME +W_2_ME*W_2_ME -amp2)
        X_2_OB = Q2_2_OB /(Q2_2_OB +W_2_OB*W_2_OB -amp2)
        Y_2_ME = (Q2_2_ME +W_2_ME*W_2_ME -amp2) /(S_ME-amp2)
        Y_2_OB = (Q2_2_OB +W_2_OB*W_2_OB -amp2) /(S(1)-amp2)

      I = 6
      Iscatt = I
        Q2_3_OB = AMASS2(2)+AMASS2(I) -2D0*( PE(4,I)*E2_lab
C    &                                      -PE(1,I)*0D0
C    &                                      -PE(2,I)*0D0
     &                                      +PE(3,I)*P2_lab )
        Q2_3_OB = max(-Q2_3_OB,0D0)
        do ii = 1, 4 !!! ME
          PE_g(ii) = PE(ii,2) - PE(ii,Iscatt)
        enddo
        PE_gp(1) = PE_g(1) + 0D0
        PE_gp(2) = PE_g(2) + 0D0
        PE_gp(3) = PE_g(3) + P1_lab
        PE_gp(4) = PE_g(4) + E1_lab
        W_3_ME = sqrt( max( 0D0,
     &                  PE_gp(4)*PE_gp(4)
     &                 -PE_gp(1)*PE_gp(1)
     &                 -PE_gp(2)*PE_gp(2)
     &                 -PE_gp(3)*PE_gp(3) ) )
        do ii = 1, 4 !!! OB
          PE_g(ii) = (PE(ii,2)+vec_isr(ii)) - PE(ii,Iscatt)
        enddo
        PE_gp(1) = PE_g(1) + 0D0
        PE_gp(2) = PE_g(2) + 0D0
        PE_gp(3) = PE_g(3) + P1_lab
        PE_gp(4) = PE_g(4) + E1_lab
        W_3_OB = sqrt( max( 0D0,
     &                  PE_gp(4)*PE_gp(4)
     &                 -PE_gp(1)*PE_gp(1)
     &                 -PE_gp(2)*PE_gp(2)
     &                 -PE_gp(3)*PE_gp(3) ) )
        X_3_ME = Q2_3_ME /(Q2_3_ME +W_3_ME*W_3_ME -amp2)
        X_3_OB = Q2_3_OB /(Q2_3_OB +W_3_OB*W_3_OB -amp2)
        Y_3_ME = (Q2_3_ME +W_3_ME*W_3_ME -amp2) /(S_ME-amp2)
        Y_3_OB = (Q2_3_OB +W_3_OB*W_3_OB -amp2) /(S(1)-amp2)


      if (l3D_cut) then
        do iL = 1, 3
          i = iL + 1
          ii = i + 2
          do iCut = 1, 3
            l3D_flag(iCut, iL) = .false.
            if ( .true.
     +        .AND.(the_d(i) .GE. the_3D_min(iCut))
     +        .AND.(the_d(i) .LE. the_3D_max(iCut))
     +        .AND.(pe(4,ii) .GT. E_3D_min(iCut))
     +        .AND.(pe(4,ii) .LT. E_3D_max(iCut))
     +        .AND.(momt(i)  .GT. P_3D_min(iCut))
     +        .AND.(momt(i)  .LT. P_3D_max(iCut))
     +        .AND.(Pt(i)    .GE. Pt_3D_min(iCut))
     +        .AND.(Pt(i)    .LT. Pt_3D_max(iCut))
     +         ) then
              l3D_flag(iCut, iL) = .true.
            endif
          enddo
        enddo
        if ( .false. !    <e+->               <l-+>               <l+->
     +    .OR. ( l3D_flag(1,1) .and. l3D_flag(2,2) .and. l3D_flag(3,3)
     +         .and. aM3_12.GT.Mass12_3D(1) .and. aM3_12.LT.Mass12_3D(2)

     +         .and. Q2_3_ME.GE.Q2_3D_ME(1) .and. Q2_3_ME.LT.Q2_3D_ME(2)
     +         .and. Q2_3_OB.GE.Q2_3D_OB(1) .and. Q2_3_OB.LT.Q2_3D_OB(2)
     +         .and. W_3_ME.GE.W_3D_ME(1) .and. W_3_ME.LT.W_3D_ME(2)
     +         .and. W_3_OB.GE.W_3D_OB(1) .and. W_3_OB.LT.W_3D_OB(2)
     +         .and. X_3_ME.GE.X_3D_ME(1) .and. X_3_ME.LT.X_3D_ME(2)
     +         .and. X_3_OB.GE.X_3D_OB(1) .and. X_3_OB.LT.X_3D_OB(2)
     +         .and. Y_3_ME.GE.Y_3D_ME(1) .and. Y_3_ME.LT.Y_3D_ME(2)
     +         .and. Y_3_OB.GE.Y_3D_OB(1) .and. Y_3_OB.LT.Y_3D_OB(2)

     +         )
     +    .OR. ( l3D_flag(1,1) .and. l3D_flag(3,2) .and. l3D_flag(2,3)
     +         .and. aM3_31.GT.Mass12_3D(1) .and. aM3_31.LT.Mass12_3D(2)
     +         .and. Q2_2_ME.GE.Q2_3D_ME(1) .and. Q2_2_ME.LT.Q2_3D_ME(2)
     +         .and. Q2_2_OB.GE.Q2_3D_OB(1) .and. Q2_2_OB.LT.Q2_3D_OB(2)
     +         .and. W_2_ME.GE.W_3D_ME(1) .and. W_2_ME.LT.W_3D_ME(2)
     +         .and. W_2_OB.GE.W_3D_OB(1) .and. W_2_OB.LT.W_3D_OB(2)
     +         .and. X_2_ME.GE.X_3D_ME(1) .and. X_2_ME.LT.X_3D_ME(2)
     +         .and. X_2_OB.GE.X_3D_OB(1) .and. X_2_OB.LT.X_3D_OB(2)
     +         .and. Y_2_ME.GE.Y_3D_ME(1) .and. Y_2_ME.LT.Y_3D_ME(2)
     +         .and. Y_2_OB.GE.Y_3D_OB(1) .and. Y_2_OB.LT.Y_3D_OB(2)

     +         )
     +    .OR. ( l3D_flag(2,1) .and. l3D_flag(1,2) .and. l3D_flag(3,3)
     +         .and. aM3_12.GT.Mass12_3D(1) .and. aM3_12.LT.Mass12_3D(2)

     +         .and. Q2_3_ME.GE.Q2_3D_ME(1) .and. Q2_3_ME.LT.Q2_3D_ME(2)
     +         .and. Q2_3_OB.GE.Q2_3D_OB(1) .and. Q2_3_OB.LT.Q2_3D_OB(2)
     +         .and. W_3_ME.GE.W_3D_ME(1) .and. W_3_ME.LT.W_3D_ME(2)
     +         .and. W_3_OB.GE.W_3D_OB(1) .and. W_3_OB.LT.W_3D_OB(2)
     +         .and. X_3_ME.GE.X_3D_ME(1) .and. X_3_ME.LT.X_3D_ME(2)
     +         .and. X_3_OB.GE.X_3D_OB(1) .and. X_3_OB.LT.X_3D_OB(2)
     +         .and. Y_3_ME.GE.Y_3D_ME(1) .and. Y_3_ME.LT.Y_3D_ME(2)
     +         .and. Y_3_OB.GE.Y_3D_OB(1) .and. Y_3_OB.LT.Y_3D_OB(2)

     +         )
     +    .OR. ( l3D_flag(2,1) .and. l3D_flag(3,2) .and. l3D_flag(1,3)
     +         .and. aM3_31.GT.Mass12_3D(1) .and. aM3_31.LT.Mass12_3D(2)

     +         .and. Q2_2_ME.GE.Q2_3D_ME(1) .and. Q2_2_ME.LT.Q2_3D_ME(2)
     +         .and. Q2_2_OB.GE.Q2_3D_OB(1) .and. Q2_2_OB.LT.Q2_3D_OB(2)
     +         .and. W_2_ME.GE.W_3D_ME(1) .and. W_2_ME.LT.W_3D_ME(2)
     +         .and. W_2_OB.GE.W_3D_OB(1) .and. W_2_OB.LT.W_3D_OB(2)
     +         .and. X_2_ME.GE.X_3D_ME(1) .and. X_2_ME.LT.X_3D_ME(2)
     +         .and. X_2_OB.GE.X_3D_OB(1) .and. X_2_OB.LT.X_3D_OB(2)
     +         .and. Y_2_ME.GE.Y_3D_ME(1) .and. Y_2_ME.LT.Y_3D_ME(2)
     +         .and. Y_2_OB.GE.Y_3D_OB(1) .and. Y_2_OB.LT.Y_3D_OB(2)

     +         )
     +    .OR. ( l3D_flag(3,1) .and. l3D_flag(1,2) .and. l3D_flag(2,3)
     +         .and. aM3_23.GT.Mass12_3D(1) .and. aM3_23.LT.Mass12_3D(2)

     +         .and. Q2_1_ME.GE.Q2_3D_ME(1) .and. Q2_1_ME.LT.Q2_3D_ME(2)
     +         .and. Q2_1_OB.GE.Q2_3D_OB(1) .and. Q2_1_OB.LT.Q2_3D_OB(2)
     +         .and. W_1_ME.GE.W_3D_ME(1) .and. W_1_ME.LT.W_3D_ME(2)
     +         .and. W_1_OB.GE.W_3D_OB(1) .and. W_1_OB.LT.W_3D_OB(2)
     +         .and. X_1_ME.GE.X_3D_ME(1) .and. X_1_ME.LT.X_3D_ME(2)
     +         .and. X_1_OB.GE.X_3D_OB(1) .and. X_1_OB.LT.X_3D_OB(2)
     +         .and. Y_1_ME.GE.Y_3D_ME(1) .and. Y_1_ME.LT.Y_3D_ME(2)
     +         .and. Y_1_OB.GE.Y_3D_OB(1) .and. Y_1_OB.LT.Y_3D_OB(2)

     +         )
     +    .OR. ( l3D_flag(3,1) .and. l3D_flag(2,2) .and. l3D_flag(1,3)
     +         .and. aM3_23.GT.Mass12_3D(1) .and. aM3_23.LT.Mass12_3D(2)

     +         .and. Q2_1_ME.GE.Q2_3D_ME(1) .and. Q2_1_ME.LT.Q2_3D_ME(2)
     +         .and. Q2_1_OB.GE.Q2_3D_OB(1) .and. Q2_1_OB.LT.Q2_3D_OB(2)
     +         .and. W_1_ME.GE.W_3D_ME(1) .and. W_1_ME.LT.W_3D_ME(2)
     +         .and. W_1_OB.GE.W_3D_OB(1) .and. W_1_OB.LT.W_3D_OB(2)
     +         .and. X_1_ME.GE.X_3D_ME(1) .and. X_1_ME.LT.X_3D_ME(2)
     +         .and. X_1_OB.GE.X_3D_OB(1) .and. X_1_OB.LT.X_3D_OB(2)
     +         .and. Y_1_ME.GE.Y_3D_ME(1) .and. Y_1_ME.LT.Y_3D_ME(2)
     +         .and. Y_1_OB.GE.Y_3D_OB(1) .and. Y_1_OB.LT.Y_3D_OB(2)

     +         )
     +     )  GOTO 777
      endif   ! (l3D_cut)

      if (l3E_cut) then
        do iL = 1, 3
          i = iL + 1
          ii = i + 2
          do iCut = 1, 3
            l3E_flag(iCut, iL) = .false.
            if ( .true.
     +        .AND.(the_d(i) .GE. the_3E_min(iCut))
     +        .AND.(the_d(i) .LE. the_3E_max(iCut))
     +        .AND.(pe(4,ii) .GT. E_3E_min(iCut))
     +        .AND.(pe(4,ii) .LT. E_3E_max(iCut))
     +        .AND.(momt(i)  .GT. P_3E_min(iCut))
     +        .AND.(momt(i)  .LT. P_3E_max(iCut))
     +        .AND.(Pt(i)    .GE. Pt_3E_min(iCut))
     +        .AND.(Pt(i)    .LT. Pt_3E_max(iCut))
     +         ) then
              l3E_flag(iCut, iL) = .true.
            endif
          enddo
        enddo
        if ( .false. !    <e+->               <l-+>               <l+->
     +    .OR. ( l3E_flag(1,1) .and. l3E_flag(2,2) .and. l3E_flag(3,3)
     +         .and. aM3_12.GT.Mass12_3E(1) .and. aM3_12.LT.Mass12_3E(2)

     +         .and. Q2_3_ME.GE.Q2_3E_ME(1) .and. Q2_3_ME.LT.Q2_3E_ME(2)
     +         .and. Q2_3_OB.GE.Q2_3E_OB(1) .and. Q2_3_OB.LT.Q2_3E_OB(2)
     +         .and. W_3_ME.GE.W_3E_ME(1) .and. W_3_ME.LT.W_3E_ME(2)
     +         .and. W_3_OB.GE.W_3E_OB(1) .and. W_3_OB.LT.W_3E_OB(2)
     +         .and. X_3_ME.GE.X_3E_ME(1) .and. X_3_ME.LT.X_3E_ME(2)
     +         .and. X_3_OB.GE.X_3E_OB(1) .and. X_3_OB.LT.X_3E_OB(2)
     +         .and. Y_3_ME.GE.Y_3E_ME(1) .and. Y_3_ME.LT.Y_3E_ME(2)
     +         .and. Y_3_OB.GE.Y_3E_OB(1) .and. Y_3_OB.LT.Y_3E_OB(2)

     +         )
     +    .OR. ( l3E_flag(1,1) .and. l3E_flag(3,2) .and. l3E_flag(2,3)
     +         .and. aM3_31.GT.Mass12_3E(1) .and. aM3_31.LT.Mass12_3E(2)

     +         .and. Q2_2_ME.GE.Q2_3E_ME(1) .and. Q2_2_ME.LT.Q2_3E_ME(2)
     +         .and. Q2_2_OB.GE.Q2_3E_OB(1) .and. Q2_2_OB.LT.Q2_3E_OB(2)
     +         .and. W_2_ME.GE.W_3E_ME(1) .and. W_2_ME.LT.W_3E_ME(2)
     +         .and. W_2_OB.GE.W_3E_OB(1) .and. W_2_OB.LT.W_3E_OB(2)
     +         .and. X_2_ME.GE.X_3E_ME(1) .and. X_2_ME.LT.X_3E_ME(2)
     +         .and. X_2_OB.GE.X_3E_OB(1) .and. X_2_OB.LT.X_3E_OB(2)
     +         .and. Y_2_ME.GE.Y_3E_ME(1) .and. Y_2_ME.LT.Y_3E_ME(2)
     +         .and. Y_2_OB.GE.Y_3E_OB(1) .and. Y_2_OB.LT.Y_3E_OB(2)

     +         )
     +    .OR. ( l3E_flag(2,1) .and. l3E_flag(1,2) .and. l3E_flag(3,3)
     +         .and. aM3_12.GT.Mass12_3E(1) .and. aM3_12.LT.Mass12_3E(2)

     +         .and. Q2_3_ME.GE.Q2_3E_ME(1) .and. Q2_3_ME.LT.Q2_3E_ME(2)
     +         .and. Q2_3_OB.GE.Q2_3E_OB(1) .and. Q2_3_OB.LT.Q2_3E_OB(2)
     +         .and. W_3_ME.GE.W_3E_ME(1) .and. W_3_ME.LT.W_3E_ME(2)
     +         .and. W_3_OB.GE.W_3E_OB(1) .and. W_3_OB.LT.W_3E_OB(2)
     +         .and. X_3_ME.GE.X_3E_ME(1) .and. X_3_ME.LT.X_3E_ME(2)
     +         .and. X_3_OB.GE.X_3E_OB(1) .and. X_3_OB.LT.X_3E_OB(2)
     +         .and. Y_3_ME.GE.Y_3E_ME(1) .and. Y_3_ME.LT.Y_3E_ME(2)
     +         .and. Y_3_OB.GE.Y_3E_OB(1) .and. Y_3_OB.LT.Y_3E_OB(2)

     +         )
     +    .OR. ( l3E_flag(2,1) .and. l3E_flag(3,2) .and. l3E_flag(1,3)
     +         .and. aM3_31.GT.Mass12_3E(1) .and. aM3_31.LT.Mass12_3E(2)

     +         .and. Q2_2_ME.GE.Q2_3E_ME(1) .and. Q2_2_ME.LT.Q2_3E_ME(2)
     +         .and. Q2_2_OB.GE.Q2_3E_OB(1) .and. Q2_2_OB.LT.Q2_3E_OB(2)
     +         .and. W_2_ME.GE.W_3E_ME(1) .and. W_2_ME.LT.W_3E_ME(2)
     +         .and. W_2_OB.GE.W_3E_OB(1) .and. W_2_OB.LT.W_3E_OB(2)
     +         .and. X_2_ME.GE.X_3E_ME(1) .and. X_2_ME.LT.X_3E_ME(2)
     +         .and. X_2_OB.GE.X_3E_OB(1) .and. X_2_OB.LT.X_3E_OB(2)
     +         .and. Y_2_ME.GE.Y_3E_ME(1) .and. Y_2_ME.LT.Y_3E_ME(2)
     +         .and. Y_2_OB.GE.Y_3E_OB(1) .and. Y_2_OB.LT.Y_3E_OB(2)

     +         )
     +    .OR. ( l3E_flag(3,1) .and. l3E_flag(1,2) .and. l3E_flag(2,3)
     +         .and. aM3_23.GT.Mass12_3E(1) .and. aM3_23.LT.Mass12_3E(2)

     +         .and. Q2_1_ME.GE.Q2_3E_ME(1) .and. Q2_1_ME.LT.Q2_3E_ME(2)
     +         .and. Q2_1_OB.GE.Q2_3E_OB(1) .and. Q2_1_OB.LT.Q2_3E_OB(2)
     +         .and. W_1_ME.GE.W_3E_ME(1) .and. W_1_ME.LT.W_3E_ME(2)
     +         .and. W_1_OB.GE.W_3E_OB(1) .and. W_1_OB.LT.W_3E_OB(2)
     +         .and. X_1_ME.GE.X_3E_ME(1) .and. X_1_ME.LT.X_3E_ME(2)
     +         .and. X_1_OB.GE.X_3E_OB(1) .and. X_1_OB.LT.X_3E_OB(2)
     +         .and. Y_1_ME.GE.Y_3E_ME(1) .and. Y_1_ME.LT.Y_3E_ME(2)
     +         .and. Y_1_OB.GE.Y_3E_OB(1) .and. Y_1_OB.LT.Y_3E_OB(2)

     +         )
     +    .OR. ( l3E_flag(3,1) .and. l3E_flag(2,2) .and. l3E_flag(1,3)
     +         .and. aM3_23.GT.Mass12_3E(1) .and. aM3_23.LT.Mass12_3E(2)

     +         .and. Q2_1_ME.GE.Q2_3E_ME(1) .and. Q2_1_ME.LT.Q2_3E_ME(2)
     +         .and. Q2_1_OB.GE.Q2_3E_OB(1) .and. Q2_1_OB.LT.Q2_3E_OB(2)
     +         .and. W_1_ME.GE.W_3E_ME(1) .and. W_1_ME.LT.W_3E_ME(2)
     +         .and. W_1_OB.GE.W_3E_OB(1) .and. W_1_OB.LT.W_3E_OB(2)
     +         .and. X_1_ME.GE.X_3E_ME(1) .and. X_1_ME.LT.X_3E_ME(2)
     +         .and. X_1_OB.GE.X_3E_OB(1) .and. X_1_OB.LT.X_3E_OB(2)
     +         .and. Y_1_ME.GE.Y_3E_ME(1) .and. Y_1_ME.LT.Y_3E_ME(2)
     +         .and. Y_1_OB.GE.Y_3E_OB(1) .and. Y_1_OB.LT.Y_3E_OB(2)

     +         )
     +     )  GOTO 777
      endif   ! (l3E_cut)

      if (l3F_cut) then
        do iL = 1, 3
          i = iL + 1
          ii = i + 2
          do iCut = 1, 3
            l3F_flag(iCut, iL) = .false.
            if ( .true.
     +        .AND.(the_d(i) .GE. the_3F_min(iCut))
     +        .AND.(the_d(i) .LE. the_3F_max(iCut))
     +        .AND.(pe(4,ii) .GT. E_3F_min(iCut))
     +        .AND.(pe(4,ii) .LT. E_3F_max(iCut))
     +        .AND.(momt(i)  .GT. P_3F_min(iCut))
     +        .AND.(momt(i)  .LT. P_3F_max(iCut))
     +        .AND.(Pt(i)    .GE. Pt_3F_min(iCut))
     +        .AND.(Pt(i)    .LT. Pt_3F_max(iCut))
     +         ) then
              l3F_flag(iCut, iL) = .true.
            endif
          enddo
        enddo
        if ( .false. !    <e+->               <l-+>               <l+->
     +    .OR. ( l3F_flag(1,1) .and. l3F_flag(2,2) .and. l3F_flag(3,3)
     +         .and. aM3_12.GT.Mass12_3F(1) .and. aM3_12.LT.Mass12_3F(2)

     +         .and. Q2_3_ME.GE.Q2_3F_ME(1) .and. Q2_3_ME.LT.Q2_3F_ME(2)
     +         .and. Q2_3_OB.GE.Q2_3F_OB(1) .and. Q2_3_OB.LT.Q2_3F_OB(2)
     +         .and. W_3_ME.GE.W_3F_ME(1) .and. W_3_ME.LT.W_3F_ME(2)
     +         .and. W_3_OB.GE.W_3F_OB(1) .and. W_3_OB.LT.W_3F_OB(2)
     +         .and. X_3_ME.GE.X_3F_ME(1) .and. X_3_ME.LT.X_3F_ME(2)
     +         .and. X_3_OB.GE.X_3F_OB(1) .and. X_3_OB.LT.X_3F_OB(2)
     +         .and. Y_3_ME.GE.Y_3F_ME(1) .and. Y_3_ME.LT.Y_3F_ME(2)
     +         .and. Y_3_OB.GE.Y_3F_OB(1) .and. Y_3_OB.LT.Y_3F_OB(2)

     +         )
     +    .OR. ( l3F_flag(1,1) .and. l3F_flag(3,2) .and. l3F_flag(2,3)
     +         .and. aM3_31.GT.Mass12_3F(1) .and. aM3_31.LT.Mass12_3F(2)

     +         .and. Q2_2_ME.GE.Q2_3F_ME(1) .and. Q2_2_ME.LT.Q2_3F_ME(2)
     +         .and. Q2_2_OB.GE.Q2_3F_OB(1) .and. Q2_2_OB.LT.Q2_3F_OB(2)
     +         .and. W_2_ME.GE.W_3F_ME(1) .and. W_2_ME.LT.W_3F_ME(2)
     +         .and. W_2_OB.GE.W_3F_OB(1) .and. W_2_OB.LT.W_3F_OB(2)
     +         .and. X_2_ME.GE.X_3F_ME(1) .and. X_2_ME.LT.X_3F_ME(2)
     +         .and. X_2_OB.GE.X_3F_OB(1) .and. X_2_OB.LT.X_3F_OB(2)
     +         .and. Y_2_ME.GE.Y_3F_ME(1) .and. Y_2_ME.LT.Y_3F_ME(2)
     +         .and. Y_2_OB.GE.Y_3F_OB(1) .and. Y_2_OB.LT.Y_3F_OB(2)

     +         )
     +    .OR. ( l3F_flag(2,1) .and. l3F_flag(1,2) .and. l3F_flag(3,3)
     +         .and. aM3_12.GT.Mass12_3F(1) .and. aM3_12.LT.Mass12_3F(2)

     +         .and. Q2_3_ME.GE.Q2_3F_ME(1) .and. Q2_3_ME.LT.Q2_3F_ME(2)
     +         .and. Q2_3_OB.GE.Q2_3F_OB(1) .and. Q2_3_OB.LT.Q2_3F_OB(2)
     +         .and. W_3_ME.GE.W_3F_ME(1) .and. W_3_ME.LT.W_3F_ME(2)
     +         .and. W_3_OB.GE.W_3F_OB(1) .and. W_3_OB.LT.W_3F_OB(2)
     +         .and. X_3_ME.GE.X_3F_ME(1) .and. X_3_ME.LT.X_3F_ME(2)
     +         .and. X_3_OB.GE.X_3F_OB(1) .and. X_3_OB.LT.X_3F_OB(2)
     +         .and. Y_3_ME.GE.Y_3F_ME(1) .and. Y_3_ME.LT.Y_3F_ME(2)
     +         .and. Y_3_OB.GE.Y_3F_OB(1) .and. Y_3_OB.LT.Y_3F_OB(2)

     +         )
     +    .OR. ( l3F_flag(2,1) .and. l3F_flag(3,2) .and. l3F_flag(1,3)
     +         .and. aM3_31.GT.Mass12_3F(1) .and. aM3_31.LT.Mass12_3F(2)

     +         .and. Q2_2_ME.GE.Q2_3F_ME(1) .and. Q2_2_ME.LT.Q2_3F_ME(2)
     +         .and. Q2_2_OB.GE.Q2_3F_OB(1) .and. Q2_2_OB.LT.Q2_3F_OB(2)
     +         .and. W_2_ME.GE.W_3F_ME(1) .and. W_2_ME.LT.W_3F_ME(2)
     +         .and. W_2_OB.GE.W_3F_OB(1) .and. W_2_OB.LT.W_3F_OB(2)
     +         .and. X_2_ME.GE.X_3F_ME(1) .and. X_2_ME.LT.X_3F_ME(2)
     +         .and. X_2_OB.GE.X_3F_OB(1) .and. X_2_OB.LT.X_3F_OB(2)
     +         .and. Y_2_ME.GE.Y_3F_ME(1) .and. Y_2_ME.LT.Y_3F_ME(2)
     +         .and. Y_2_OB.GE.Y_3F_OB(1) .and. Y_2_OB.LT.Y_3F_OB(2)

     +         )
     +    .OR. ( l3F_flag(3,1) .and. l3F_flag(1,2) .and. l3F_flag(2,3)
     +         .and. aM3_23.GT.Mass12_3F(1) .and. aM3_23.LT.Mass12_3F(2)

     +         .and. Q2_1_ME.GE.Q2_3F_ME(1) .and. Q2_1_ME.LT.Q2_3F_ME(2)
     +         .and. Q2_1_OB.GE.Q2_3F_OB(1) .and. Q2_1_OB.LT.Q2_3F_OB(2)
     +         .and. W_1_ME.GE.W_3F_ME(1) .and. W_1_ME.LT.W_3F_ME(2)
     +         .and. W_1_OB.GE.W_3F_OB(1) .and. W_1_OB.LT.W_3F_OB(2)
     +         .and. X_1_ME.GE.X_3F_ME(1) .and. X_1_ME.LT.X_3F_ME(2)
     +         .and. X_1_OB.GE.X_3F_OB(1) .and. X_1_OB.LT.X_3F_OB(2)
     +         .and. Y_1_ME.GE.Y_3F_ME(1) .and. Y_1_ME.LT.Y_3F_ME(2)
     +         .and. Y_1_OB.GE.Y_3F_OB(1) .and. Y_1_OB.LT.Y_3F_OB(2)

     +         )
     +    .OR. ( l3F_flag(3,1) .and. l3F_flag(2,2) .and. l3F_flag(1,3)
     +         .and. aM3_23.GT.Mass12_3F(1) .and. aM3_23.LT.Mass12_3F(2)

     +         .and. Q2_1_ME.GE.Q2_3F_ME(1) .and. Q2_1_ME.LT.Q2_3F_ME(2)
     +         .and. Q2_1_OB.GE.Q2_3F_OB(1) .and. Q2_1_OB.LT.Q2_3F_OB(2)
     +         .and. W_1_ME.GE.W_3F_ME(1) .and. W_1_ME.LT.W_3F_ME(2)
     +         .and. W_1_OB.GE.W_3F_OB(1) .and. W_1_OB.LT.W_3F_OB(2)
     +         .and. X_1_ME.GE.X_3F_ME(1) .and. X_1_ME.LT.X_3F_ME(2)
     +         .and. X_1_OB.GE.X_3F_OB(1) .and. X_1_OB.LT.X_3F_OB(2)
     +         .and. Y_1_ME.GE.Y_3F_ME(1) .and. Y_1_ME.LT.Y_3F_ME(2)
     +         .and. Y_1_OB.GE.Y_3F_OB(1) .and. Y_1_OB.LT.Y_3F_OB(2)

     +         )
     +     )  GOTO 777
      endif   ! (l3F_cut)


      l2e_visiA = .false.
      if (l2e_visiA_flag) then   
      l2e_visi4 = (the_d(2).GT.the_2e(1)).and.(the_d(2).LT.the_2e(3))
     +.and.(  ( (pe(4,4).GT.E_min_2e(1)).and.(the_d(2).LT.the_2e(2)) )
     +    .or.( (pe(4,4).GT.E_min_2e(2)).and.(the_d(2).GT.the_2e(2)) )
     +     )
      l2e_visi5 = (the_d(3).GT.the_2e(1)).and.(the_d(3).LT.the_2e(3))
     +.and.(  ( (pe(4,5).GT.E_min_2e(1)).and.(the_d(3).LT.the_2e(2)) )
     +    .or.( (pe(4,5).GT.E_min_2e(2)).and.(the_d(3).GT.the_2e(2)) )
     +     )
      l2e_visi6 = (the_d(4).GT.the_2e(1)).and.(the_d(4).LT.the_2e(3))
     +.and.(  ( (pe(4,6).GT.E_min_2e(1)).and.(the_d(4).LT.the_2e(2)) )
     +    .or.( (pe(4,6).GT.E_min_2e(2)).and.(the_d(4).GT.the_2e(2)) )
     +     )
      if (  (      (l2e_visi4.AND.l2e_visi5)
     +         .OR.(l2e_visi4.AND.l2e_visi6)
     +         .OR.(l2e_visi5.AND.l2e_visi6)
     +       )  ) then
        l2e_visiA = .true.
      endif
      if (l2e_visiA)  GOTO 777
      endif   
      l2e_visiB = .false.
      if (l2e_visiB_flag) then   
      Ivisi1 = 0
      Ivisi2 = 0
      do 111 i=2,4
        if ( (the_d(i).GT.the_2eB(2)).and.(the_d(i).LT.the_2eB(4))
     + .and.( ((pe(4,i+2).GT.E_min_2eB(1)).and.(the_d(i).LT.the_2eB(3)))
     +    .or.((pe(4,i+2).GT.E_min_2eB(2)).and.(the_d(i).GT.the_2eB(3)))
     +      )
     +     ) then
            Ivisi1 = Ivisi1 + 1
         elseif ( (the_d(i).GT.the_2eB(1)).and.(the_d(i).LT.the_2eB(5))
     + .and.( ((pe(4,i+2).GT.E_min_2eB(3)).and.(the_d(i).LT.the_2eB(3)))
     +    .or.((pe(4,i+2).GT.E_min_2eB(4)).and.(the_d(i).GT.the_2eB(3)))
     +      )
     +     ) then
            Ivisi2 = Ivisi2 + 1
        endif
 111  continue
      if (    (Ivisi1.GE.2)
     +    .or.( (Ivisi1.EQ.1).and.(Ivisi2.GE.1) )
     +   ) then
         l2e_visiB = .true.
      endif
      if (l2e_visiB)  GOTO 777
      endif   
      l3e_visi = .false.
      if (l3e_visi_flag) then   
      l3e_visi4= (the_d(2).GT.the_3e(2)).and.(the_d(2).LT.the_3e(3))
     +      .and.(the_d(3).GT.the_3e(1)).and.(the_d(3).LT.the_3e(4))
     +      .and.(the_d(4).GT.the_3e(1)).and.(the_d(4).LT.the_3e(4))
      l3e_visi5= (the_d(2).GT.the_3e(1)).and.(the_d(2).LT.the_3e(4))
     +      .and.(the_d(3).GT.the_3e(2)).and.(the_d(3).LT.the_3e(3))
     +      .and.(the_d(4).GT.the_3e(1)).and.(the_d(4).LT.the_3e(4))
      l3e_visi6= (the_d(2).GT.the_3e(1)).and.(the_d(2).LT.the_3e(4))
     +      .and.(the_d(3).GT.the_3e(1)).and.(the_d(3).LT.the_3e(4))
     +      .and.(the_d(4).GT.the_3e(2)).and.(the_d(4).LT.the_3e(3))
      if (  ( l3e_visi4 .or. l3e_visi5 .or. l3e_visi6 )
     +     .and.( min( min(pe(4,4),pe(4,5)), pe(4,6) ) .GT. E_min_3e )
     +   ) then
        l3e_visi = .true.
      endif
      if (l3e_visi)  GOTO 777
      endif   
      lPt_visi = .false.
      if (lPt_visi_flag) then   
      do i=2,4   
      do j=1,3   
        if (
     +     (the_d(i).GT.the_min_pt(j)).and.(the_d(i).LT.the_max_pt(j))
     +                                .and.(Pt(i).GT.Pt_min_pt(j))
     +     ) then
          lPt_visi = .true.
        endif
        if (lPt_visi)  GOTO 777
      enddo
      enddo
      endif   
      lemu_visi = .false.
      if (lscattL_flag) then   
        lscattL_visi =  (the_d(2).GT.theta_scattL(1))
     +             .and.(the_d(2).LT.theta_scattL(2))
     +             .and.(pe(4,4).GT.E_scattL(1))
     +             .and.(pe(4,4).LT.E_scattL(2))
        lprodL_visi(1) = (the_d(3).GT.theta_prodL(1))
     +              .and.(the_d(3).LT.theta_prodL(2))
     +              .and.(momt(3).GT.P_prodL(1))
     +              .and.(momt(3).LT.P_prodL(2))
        lprodL_visi(2) = (the_d(4).GT.theta_prodL(1))
     +              .and.(the_d(4).LT.theta_prodL(2))
     +              .and.(momt(4).GT.P_prodL(1))
     +              .and.(momt(4).LT.P_prodL(2))
        lemu_visi = lscattL_visi
     +          .and. (lprodL_visi(1) .or. lprodL_visi(2))
        if (lemu_visi)  GOTO 777
      endif   
      if (   lvisi .or. l2e_visiA .or. l2e_visiB .or. l3e_visi
     +  .or. lPt_visi  .or. lemu_visi
     +   ) then
         GOTO 777
       else
         GOTO 999
      endif
********************************************************************
      endif
* ----------------------------------------------
 777  continue
      if(isr.eq.1) then
       nn=NN_ISR
       if     (ISR_scale .EQ. 1) then
         if (lee_int) then
           QED_scale = min( abs(vn20), abs(vn58) )
         else
           QED_scale = abs(TTT2)
         endif
       elseif (ISR_scale .EQ. 2) then
          QED_scale = S(2)
       else
          write(6,*) '!!!Error in KINEM_4f!!!'
          write(6,*) '  ---> Invalid ISRSCALE =', ISR_scale
          write(6,*) '  ---> Good-bye!'
          STOP
       endif
       QED_scale = max(0.001D0,QED_scale)
       QED_scale = QED_scale *Factor_ISR*Factor_ISR
       al     = log(QED_scale/amass2(2))
       if ((al-1.d0).LE.0) then
         GOTO 999
       endif
       abeta  = 2.d0*alpha0/pi*(al-1.d0)
       radi = hradiO(yy2,al,abeta/2D0) 
        radi = radi *dble(nn) *(y2max/yy2)**(1d0/dble(nn))
cccccccc
      else
       radi=1
      endif
cccccccc
 888  yacob = yacob * AJACOB  *radi
*>>> Boost from LAB to gp...
      if (Frame_amp .EQ. 2) then
         do k=1,NEXTRN
           do i=1,4
             PE_str(i,k) = PE(i,k)
           enddo
         enddo
         PP_gp_LAB(1) = -(PE(1,1) + PE(1,2) - PE(1,4))
         PP_gp_LAB(2) = -(PE(2,1) + PE(2,2) - PE(2,4))
         PP_gp_LAB(3) = -(PE(3,1) + PE(3,2) - PE(3,4))
         PP_gp_LAB(4) =   PE(4,1) + PE(4,2) - PE(4,4)
         if (   abs(PP_gp_LAB(1)/W(1)) .LE. 1D-12
     +    .AND. abs(PP_gp_LAB(2)/W(1)) .LE. 1D-12
     +    .AND. abs(PP_gp_LAB(3)/W(1)) .LE. 1D-12  ) then !!!LAB.eq.CMS_gp!!!
           continue
         else
           do k=1,NEXTRN
             call PBoost( PE(1,k), PP_gp_LAB, PE(1,k) )
           enddo
         endif
      endif
      return
999   continue
      yacob=0
      jump=1
      return
      end
      SUBROUTINE K6NEM0(NEXTRN,X,PE,PP,YACOB,NREG,IREG,JUMP)
      IMPLICIT REAL* 8(A-H,O-Z)
************************************************************************
      INTEGER NEXTRN
      PARAMETER ( MXDIM = 50 )
      COMMON / BPARM1 / XL(MXDIM),XU(MXDIM),NDIM,NWILD,
     &                 IG(MXDIM),NCALL
      REAL*8  X(MXDIM)
      REAL*8  PE(4,NEXTRN), PP(NEXTRN,NEXTRN)
      REAL*8  YACOB
      INTEGER NREG, IREG
      INTEGER JUMP
      INCLUDE 'inclk.h'
      common/ktest/qx56
      common/atest/ttt1,uuu1,ttt2,uuu2,tttt,uuuu,sss2,sss1,q56,q12,q22
     .,d1,g1,t1,u1,t2,u2
      common/kmcntl_4f/iresns(4),icos3,icosq3,icos5,isr,iswap,ident
     &                ,iphi6,ieeee,i34,itag,isym
*--------------------------------------------------------------------
      COMMON/KINEM1/S(3),W(3),FACT
      COMMON/CUT001_4f/COSCUT(2,4),ENGYCT(2,4),AMASCT(2,6),ARESNS(2,4)
     .,opncut,swapm2
      REAL*8 PBOST1(4),PBOST2(4),PK1(4),TT1(4),TT2(4),PDUM(4)
      DATA PDUM/4*1.D0/
      common/TAISR/yy2,y2max
*-----------------------------------------------------------------------
      data nev/0/
*-----------------------------------------------------------------------
      double precision  P1_lab,P2_lab, E1_lab,E2_lab  
     &                 ,Pcms_lab(3),Ecms_lab(3)
     &                 ,GAMMAcms_lab(3),BETGAMcms_lab(3)
     &                 ,vec_isr(4)
       common /GEP_LAB/ P1_lab,P2_lab, E1_lab,E2_lab
     &                 ,Pcms_lab,Ecms_lab
     &                 ,GAMMAcms_lab,   BETGAMcms_lab
     &                 ,vec_isr
*-----------------------------------------------------------------------
      include './inc/graepia.h'
      double precision  DRN, Flux_factor
       external         DRN, Flux_factor
*-----------------------------------------------------------------------
      JUMP   =   0
      YACOB  =   0.D0
      AJACOB=    1
      ipoint = 1
      csut11=coscut(1,1)
      csut21=coscut(2,1)
      csut12=coscut(1,2)
      csut22=coscut(2,2)
      csut13=coscut(1,3)
      csut23=coscut(2,3)
      csut14=coscut(1,4)
      csut24=coscut(2,4)
      if(isr.eq.1) then
         nn=NN_ISR
       xx2min = ((amass1(3)+amass1(4)+amass1(5)+amass1(6))**2-amass2(1))
     &          /(s(2)-amass2(1))
       y2max  = 1.d0-xx2min  
       yy2    = y2max * X(NDIM-ipol_Ebeam)**nn
       if (yy2.lt.1D-300) yy2=1D-300
       S(3)    = (S(2)-amass2(1))*(1.d0-yy2) +amass2(1)
       W(3)    = sqrt(  max( S(3), 0.D0 )  )
       if ( W(3) .LE. (amass1(1)+amass1(2)) )  GOTO 999
       vec_isr(1) = 0.D0            
       vec_isr(2) = 0.D0
       vec_isr(3) = -(S(2)-amass2(1))/2.D0/W(2) *yy2
       vec_isr(4) = abs(vec_isr(3))
       Pz_cms = -vec_isr(3)         
       E_cms  = sqrt(  vec_isr(4)*vec_isr(4)  +  S(3)  )
         Pcms_lab(3) = GAMMAcms_lab(2)  * Pz_cms
     &               + BETGAMcms_lab(2) * E_cms
         Ecms_lab(3) = GAMMAcms_lab(2)  * E_cms
     &               + BETGAMcms_lab(2) * Pz_cms
        GAMMAcms_lab(3)  = Ecms_lab(3)/W(3)
        BETGAMcms_lab(3) = Pcms_lab(3)/W(3)
        Pz_cms = vec_isr(3)         
        E_cms  = vec_isr(4)
          vec_isr(3) = GAMMAcms_lab(2)  * Pz_cms
     &               + BETGAMcms_lab(2) * E_cms
          vec_isr(4) = abs(vec_isr(3))
cccccccc
      else
       S(3) = S(2)
       W(3) = W(2)
       radi=1
       do i=1,4
         vec_isr(i) = 0.D0
       enddo
       Pcms_lab(3) = Pcms_lab(2)
       Ecms_lab(3) = Ecms_lab(2)
       GAMMAcms_lab(3)  = GAMMAcms_lab(2)
       BETGAMcms_lab(3) = BETGAMcms_lab(2)
      endif
cccccccc
      ipoint = 2
      E1      = (S(3)+AMASS2(1)-AMASS2(2))/2.D0/W(3)
      P1      = SQRT( max( E1-AMASS1(1), 0.D0 ) * (E1+AMASS1(1)) )
      E2      = (S(3)+AMASS2(2)-AMASS2(1))/2.D0/W(3)
      P2      = P1
*Particle 1
      PE(1,1) = 0
      PE(2,1) = 0
      PE(3,1) = P1
      PE(4,1) = E1
*Particle 2
      PE(1,2) = 0
      PE(2,2) = 0
      PE(3,2) =-P2
      PE(4,2) = E2
      IF(ISYM.EQ.1) THEN
       IF(X(4).LT.0.5D0) THEN
        X4=X(4)*2.D0
        I34=3
       ELsE
        X4=2.D0-X(4)*2.D0
        I34=4
       END IF
                      AJACOB=AJACOB*2.D0
      ELSE
       X4=X(4)
      END IF
***************
*    QX56     *
*    X=3 or 4 *
***************
       QX56MN=(     AMASS1(7-I34)+AMASS1(5)+AMASS1(6))**2
       QX56MX=(W(3)-AMASS1(  I34)                    )**2
       IF (QX56MX.LE.QX56MN) GOTO 999
       IF (QX56MN.LE.0)      GOTO 999
       IF     (IRESNS(3).EQ. 0) THEN
        QX56=QX56MN+(QX56MX-QX56MN)*X(1)
                      AJACOB=AJACOB*(QX56MX-QX56MN)
       ELSE IF(IRESNS(3).EQ. 1) THEN
        ZM    =ARESNS(1,3)
        ZM2   =ZM*ZM
        ZMG   =ARESNS(1,3)*ARESNS(2,3)
        THEMIN=ATAN((QX56MN-ZM2)/ZMG)
        THEMAX=ATAN((QX56MX-ZM2)/ZMG)
        THE   =THEMIN+(THEMAX-THEMIN)*X(1)
        QX56  =ZMG*TAN(THE)+ZM2
                      AJACOB=AJACOB*(THEMAX-THEMIN)
     .*                 ((QX56-ZM2)**2+ZMG**2)/ZMG
       ELSE IF(IRESNS(3).EQ.-1) THEN
        IF(QX56MN.LE.0) GOTO 999
        QX56=QX56MN*(QX56MX/QX56MN)**X(1)
                      AJACOB=AJACOB*QX56*LOG(QX56MX/QX56MN)
       ELSE IF(IRESNS(3).EQ.-2) THEN
        IF(QX56MN.LE.0) GOTO 999
        t_min = sqrt(QX56MN)
        t_max = sqrt(QX56MX)
        t_diff = t_max - t_min
        t_val = t_min + t_diff * X(1)
        QX56 = t_val**2
                      AJACOB=AJACOB* 2.D0*t_val*t_diff
       ELSE IF(IRESNS(3).EQ. 40) THEN
         DD = QX56MX - QX56MN
         QX56 = QX56MN + DD*X(1)**Rnn_456
             AJACOB = AJACOB *DD *Rnn_456 *X(1)**(Rnn_456-1D0)
       ELSE IF(IRESNS(3).EQ. 2) THEN
	IF(X(1).LT.0.5D0) THEN
	X1=X(1)*2.D0
                      AJACOB=AJACOB*2.D0
        IF(QX56MN.LE.0) GOTO 999
        QX56MX=(QX56MX+QX56MN)/2.D0
        QX56=QX56MN*(QX56MX/QX56MN)**X1
                      AJACOB=AJACOB*QX56*LOG(QX56MX/QX56MN)
        ELSE
	X1=(1.D0-X(1))*2.D0
                      AJACOB=AJACOB*2.D0
        QX56MN=(QX56MX+QX56MN)/2.D0
        QX56=QX56MN+(QX56MX-QX56MN)*X1
                      AJACOB=AJACOB*(QX56MX-QX56MN)
        END IF
       ELSE
        WRITE(6,*)' IRESNS(3) =',IRESNS(3),' is NOT supported. '
        stop
       END IF
       AX56=QX56
       EX   = (S(3)+AMASS2(I34)-QX56)/2.D0/W(3)
       PX   = SQRT( max( EX-AMASS1(I34), 0.D0 ) * (EX+AMASS1(I34)) )
       EX56 = (S(3)+QX56-AMASS2(I34))/2.D0/W(3)
       DE=(AMASS2(I34-2)-AMASS2(5-I34)+QX56-AMASS2(I34))/2.D0/W(3)
       DP=(DE*(EX+PE(4,I34-2))-AMASS2(I34-2)
     .+AMASS2(I34))/(PX+ABS(PE(3,I34-2)))
***************
*    COSX     *
*    X=3 or 4 *
***************
      ipoint = 3
       IF     (ICOS3.EQ.0) THEN
        COSX = COSCUT(1,I34-2)+(COSCUT(2,I34-2)-COSCUT(1,I34-2))*X(3)
                      AJACOB=AJACOB*(COSCUT(2,I34-2)-COSCUT(1,I34-2))
       ELSE IF(ICOS3.EQ.2) THEN
        IF(I34.EQ.3) THEN
         ETA=1.D0-COSCUT(2,I34-2)
	 DMIN=2.D0*
     . (AMASS2(I34  )*(DE/(EX+PX)+AMASS2(I34  )/2.D0/(EX+PX)**2)
     .- AMASS2(I34-2)*(DP/(E1+P1)+AMASS2(I34-2)/2.D0/(E1+P1)**2) )
     .+   2.D0*ETA*P1*PX
         ETA=1.D0-COSCUT(1,I34-2)
	 DMAX=2.D0*
     . (AMASS2(I34  )*(DE/(EX+PX)+AMASS2(I34  )/2.D0/(EX+PX)**2)
     .- AMASS2(I34-2)*(DP/(E1+P1)+AMASS2(I34-2)/2.D0/(E1+P1)**2) )
     .+   2.D0*ETA*P1*PX
         IF(DMIN.GT.DMAX) GOTO 999
         nev=nev+1
         if(dmin.le.0) goto 999
         D=DMIN*(DMAX/DMIN)**X(3)
                      AJACOB=AJACOB*D*LOG(DMAX/DMIN)/P1/PX/2.D0
         PP13=(D+AMASS2(1)+AMASS2(I34))/2.D0
         XXXX=AMASS2(I34  )*(DE/(EX+PX)+AMASS2(I34  )/2.D0/(EX+PX)**2)
     .-       AMASS2(I34-2)*(DP/(E1+P1)+AMASS2(I34-2)/2.D0/(E1+P1)**2)
     .-D/2.D0
         COSX=1.D0+XXXX/P1/PX
	 SINX=SQRT(-XXXX/P1/PX*(2.D0+XXXX/P1/PX))
	 TTT1=-D
        ELSE
         ETA=1.D0+COSCUT(1,I34-2)
	 DMIN=2.D0*
     . (AMASS2(I34  )*(DE/(EX+PX)+AMASS2(I34  )/2.D0/(EX+PX)**2)
     .- AMASS2(I34-2)*(DP/(E2+P2)+AMASS2(I34-2)/2.D0/(E2+P2)**2) )
     .+   2.D0*ETA*P2*PX
         ETA=1.D0+COSCUT(2,I34-2)
	 DMAX=2.D0*
     . (AMASS2(I34  )*(DE/(EX+PX)+AMASS2(I34  )/2.D0/(EX+PX)**2)
     .- AMASS2(I34-2)*(DP/(E2+P2)+AMASS2(I34-2)/2.D0/(E2+P2)**2) )
     .+   2.D0*ETA*P2*PX
         IF(DMIN.GT.DMAX) GOTO 999
         D=DMIN*(DMAX/DMIN)**X(3)
                      AJACOB=AJACOB*D*LOG(DMAX/DMIN)/P2/PX/2.D0
         PP24=(D+AMASS2(2)+AMASS2(I34))/2.D0
         XXXX=AMASS2(I34  )*(DE/(EX+PX)+AMASS2(I34  )/2.D0/(EX+PX)**2)
     .-       AMASS2(I34-2)*(DP/(E2+P2)+AMASS2(I34-2)/2.D0/(E2+P2)**2)
     .-D/2.D0
         COSX=1.D0+XXXX/P2/PX
	 SINX=SQRT(-XXXX/P1/PX*(2.D0+XXXX/P1/PX))
         COSX=-COSX
	 TTT2=-D
        END IF
       ELSE
        WRITE(6,*)' ICOS3 =',ICOS3,' is NOT supported. '
        stop
       END IF
       IF(ABS(COSX).GT.1) GOTO 999
***************
*    PHIX     *
*    X=3 or 4 *
***************
      if (Lpol_Ebeam .AND. ipol_Ebeam.GT.0) then
         PHIX = X(NDIM)*2D0*PI
      else
         if (ROT_kinem_flag) then
            PHIX = DRN(idummy) *2D0*PI
          else
            PHIX = 0D0
         endif
      endif
                      AJACOB=AJACOB*2.D0*PI
      PHI1 = PHIX
*Particle X
      PE(1,I34) = PX*SINX*COS(PHIX)
      PE(2,I34) = PX*SINX*SIN(PHIX)
      PE(3,I34) = PX*COSX
      PE(4,I34) = EX
*Particle 456
      PBOST1(1) =-PX*SINX*COS(PHIX)
      PBOST1(2) =-PX*SINX*SIN(PHIX)
      PBOST1(3) =-PX*COSX
      PBOST1(4) = EX56
***************
*    Q56      *
***************
      ipoint = 4
       Q56MN=MAX((           AMASS1(5)+AMASS1(6))**2,AMASCT(1,2)**2)
       Q56MX=MIN((SQRT(aX56)-AMASS1(7-I34)      )**2,AMASCT(2,2)**2)
       IF(Q56MX.LE.Q56MN) GOTO 999
       IF     (IRESNS(2).EQ. 0) THEN
        Q56=Q56MN+(Q56MX-Q56MN)*X(2)
                      AJACOB=AJACOB*(Q56MX-Q56MN)
       ELSE IF(IRESNS(2).EQ. 40) THEN
         DD = Q56MX - Q56MN
         if (num_reso_mj .EQ. 0) then   
           if (Diff_nn .EQ. 0) then
             Q56 = Q56MN + DD*X(2)**Inn
             AJACOB = AJACOB *DD*Dnn*X(2)**(Inn-1)
           else
             Q56 = Q56MN + DD*X(2)**Dnn
             AJACOB = AJACOB *DD*Dnn*X(2)**(Dnn-1D0)
           endif
         else
           RN_mj = DRN(idummy)   
           do j = 0, num_reso_mj   
             if (RN_mj .LT. weisel_mj(j)) GOTO 4100
           enddo
           write(6,*) '!!!Warning in KINEM_4f!!!'
           write(6,*) '  ---> Inconsistent Weight_mj members!'
           write(6,*) '  ---> Weights_mj ='
     &                ,(weight_mj(i),i=0,num_reso_mj)
           write(6,*) '  ---> Good-bye!'
           STOP
 4100      continue
           if (j .EQ. 0) then
             if (Diff_nn .EQ. 0) then
               Q56 = Q56MN + DD*X(2)**Inn
             else
               Q56 = Q56MN + DD*X(2)**Dnn
             endif
           else
             AMG = mass_mj(j)*width_mj(j)
             AM2 = mass2_mj(j)
             FFmin = atan((Q56MN-AM2)/AMG)  
             FFmax = atan((Q56MX-AM2)/AMG)  
C             dFF = FFmax - FFmin
             Q56 = AM2 + AMG*tan(FFmin+(FFmax-FFmin)*X(2))
           endif
           Density = 0D0
           if (weight_mj(0) .GT. 0) then
             XX_inv = ((Q56-Q56MN)/DD)**(1D0/Dnn)
             if (Diff_nn .EQ. 0) then
               Y_inv = DD*Dnn*XX_inv**(Inn-1D0)
             else
               Y_inv = DD*Dnn*XX_inv**(Dnn-1D0)
             endif
             Density = Density + weight_mj(0)/Y_inv
           endif
           do j = 1, num_reso_mj
             AMG = mass_mj(j)*width_mj(j)
             AM2 = mass2_mj(j)
             FFmin = atan((Q56MN-AM2)/AMG)  
             FFmax = atan((Q56MX-AM2)/AMG)  
C             dFF = FFmax - FFmin
             Y_inv = ((Q56-AM2)**2 + AMG*AMG) *(FFmax-FFmin)/AMG
             Density = Density + weight_mj(j)/Y_inv
           enddo
           AJACOB = AJACOB *(1D0/Density)
         endif
       ELSE IF(IRESNS(2).EQ. 1) THEN
        ZM    =ARESNS(1,2)
        ZM2   =ZM*ZM
        ZMG   =ARESNS(1,2)*ARESNS(2,2)
        THEMIN=ATAN((Q56MN-ZM2)/ZMG)
        THEMAX=ATAN((Q56MX-ZM2)/ZMG)
        THE   =THEMIN+(THEMAX-THEMIN)*X(2)
        Q56   =ZMG*TAN(THE)+ZM2
                      AJACOB=AJACOB*(THEMAX-THEMIN)
     .*                 ((Q56-ZM2)**2+ZMG**2)/ZMG
       ELSE IF(IRESNS(2).EQ.-1) THEN
        IF(Q56MN.LE.0) goto 999
        Q56=Q56MN*(Q56MX/Q56MN)**X(2)
                      AJACOB=AJACOB*Q56*LOG(Q56MX/Q56MN)
       ELSE IF(IRESNS(2).EQ.-2) THEN
        IF(Q56MN.LE.0) GOTO 999
        t_min = sqrt(Q56MN)
        t_max = sqrt(Q56MX)
        t_diff = t_max - t_min
        t_val = t_min + t_diff * X(2)
        Q56 = t_val**2
                      AJACOB=AJACOB* 2.D0*t_val*t_diff
       ELSE IF(IRESNS(2).EQ. 2) THEN
        IF(X(2).LT.0.5D0) THEN
         X2=X(2)*2.D0
         q56MX=(Q56MX+Q56MN)/2.D0
         IF(Q56MN.LE.0) GOTO 999
         Q56=Q56MN*(Q56MX/Q56MN)**X2
                      AJACOB=AJACOB*Q56*LOG(Q56MX/Q56MN)*2.d0
        ELSE
         X2=2.D0-X(2)*2.D0
         Q56MN=(Q56MX+Q56MN)/2.D0
         Q56=Q56MN+(Q56MX-Q56MN)*X2
                      AJACOB=AJACOB*(Q56MX-Q56MN)*2.d0
        END IF
       ELSE IF(IRESNS(2).EQ. 3) THEN
	thre = thresh56   
        IF(X(2).LT.thre) THEN
         X2=X(2)/thre
         Q56MX=MIN(Q56MX,(ARESNS(1,2)-3*ARESNS(2,2))**2)
         IF(Q56MX.LT.Q56MN) GOTO 999
         IF(Q56MN.LE.0)     GOTO 999
         Q56=Q56MN*(Q56MX/Q56MN)**X2
                      AJACOB=AJACOB*Q56*LOG(Q56MX/Q56MN)/thre
        ELSE
         X2=(x(2)-thre)/(1.d0-thre)
         Q56MN=MAX(Q56MN,(ARESNS(1,2)-3*ARESNS(2,2))**2)
         IF(Q56MX.LT.Q56MN) GOTO 999
         ZM    =ARESNS(1,2)
         ZM2   =ZM*ZM
         ZMG   =ARESNS(1,2)*ARESNS(2,2)
         THEMIN=ATAN((Q56MN-ZM2)/ZMG)
         THEMAX=ATAN((Q56MX-ZM2)/ZMG)
         THE   =THEMIN+(THEMAX-THEMIN)*X2
         Q56   =ZMG*TAN(THE)+ZM2
                      AJACOB=AJACOB*(THEMAX-THEMIN)
     .*                 ((Q56-ZM2)**2+ZMG**2)/ZMG/(1.d0-thre)
        END IF
       ELSE IF(IRESNS(2).EQ. 4) THEN
	thre = thresh56   
        IF(X(2).LT.thre) THEN
          X2=X(2)/thre
          Q56MX=MIN(Q56MX,(ARESNS(1,2)-3*ARESNS(2,2))**2)
          IF(Q56MX.LT.Q56MN) GOTO 999
          IF(Q56MN.LE.0)     GOTO 999
          t_min = sqrt(Q56MN)
          t_max = sqrt(Q56MX)
          t_diff = t_max - t_min
          t_val = t_min + t_diff * X2
          Q56 = t_val**2
                         AJACOB=AJACOB* 2.D0*t_val*t_diff/thre
        ELSE
          X2=(x(2)-thre)/(1.d0-thre)
          Q56MN=MAX(Q56MN,(ARESNS(1,2)-3*ARESNS(2,2))**2)
          IF(Q56MX.LT.Q56MN) GOTO 999
          ZM    =ARESNS(1,2)
          ZM2   =ZM*ZM
          ZMG   =ARESNS(1,2)*ARESNS(2,2)
          THEMIN=ATAN((Q56MN-ZM2)/ZMG)
          THEMAX=ATAN((Q56MX-ZM2)/ZMG)
          THE   =THEMIN+(THEMAX-THEMIN)*X2
          Q56   =ZMG*TAN(THE)+ZM2
                       AJACOB=AJACOB*(THEMAX-THEMIN)
     .*                  ((Q56-ZM2)**2+ZMG**2)/ZMG/(1.d0-thre)
        END IF
       ELSE IF(IRESNS(2).EQ. 5) THEN
	thre = thresh56   
        IF(X(2).LT.thre) THEN
         X2=X(2)/thre
         Q56MX=MIN(Q56MX,(ARESNS(1,2)-3*ARESNS(2,2))**2)
         IF(Q56MX.LT.Q56MN) GOTO 999
         IF(Q56MN.LE.0)     GOTO 999
         Q56=Q56MN+(Q56MX-Q56MN)*X2
                       AJACOB=AJACOB*(Q56MX-Q56MN)/thre
        ELSE
         X2=(x(2)-thre)/(1.d0-thre)
         Q56MN=MAX(Q56MN,(ARESNS(1,2)-3*ARESNS(2,2))**2)
         IF(Q56MX.LT.Q56MN) GOTO 999
         ZM    =ARESNS(1,2)
         ZM2   =ZM*ZM
         ZMG   =ARESNS(1,2)*ARESNS(2,2)
         THEMIN=ATAN((Q56MN-ZM2)/ZMG)
         THEMAX=ATAN((Q56MX-ZM2)/ZMG)
         THE   =THEMIN+(THEMAX-THEMIN)*X2
         Q56   =ZMG*TAN(THE)+ZM2
                      AJACOB=AJACOB*(THEMAX-THEMIN)
     .*                 ((Q56-ZM2)**2+ZMG**2)/ZMG/(1.d0-thre)
        END IF
       ELSE
        WRITE(6,*)' IRESNS(2) =',IRESNS(2),' is NOT supported. '
        stop
       END IF
       EX  = (QX56+AMASS2(7-I34)-Q56)/2.D0/SQRT(QX56)
       PX  = SQRT( max( EX-AMASS1(7-I34), 0.D0 ) * (EX+AMASS1(7-I34)) )
       E56 = (QX56+Q56-AMASS2(7-I34))/2.D0/SQRT(QX56)
       if (EX-AMASS1(7-I34) .LE. 0) then
C         write(6,*) '!!!Warning in K6NEM0!!!'
C         write(6,*) '  ---> EX            =', EX
C         write(6,*) '  ---> I34           =', I34
C         write(6,*) '  ---> AMASS1(7-I34) =', AMASS1(7-I34)
C         write(6,*) '  ---> Calculated PX =', PX
C         write(6,*) '  ---> This event is rejected.'
C         write(6,*) '  (User need not to worry about this.)'
         GOTO 999
       endif
***************
*    COSX     *
***************
      ipoint = 5
        X7=X(7)
        CALL LABTOK(PE(1,5-I34),PBOST1,pe(1,5-I34),PDUM,
     .             PK1,PBOST2)
	EE1=PK1(4)
        PP1=SQRT( (EE1-AMASS1(5-I34))*(EE1+AMASS1(5-I34)) )
	DE=EE1-EX
	DP=(DE*(EX+EE1)-AMASS2(5-I34)+AMASS2(7-I34))
     ./(PP1+PX)
        TTMIN=2.D0*
     . (AMASS2(7-I34)*(DE/(EX +PX )+AMASS2(7-I34)/2.D0/(EX +PX )**2)
     .- AMASS2(5-I34)*(DP/(EE1+PP1)+AMASS2(5-I34)/2.D0/(EE1+PP1)**2) )
        TTMAX=2.D0*
     . (AMASS2(7-I34)*(DE/(EX +PX )+AMASS2(7-I34)/2.D0/(EX +PX )**2)
     .- AMASS2(5-I34)*(DP/(EE1+PP1)+AMASS2(5-I34)/2.D0/(EE1+PP1)**2) )
     .+ 2.D0*2.D0*PX*PP1
       TMAX=TTMAX
       TMIN=TTMIN
       IF(TMAX.LT.TMIN) GOTO 999
       DTT=TMAX-TMIN
      IF(abs(coscut(i34-2,5-i34)).lt.cos(4.d0*rad)) THEN
       TT = TMIN+DTT*X7
      ELSE
       EPS=10D0**(-Neps_p)
       YMIN=-LOG( 1.D0+2.D0/EPS )/2
       YMAX= LOG( 1.D0+2.D0/EPS )/2
       Y=YMIN+(YMAX-YMIN)*X7
       Y7 = ((2.0D0+EPS)*EXP(2.0D0*Y) - EPS)/
     &          (EXP(2.0D0*Y) + 1.0D0 )/2.0D0
                  AJACOB=AJACOB*(YMAX-YMIN)/COSH(Y)**2*(1+EPS)/2.D0
       TT = TMIN+DTT*Y7
      END IF
                      AJACOB=AJACOB*DTT/PX/PP1/2.D0
      PPPP=(TT+AMASS2(5-I34)+AMASS2(7-I34))/2.D0
      XXXX=AMASS2(7-I34)*(DE/(EX +PX )+AMASS2(7-I34)/2.D0/(EX +PX )**2)
     .-    AMASS2(5-I34)*(DP/(EE1+PP1)+AMASS2(5-I34)/2.D0/(EE1+PP1)**2)
     .-TT/2.D0
	IF(-XXXX/PP1/PX*(2.D0+XXXX/PP1/PX).GT.0.D0) THEN
         COSX=1.D0+XXXX/PX/PP1
	 SINX=SQRT(-XXXX/PP1/PX*(2.D0+XXXX/PP1/PX))
        ELSE
         COSX=1.D0
	 SINX=0.D0
	END IF
        IF(I34.EQ.3) THEN
	 PP24=PPPP
	 TTT2=-TT
        ELSE
	 PP13=PPPP
	 TTT1=-TT
        END IF
***************
*    PHIX     *
***************
       PHIX = X4*2*PI+PHI1
                      AJACOB=AJACOB*2.D0*PI
*Particle X in px-p5-p6 CM frame
      PE(1,7-I34) = PX*SINX*COS(PHIX)
      PE(2,7-I34) = PX*SINX*SIN(PHIX)
      PE(3,7-I34) = PX*COSX
      PE(4,7-I34) = EX
*Particle 56 in px-p5-p6 CM frame
      PBOST2(1) =-PX*SINX*COS(PHIX)
      PBOST2(2) =-PX*SINX*SIN(PHIX)
      PBOST2(3) =-PX*COSX
      PBOST2(4) = E56
* boost from x-5-6 CM frame to 1-2 CM frame (Lab-frame)
        CALL KTOLAB(PE(1,5-I34),PBOST1,pe(1,7-I34),PBOST2,
     .              PE(1,7-I34),PBOST2)
         PP34=PE(4,4)*PE(4,3)-PE(3,4)*PE(3,3)
     .                       -PE(1,4)*PE(1,3)
     .                       -PE(2,4)*PE(2,3)
        Q34=AMASS2(3)+AMASS2(4)+2.D0*PP34
**********************************************************************
      IF(ISYM.EQ.1) THEN
       if(I34.eq.3 .and. pbost2(3).gt.0) goto 999
       if(I34.eq.4 .and. pbost2(3).lt.0) goto 999
      END IF
**********************************************************************
      IF     (ICOS5.EQ.2) THEN
       TT1(1)=PE(1,I34-2)-PE(1,I34)
       TT1(2)=PE(2,I34-2)-PE(2,I34)
       TT1(3)=PE(3,I34-2)-PE(3,I34)
       TT1(4)=PE(4,I34-2)-PE(4,I34)
       TT2(1)=PE(1,I34-2)-PE(1,I34)
       TT2(2)=PE(2,I34-2)-PE(2,I34)
       TT2(3)=PE(3,I34-2)-PE(3,I34)
       TT2(4)=PE(4,I34-2)-PE(4,I34)
       PK1(1)=PE(1,5-I34)-PE(1,7-I34)
       PK1(2)=PE(2,5-I34)-PE(2,7-I34)
       PK1(3)=PE(3,5-I34)-PE(3,7-I34)
       PK1(4)=PE(4,5-I34)-PE(4,7-I34)
       CALL LABTOK(TT2,PBOST2,TT1,PK1, TT1,PK1)
      ELSE IF(ICOS5.EQ.1) THEN
       TT1(1)=PE(1,5-I34)
       TT1(2)=PE(2,5-I34)
       TT1(3)=PE(3,5-I34)
       TT1(4)=PE(4,5-I34)
       TT2(1)=PE(1,5-I34)
       TT2(2)=PE(2,5-I34)
       TT2(3)=PE(3,5-I34)
       TT2(4)=PE(4,5-I34)
       CALL LABTOK(TT2,PBOST2,TT1,PDUM, TT1,PK1)
      END IF
***************
*    COS5     *
***************
      ipoint = 6
      E5      = (Q56+AMASS2(5)-AMASS2(6))/2.D0/SQRT(Q56)
      P5      = SQRT( max( E5-AMASS1(5), 0.D0 ) * (E5+AMASS1(5)) )
      E6      = (Q56+AMASS2(6)-AMASS2(5))/2.D0/SQRT(Q56)
      P6      = P5
      ITU=0
      IF(ICOS5.EQ.0) THEN
       COS5 = -1+2*X(5)
                      AJACOB=AJACOB*2
       sin5=sqrt( max( 1.d0-cos5, 0.D0 ) * (1.d0+cos5) )
      ELSE IF(ICOS5.EQ.1) THEN
        EPS=1.D-2
        YMIN=-LOG( 1.D0+2.D0/EPS )/2
        YMAX= LOG( 1.D0+2.D0/EPS )/2
        Y=YMIN+(YMAX-YMIN)*X(5)
        COS5=(1+EPS)*TANH(Y)
                  AJACOB=AJACOB*(YMAX-YMIN)/COSH(Y)**2*(1+EPS)
       sin5=sqrt( max( 1.d0-cos5, 0.D0 ) * (1.d0+cos5) )
      ELSE IF(ICOS5.EQ.2) THEN
       IF(X(5).LT.0.5D0) THEN
        X5=X(5)*2.D0
        ET1=TT1(4)
        PT1=SQRT( ET1**2+D)
        T2MIN=D-AMASS2(5)
     .+2.D0*(ET1*AMASS2(5)/(E5+P5)-P5*D/(ET1+PT1))
        T2MAX=D-AMASS2(5)+2.D0*(ET1*E5       )
        IF(T2MAX.Le.T2MIN) GOTO 999
        IF(T2MIN.LE.0) GOTO 999
        T2=T2MIN*(T2MAX/T2MIN)**X5
                     AJACOB=AJACOB*T2*LOG(T2MAX/T2MIN)/PT1/P5
        XXXX = T2+AMASS2(5)-D
     .-2.D0*(ET1*AMASS2(5)/(E5+P5)-P5*D/(ET1+PT1))
        COS5=1.D0-XXXX/P5/PT1/2.D0
        SIN5=SQRT(XXXX/P5/PT1/2.D0*(2.D0-XXXX/P5/PT1/2.D0))
        U2=-(AMASS2(5)+AMASS2(6)-D-TT-Q56+T2)
        ITU=1
       ELSE
        X5=2.D0-X(5)*2.D0
        ET1=TT1(4)
        PT1=SQRT( ET1**2+D)
        U2MIN=D-AMASS2(6)
     .+2.D0*(ET1*AMASS2(6)/(E6+P6)-P6*D/(ET1+PT1))
        U2MAX=D-AMASS2(6)+2.D0*(ET1*E6       )
        IF(U2MAX.LT.U2MIN) GOTO 999
        IF(U2MIN.LE.0) GOTO 999
        U2=U2MIN*(U2MAX/U2MIN)**X5
                     AJACOB=AJACOB*U2*LOG(U2MAX/U2MIN)/PT1/P6
        XXXX = U2+AMASS2(6)-D
     .-2.D0*(ET1*AMASS2(6)/(E6+P6)-P6*D/(ET1+PT1))
        COS6=1.D0-XXXX/P6/PT1/2.D0
        SIN5=SQRT(XXXX/P6/PT1/2.D0*(2.D0-XXXX/P6/PT1/2.D0))
        COS5=-COS6
        T2=-(AMASS2(5)+AMASS2(6)-D-TT-Q56+U2)
        ITU=2
       END IF
      ELSE
        WRITE(6,*)' ICOS5 =',ICOS5,' is NOT supported. '
        stop
      END IF
      IF(ABS(COS5).GT.1) GOTO 999
***************
*   PHI_P5    *
***************
       PHI5 = X(6)*2*PI+PHI1
                      AJACOB=AJACOB*2*PI
*Particle 5 in 5-6 CM frame
      PE(1,5) = P5*SIN5*COS(PHI5)
      PE(2,5) = P5*SIN5*SIN(PHI5)
      PE(3,5) = P5*COS5
      PE(4,5) = E5
*Particle 6 in 5-6 CM frame
      PE(1,6) =-PE(1,5)
      PE(2,6) =-PE(2,5)
      PE(3,6) =-PE(3,5)
      PE(4,6) = E6
* boost from 5-6 CM frame to 1-2 CM frame (Lab-frame)
      IF(ICOS5.le.1) THEN
       CALL WTOLAB(PE(1,5),PE(1,6),PBOST2, PE(1,5),PE(1,6))
      ELSE
       CALL KTOLAB(TT2,PBOST2,PE(1,5),PE(1,6), PE(1,5),PE(1,6))
      END IF
      ipoint = 7
      ipoint = 71
13    continue
12    continue
*Set invariants
      DO 1 I1=1,6
      DO 1 J1=1,6
1      PP(I1,J1)=PE(4,I1)*PE(4,J1)
     .-          PE(3,I1)*PE(3,J1)
     .-          PE(1,I1)*PE(1,J1)
     .-          PE(2,I1)*PE(2,J1)
      PP(1,2) = (S(3) -AMASS2(1)-AMASS2(2))/2.D0
      PP(2,1) = (S(3) -AMASS2(1)-AMASS2(2))/2.D0
      PP(1,3)=PP13
      PP(3,1)=PP13
      PP(2,4)=PP24
      PP(4,2)=PP24
      PP(3,4)=(Q34-AMASS2(3)-AMASS2(4))/2.D0
      PP(4,3)=(Q34-AMASS2(3)-AMASS2(4))/2.D0
      PP(5,6)=(Q56-AMASS2(5)-AMASS2(6))/2.D0
      PP(6,5)=(Q56-AMASS2(5)-AMASS2(6))/2.D0
      PP(7-I34,5)=(-Q56-AMASS2(7-I34)+QX56-2*PP(7-I34,6))/2.D0
      PP(5,7-I34)=(-Q56-AMASS2(7-I34)+QX56-2*PP(7-I34,6))/2.D0
      IF(I34.EQ.3) THEN
       SSS1=AMASS2(3)+Q56+2.D0*(PP(3,5)+PP(3,6))
       SSS2=QX56
      ELSE
       SSS1=QX56
       SSS2=AMASS2(4)+Q56+2.D0*(PP(4,5)+PP(4,6))
      END IF
      UUU1=AMASS2(1)+AMASS2(3)+TTT2+Q56-SSS1-TTT1
      UUU2=AMASS2(2)+AMASS2(4)+TTT1+Q56-SSS2-TTT2
      IF(ICOS5.EQ.2) THEN
      IF(I34.EQ.3) THEN
       TTTT=-T2
       UUUU=-U2
      ELSE
       TTTT=-U2
       UUUU=-T2
      END IF
      END IF
        colmbf=1.d0
      FACT = GEVPB * Flux_factor( S(3), amass2(1), amass2(2), jump )
      if (jump.EQ.1) then
        YACOB = 0
        RETURN
      endif
*Set jacobian
      YACOB = FACT*AJACOB    
     ./(2*PI)**2/(32*PI2)**3
     .*BETA(AMASS2(I34)/S(3),aX56/S(3))
     .*BETA(AMASS2(7-I34)/aX56,Q56/aX56)
     .*BETA(AMASS2(5)/Q56,AMASS2(6)/Q56)
      RETURN
999   CONTINUE
      JUMP=1
      yacob=0
      coscut(1,1)=csut11
      coscut(2,1)=csut21
      coscut(1,2)=csut12
      coscut(2,2)=csut22
      coscut(1,3)=csut13
      coscut(2,3)=csut23
      coscut(1,4)=csut14
      coscut(2,4)=csut24
      RETURN
      END
      Subroutine k6exch(i1,i2,pe,pp)
      IMPLICIT REAL* 8(A-H,O-Z)
      REAL*8  PE(4,6), PP(6,6)
      INCLUDE 'inclk.h'
*--------------------------------------------------------------------
      COMMON/CUT001_4f/COSCUT(2,4),ENGYCT(2,4),AMASCT(2,6),ARESNS(2,4)
     .,opncut,swapm2
*-----------------------------------------------------------------------
       cstemp1=coscut(1,i1-2)
       cstemp2=coscut(2,i1-2)
       coscut(1,i1-2)=coscut(1,i2-2)
       coscut(2,i1-2)=coscut(2,i2-2)
       coscut(1,i2-2)=cstemp1
       coscut(2,i2-2)=cstemp2
       egtemp1=engyct(1,i1-2)
       egtemp2=engyct(2,i1-2)
       engyct(1,i1-2)=engyct(1,i2-2)
       engyct(2,i1-2)=engyct(2,i2-2)
       engyct(1,i2-2)=egtemp1
       engyct(2,i2-2)=egtemp2
      if(i1.eq.4 .or. i1.eq.6) then
       amtemp1=amasct(1,1)
       amtemp2=amasct(2,1)
       amtemp3=amasct(1,2)
       amtemp4=amasct(2,2)
       amasct(1,1)=amasct(1,5)
       amasct(2,1)=amasct(2,5)
       amasct(1,2)=amasct(1,6)
       amasct(2,2)=amasct(2,6)
       amasct(1,5)=amtemp1
       amasct(2,5)=amtemp2
       amasct(1,6)=amtemp3
       amasct(2,6)=amtemp4
      else
       amtemp1=amasct(1,1)
       amtemp2=amasct(2,1)
       amtemp3=amasct(1,2)
       amtemp4=amasct(2,2)
       amasct(1,1)=amasct(1,6)
       amasct(2,1)=amasct(2,6)
       amasct(1,2)=amasct(1,5)
       amasct(2,2)=amasct(2,5)
       amasct(1,6)=amtemp1
       amasct(2,6)=amtemp2
       amasct(1,5)=amtemp3
       amasct(2,5)=amtemp4
      end if
*******
       if(i1.lt.i2) return
*******
       do 1 i=1,4
        ptemp=pe(i,i1)
        pe(i,i1)=pe(i,i2)
        pe(i,i2)=ptemp
1      continue
       do 2 i=1,6
        ptemp=pp(i,i1)
        pp(i,i1)=pp(i,i2)
        pp(i,i2)=ptemp
2      continue
       do 3 i=1,6
        ptemp=pp(i1,I)
        pp(i1,I)=pp(i2,I)
        pp(i2,I)=ptemp
3      continue
      return
      end
c---------------------------------------------------------------------
      subroutine kinemp_4f(ireg)
      IMPLICIT REAL* 8(A-H,O-Z)
      INCLUDE 'inclk.h'
      include 'inclp_4f.h'
      common/atest/ttt1,uuu1,ttt2,uuu2,tttt,uuuu,sss2,sss1,q56,q12,q22
     .,d1,g1,t1,u1,t2,u2
************************************************************************
      if(ireg.eq.1) then
       vn10 = TTT1
       vn20 = TTT2
       vn42 = TTTT
       vn52 = UUUU
       vn22 = SSS1
       vn14 = SSS2
       vn28 = UUU1
       vn26 = UUU2
       vn30 = Q56
       vn11 = TTT1
       vn21 = TTT2
       vn43 = TTTT
       vn53 = UUUU
       vn23 = SSS1
       vn15 = SSS2
       vn29 = UUU1
       vn27 = UUU2
       vn31 = Q56
      else
       vn10 = TTT1
       vn58 = TTT2
       vn42 = TTTT
       vn26 = UUUU
       vn56 = SSS1
       vn14 = SSS2
       vn50 = UUU1
       vn52 = UUU2
       vn48 = Q56
       vn11 = TTT1
       vn59 = TTT2
       vn43 = TTTT
       vn27 = UUUU
       vn57 = SSS1
       vn15 = SSS2
       vn51 = UUU1
       vn53 = UUU2
       vn49 = Q56
      end if
      return
      end
*-----------------------------------------------------------------------
      function beta(z1,z2)
      implicit real* 8(a-h,o-z)
      beta=sqrt( max(1-2*(z1+z2)+(z1-z2)**2,0D0) )
      return
      end
*-----------------------------------------------------------------------
      function hradi(x1,al,beta)
      implicit real* 8(a-h,o-z)
      include 'inclk.h'
      hradi =1.d0+0.75d0*beta+beta**2/4.d0*(9.d0/8.d0-pi2/3.d0)
     .+(-beta*(1.d0-x1/2.d0)
     .+beta*beta/8.d0*(-4.d0*(2.d0-x1)*log(x1)-(1.d0+3.d0*(1.d0-x1)**2
     .)/x1*log(1.d0-x1)-6+x1)
     .)*x1**(1.d0-beta)/beta      
      return
      end
*-----------------------------------------------------------------------
      function hradiO(x1,al,beta)
      implicit real* 8(a-h,o-z)
      include 'inclk.h'
      hradiO =
     .(1.d0+0.75d0*beta+beta**2/4.d0*(9.d0/8.d0-pi2/3.d0))
C     .   * beta*x1**(beta-1.d0)
     .   * beta*x1**(beta)
     .+x1*(
     .-beta*(1.d0-x1/2.d0)
     .+beta*beta/8.d0*(-4.d0*(2.d0-x1)*log(x1)-(1.d0+3.d0*(1.d0-x1)**2
     .)/x1*log(1.d0-x1)-6+x1)           
     .)
C      write(6,*) ' '
C      write(6,*) 'arguments =', x1,al,beta
C      write(6,*) (1.d0+0.75d0*beta+beta**2/4.d0*(9.d0/8.d0-pi2/3.d0))
C     .   * beta*x1**(beta-1.d0)
C      write(6,*) beta
C      write(6,*) x1
C      write(6,*) -beta*(1.d0-x1/2.d0)
C      write(6,*) beta*beta/8.d0*(-4.d0*(2.d0-x1)*log(x1))
C      write(6,*) beta*beta/8.d0*(-(1.d0+3.d0*(1.d0-x1)**2)/x1*log(1.d0-x1))
C      write(6,*) beta*beta/8.d0*(-6+x1)
      return
      end
*#######################################################################
      double precision function  Flux_factor(S, Mass2_1, Mass2_2, jump)
      implicit NONE
*****************************************************************
* S       : CM energy squared                                   *
* Mass2_* : masses squared of colliding particles               *
* jump    : error flag (if one, that point should be excluded.) *
*              (written by T.Abe on Aug. 22, 1998)              *
*****************************************************************
* ------ Argument ------
      double precision  S, Mass2_1, Mass2_2
      integer           jump
* ----------------------
* --- Local variable ---
      double precision  FF, FFF
* ----------------------
      FF  = S - Mass2_1 - Mass2_2
      FFF = (FF)**2/4.D0 - Mass2_1*Mass2_2
      if ((FF .LE. 0D0).or.(FFF .LE. 0D0)) then
        jump = 1
        Flux_factor = 0.D0
        RETURN
      endif
      Flux_factor = 4.D0 * sqrt( FFF )
      Flux_factor = 1.D0 / Flux_factor
      return
      end
*#######################################################################
