      SUBROUTINE KINIT_4f
C---------------------------------------------------------------------
C   GRACE System Library File
C   KINEM No. : 4005
C   Date      : 1995.01.25
C   Author    : Y.Kurihara
C---------------------------------------------------------------------
C   Interfaced to ep scatterings
C   Date      : 1998.02.28
C   Author    : T.Abe
C   Update on 1998.07.24
C---------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER ( MXDIM = 50 )
      COMMON / LOOP0 / LOOP
      COMMON / BPARM1 / XL(MXDIM),XU(MXDIM),NDIM,NWILD,
     &                 IG(MXDIM),NCALL
      COMMON / BPARM2 / ACC1,ACC2,ITMX1,ITMX2
      COMMON / BASE3 / SI,SI2,SWGT,SCHI,SCALLS,ATACC,NSU,IT,WGT
      CHARACTER XSTR*14
      INCLUDE 'inclk.h'      
                             
      common/kmcntl_4f/iresns(4),icos3,icosq3,icos5,isr,iswap,ident
     &                ,iphi6,ieeee,i34,itag,isym
      common/CMN_ISRPOL/isrpol
      integer          ISPRING
       common /TA_USR/ ISPRING
*-----------------------------------------------------------------------
      COMMON/KINEM1/S(3),W(3),FACT
      COMMON/CUT001_4f/COSCUT(2,4),ENGYCT(2,4),AMASCT(2,6),ARESNS(2,4)
     .,opncut,swapm2
*-----------------------------------------------------------------------
* -------------------- PDFLIB stuff --------------------
      character*20      PARM(20)
      double precision  VALUE(20)
      integer           IFLPRT
       common /W50510/  IFLPRT
      double precision  xmin,xmax, Q2min,Q2max
       common /W50513/  xmin,xmax, Q2min,Q2max
      logical           FIRST
       common /W50516/  FIRST
* ------------------------------------------------------
      include './inc/graepia.h'
      integer         jproc      
       common /amjprc/jproc
      double precision  Flux_factor
       external         Flux_factor
* ========== Initialization of parameters for kinematics ==========
      double precision  P1_lab,P2_lab, E1_lab,E2_lab  ! 4-momentum of beams
     &                 ,Pcms_lab(3),Ecms_lab(3)
     &                 ,GAMMAcms_lab(3),BETGAMcms_lab(3)
     &                 ,vec_isr(4)
       common /GEP_LAB/ P1_lab,P2_lab, E1_lab,E2_lab
     &                 ,Pcms_lab,Ecms_lab
     &                 ,GAMMAcms_lab,   BETGAMcms_lab
     &                 ,vec_isr
      write(6,*) ' '
      write(6,*) ' '
      write(6,*) '========> Start of Kinematics Initialization'
C      write(6,*) ' '
      P1_lab = dble(P_p_beam)*1.D-3   
      P2_lab = dble(P_e_beam)*1.D-3   
       E1_lab = sqrt( P1_lab**2 + amp**2 )  
       E2_lab = sqrt( P2_lab**2 + amlp(1)**2 ) 
      S(1) = amp**2 + amlp(1)**2 + 2.D0*(E1_lab*E2_lab + P1_lab*P2_lab)
       W(1) = sqrt(  max( S(1), 0.D0 )  )   
      Pcms_lab(1) = P1_lab - P2_lab         
      Ecms_lab(1) = E1_lab + E2_lab         
       GAMMAcms_lab(1)   = Ecms_lab(1)/W(1)
       BETGAMcms_lab(1)  = Pcms_lab(1)/W(1)
      do i=2,3
        S(i) = S(1)
        W(i) = W(1)
        Pcms_lab(i) = Pcms_lab(1)
        Ecms_lab(i) = Ecms_lab(1)
        GAMMAcms_lab(i)  = GAMMAcms_lab(1)
        BETGAMcms_lab(i) = BETGAMcms_lab(1)
      enddo
      write(6,*) ' '
      write(6,*) '********** Information (in Lab. frame) **********'
      write(6,*) '                              (in unit of GeV)'
      write(6,*) '  P of electrons    =',real(P2_lab)
      write(6,*) '  P of protons      =',real(P1_lab)
      write(6,*) '  Mass of electron  =',real(amlp(1))
      write(6,*) '  Mass of proton    =',real(amp)
      write(6,*) '  sqrt(S)           =',real(W(1))
      write(6,*) '  P of CMS          =',real(Pcms_lab(1))
      write(6,*) '  E of CMS          =',real(Ecms_lab(1))
      write(6,*) '  gamma of CMS      =',real(GAMMAcms_lab(1))
      write(6,*) '  beta*gamma of CMS =',real(BETGAMcms_lab(1))
      write(6,*) '***********************************************'
C      write(6,*) ' '
      totmas = amass1(3)+amass1(4)+amass1(5)+amass1(6)
      if (W(1) .LT. totmas) then
        write(6,*) '!!!Error in kinit_4f!!!'
        write(6,*) '  ---> CM energy(=', W(1), ' )'
        write(6,*) '       < Sum of final particle masses'
        write(6,*) '  ---> Good-bye!'
        STOP
      endif
      if (process.EQ.2 .OR. process.EQ.3) then
        Wmin = W_cut(1)
        Wmax = W_cut(2)
         write(6,*) ' '
         write(6,*) ' << Mass range for the hadronic system >> '
         write(6,*) '      Min. =', Wmin, ' GeV'
         write(6,*) '      Max. =', Wmax, ' GeV'
        if (Wmin.GE.Wmax) then
           write(6,*) ' !!!Error: Min. >= Max.!!! '
           write(6,*) '   ---> Good-bye!'
           STOP
        endif
        if ( Wmin .LT. 1.0799 ) then
           write(6,*) ' !!!Warning: Min. is too small!!! '
           write(6,*) '   ---> Min. is set to 1.08 GeV.'
           Wmin = 1.08
        endif
        totmas = amass1(4)+amass1(5)+amass1(6)
        if ( Wmax .GT. (W(1)-totmas) ) then
           Wmax = W(1)-totmas
           write(6,*) ' !!!Warning: Max. is too large!!! '
           write(6,*) '   ---> Max. is set to (CM energy)'
     &                   //'-(sum of lepton masses) =',real(W(1)-totmas)
        endif
C        write(6,*) ' '
      endif
      icount = 0            
      xmin  = 0.D0          
      xmax  = 1.D0          
      Q2min = 0.D0          
      Q2max = 1.D40         
* ==================== PDFLIB Initialization ====================
      write(6,*) ' '
      write(6,*) '------> PDFLIB Initialization started'
      IFLPRT = 2         
* ----- Initialization of PDFLIB common -----
      PARM(1) = 'Init0'
       VALUE(1) = 0.D0
      call PDFSET(PARM, VALUE)
* -------------------------------------------
      PARM(1) = 'Nptype'
       VALUE(1) = 1        
      PARM(2) = 'Ngroup'
       VALUE(2) = Ngroup   
      PARM(3) = 'Nset'
       VALUE(3) = Nset     
      call PDFSET(PARM, VALUE)
      write(6,*) ' '
      write(6,*) '------> PDFLIB Initialization finished'
      write(6,*) ' '
C ===============================================================
************************
* Angular cuts in CMS  *
************************
      opncut=0
      angcut=0
* particle 3
* minimum cos  cut
      COSCUT(1,1)= -cos(angcut*rad)
* maximum cos  cut
      COSCUT(2,1)=  cos(angcut*rad)
* particle 4
* minimum cos  cut
      COSCUT(1,2)= -cos(angcut*rad)
* maximum cos  cut
      COSCUT(2,2)=  cos(angcut*rad)
      angcut=10
      angcut=0
* particle 5
* minimum cos  cut
      COSCUT(1,3)= -cos(angcut*rad)
* maximum cos  cut
      COSCUT(2,3)=  cos(angcut*rad)
* particle 6
* minimum cos  cut
      COSCUT(1,4)= -cos(angcut*rad)
* maximum cos  cut
      COSCUT(2,4)=  cos(angcut*rad)
************************
* Energy  cuts in CMS  *
************************
* particle 3
* minimum energy cut
      ENGYCT(1,1) = AMASS1(3)
* maximum energy cut
      ENGYCT(2,1) = max(E1_lab,E2_lab)  
* particle 4
* minimum energy cut
      ENGYCT(1,2) = AMASS1(4)
* maximum energy cut
      ENGYCT(2,2) = max(E1_lab,E2_lab)  
* particle 5
* minimum energy cut
      ENGYCT(1,3) = AMASS1(5)
* maximum energy cut
      ENGYCT(2,3) = max(E1_lab,E2_lab)  
* particle 6
* minimum energy cut
      ENGYCT(1,4) = AMASS1(6)
* maximum energy cut
      ENGYCT(2,4) = max(E1_lab,E2_lab)  
*************************
* Cut on invariant mass *
*************************
* Paartile 3-4
* minimum
      AMASCT(1,1)=   AMASS1(3)+AMASS1(4)
* maximum
      AMASCT(2,1)= W(1)-AMASS1(5)-AMASS1(6)
* << Partile 5-6 >>
* minimum
         AMASCT(1,2) = max( dble(Mass56_cut(1)), AMASS1(5)+AMASS1(6) )
* maximum
      if (lpair .EQ. 1) then
         AMASCT(2,2) = min( dble(Mass56_cut(4)), W(1) )
       else
         AMASCT(2,2) = min( dble(Mass56_cut(2)), W(1) )
      endif
* ---> Checking upper and lower values...
      if (AMASCT(1,2) .GE. AMASCT(2,2)) then
        write(6,*) '!!!Error in kinit_4f!!!'
        write(6,*) '  ---> Lower mass cut for particle 5-6(='
     &                                      ,real(AMASCT(1,2)), ' )'
        write(6,*) '       >= upper one(=', real(AMASCT(2,2)), ' )'
        write(6,*) '  ---> Good-bye!'
        STOP
      endif
* Particle 3-5
* minimum
      AMASCT(1,3)=   AMASS1(3)+AMASS1(5)
* maximum
      AMASCT(2,3)= W(1)-AMASS1(4)-AMASS1(6)
* Particle 4-6
* minimum
      AMASCT(1,4)=   AMASS1(4)+AMASS1(6)
* maximum
      AMASCT(2,4)= W(1)-AMASS1(3)-AMASS1(5)
* Particle 3-6
* minimum
      AMASCT(1,5)=   AMASS1(3)+AMASS1(6)
* maximum
      AMASCT(2,5)= W(1)-AMASS1(4)-AMASS1(5)
* Particle 4-5
* minimum
      AMASCT(1,6)=   AMASS1(4)+AMASS1(5)
* maximum
      AMASCT(2,6)= W(1)-AMASS1(3)-AMASS1(6)
C      write(6,*) ' '
      if      (Isr_flag .EQ. 0) then
         isr = 0
         write(6,*) '------> No ISR for incoming lepton '
       elseif (Isr_flag .EQ. 1) then
         isr = 1
         write(6,*) '------> ISR for incoming lepton '
     &            //' using Structure Function method '
       elseif (Isr_flag .EQ. 2) then
         isr = 0
         write(6,*) '------> ISR for incoming lepton '
     &            //' using QED Parton Shower method '
       else
         write(6,*) ' '
         write(6,*) ' !!!Error : unknown ISR_flag(=',Isr_flag,')!!! '
         write(6,*) '   ---> Good-bye!'
         STOP
      endif
      isrpol = isr   
      write(6,*) ' '
      call k2nit
      FACT = GEVPB * Flux_factor( S(1), amp**2, amlp(1)**2, jump )
C ================= BASES RELATED INITIALIZATIONS =================
      if     (Ebeam_pol(1) .EQ. 0.) then   
         ipol_Ebeam = 0
         Lpol_Ebeam = .false.
      elseif (Ebeam_pol(2) .EQ. 0.) then   
         ipol_Ebeam = 0
         Lpol_Ebeam = .true.
      else
         ipol_Ebeam = 1
         Lpol_Ebeam = .true.
         ROT_pyrand_flag = .false.   
      endif
      if (Lpol_Ebeam) then
         if (Ebeam_pol(1).LT.-1. .OR. Ebeam_pol(1).GT.+1.) then
           write(6,*) '!!!Error in Kinem_4f!!!'
           write(6,*) '  ---> Invalid EPOL(1) =', Ebeam_pol(1)
           write(6,*) '  ---> It should be in [-1,+1].'
           write(6,*) '  ---> Good-bye!'
           STOP
         endif
         if (Ebeam_pol(2).LT.0. .OR. Ebeam_pol(2).GT.180.) then
           write(6,*) '!!!Error in Kinem_4f!!!'
           write(6,*) '  ---> Invalid EPOL(2) =', Ebeam_pol(2)
           write(6,*) '  ---> It should be in [0,180].'
           write(6,*) '  ---> Good-bye!'
           STOP
         endif
         if (Ebeam_pol(3).LT.0. .OR. Ebeam_pol(3).GT.360.) then
           write(6,*) '!!!Error in Kinem_4f!!!'
           write(6,*) '  ---> Invalid EPOL(3) =', Ebeam_pol(3)
           write(6,*) '  ---> It should be in [0,360].'
           write(6,*) '  ---> Good-bye!'
           STOP
         endif
      endif   
      if      (process .EQ. 1) then      
          NDIM =   7 + abs( isr/max(isr,1) ) +ipol_Ebeam
          NWILD=   7 + abs( isr/max(isr,1) ) -2
       elseif (process .EQ. 2) then      
          NDIM =   8 + abs( isr/max(isr,1) ) +ipol_Ebeam
          NWILD=   8 + abs( isr/max(isr,1) ) -2
       elseif (process .EQ. 3) then      
          NDIM =   8 + abs( isr/max(isr,1) ) +ipol_Ebeam
          NWILD=   8 + abs( isr/max(isr,1) ) -2
       else
          write(6,*) ' '
          write(6,*) ' !!!Error : unknown process(=',process,')!!! '
          write(6,*) '   ---> Good-bye!'
          STOP
      endif
C-----------------------------------------------------------------------
        DO 1 I=1,NDIM
         XL(I)  =   0.D0
         XU(I)  =   1.D0
         IG(I)  =   0
 1      CONTINUE
        DO 2 I=1,NWILD   
         IG(I)  =   1
 2      CONTINUE
C-----------------------------------------------------------------------
      ITMX1  =  num_it_grid
      ITMX2  =  num_it_integ
      ACC1 = dble(acc1_card*1000.)/1000D0
      ACC2 = dble(acc2_card*1000.)/1000D0
      NCALL  =  num_call
C-----------------------------------------------------------------------
      MXREG =   2      ! Two regions are allowed to use.
*--- 6. Set histograms
      if (ISPRING .NE. 1)  then
      NX = 50
      DO 100 I = 1, NDIM
         WRITE(XSTR, 110) I
  110    FORMAT('X(',I2,') SPECTRUM')
         CALL XHINIT(I, XL(I), XU(I), NX, XSTR)
  100 CONTINUE
      E1_max = max( E1_lab, W(1) )
      E2_max = max( E2_lab, W(1)/2.D0 )
      CALL XHINIT(ndim+1, 0.d0, E1_max, NX, 'Energy of Particle 3')
      CALL XHINIT(ndim+2, 0.d0, E2_max, NX, 'Energy of Particle 4')
      CALL XHINIT(ndim+3, 0.d0, E2_max, NX, 'Energy of Particle 5')
      CALL XHINIT(ndim+4, 0.d0, E2_max, NX, 'Energy of Particle 6')
      CALL XHINIT(ndim+5,-1.d0,   1.d0, NX, 'cos_the of Particle 3')
      CALL XHINIT(ndim+6,-1.d0,   1.d0, NX, 'cos_the of Particle 4')
      CALL XHINIT(ndim+7,-1.d0,   1.d0, NX, 'cos_the of Particle 5')
      CALL XHINIT(ndim+8,-1.d0,   1.d0, NX, 'cos_the of Particle 6')
      Dmax = W(1)/2D0
      CALL XHINIT(ndim+9, 0.d0, Dmax, NX, 'Mass 4-5')
      CALL XHINIT(ndim+10,0.d0, Dmax, NX, 'Mass 5-6')
      endif   ! (ISPRING .NE. 1)
      write(6,*) ' '
      write(6,*) '========> End of Kinematics Initialization'
      write(6,*) ' '
      RETURN
      END
C##########################################################################
      SUBROUTINE K2NIT
C---------------------------------------------------------------------
C   GRACE System Library File
C   KINEM No. : 4009
C   Date      : 1996.04.14  ; Modified from 4007 to fit for two-photon
C                            processes.
C             : 1996.07.23  ; To fit SF.
C   Author    : Y.Kurihara
C---------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'inclk.h'
      common/kmcntl_4f/iresns(4),icos3,icosq3,icos5,isr,iswap,ident
     &                ,iphi6,ieeee,i34,itag,isym
      data iprint/0/
*-----------------------------------------------------------------------
      COMMON/KINEM1/S(3),W(3),FACT
      COMMON/CUT001_4f/COSCUT(2,4),ENGYCT(2,4),AMASCT(2,6),ARESNS(2,4)
     .,opncut,swapm2
*-----------------------------------------------------------------------
      include './inc/graepia.h'
      integer         jproc      
       common /amjprc/jproc
*-----------------------------------------------------------------------
      cstag1=cos(0.5d0*rad)
      cstag2=cos( 4.d0*rad)
      cstag2=cos( 6.d0*rad)
      if     (coscut(2,1).ge. cstag1        .and.
     .        coscut(1,2).gt.-cstag1 ) then
       itag=1
       isym=0
       i34=3
      else if(coscut(2,1).lt. cstag1        .and.
     .        coscut(1,2).le.-cstag1 ) then
       itag=1
       isym=0
       i34=4
      else if(coscut(2,1).gt. cstag1        .and.
     .        coscut(1,2).lt.-cstag1 ) then
       itag=0
       isym=1
       i34=-999
      else
       if(coscut(2,1).lt. cstag2        .and.
     .    coscut(1,2).gt.-cstag2 ) then
        itag=2
       else
        itag=1
       end if
        isym=0
       if     (coscut(2,1).eq.-coscut(1,2)) then
        isym=1
        i34=-999
       else if(coscut(2,1).gt.-coscut(1,2)) then
        i34=3
       else
        i34=4
       end if
      end if
      isym = isym_34      !!!!!!
      i34  = ii34         !!!!!!
      write(6,*)'   itag,isym,i34 =',itag,isym,i34
      if (Neps_p.LT.1 .OR. Neps_p.GT.100) then
         write(6,*) '!!!Error in KINIT_4f!!!'
         write(6,*) '  ---> Invalid NEPSP =', Neps_p
         write(6,*) '  ---> It should be in [1,100].'
         write(6,*) '  ---> Good-bye!'
         STOP
      endif
* << Particle 5-6 >>
      IRESNS(2)  = Ireso56
      if (IRESNS(2).GE.1) then
        ARESNS(1,2)=amz      
        ARESNS(2,2)=agz      
      endif
* << Particle 4-5-6 >>
      IRESNS(3)  = Ireso456
      if (IRESNS(3).GE.1) then
        ARESNS(1,3)=Mas_reso456     ! If you want to treat narrow resonance,
        ARESNS(2,3)=Wid_reso456     ! set resonance mass and width.
      endif
*     -------
      icos3 = IcosP3   
*     -------
      icos5 = IcosP5   
*     -------
      iphi6=0
*     -------
      icosq3 = 2
*     ----------
      ident = 0
*     ---------
      write(6,*) '   Ireso56,Ireso456 =', IRESNS(2),IRESNS(3)
      if (Rnn_MJ56 .LE. 0) then
        write(6,*) '!!!Error in KINIT_4f!!!'
        write(6,*) '  ---> NNMJ56 =', Rnn_MJ56
        write(6,*) '  ---> It should be >0.'
        write(6,*) '  ---> Good-bye!'
        STOP
      endif
      Inn =  int(Rnn_MJ56) +1
      Dnn = dble(Rnn_MJ56) +1D0
      Diff_nn = Dnn - dble(Inn)
      if (Rnn_456 .LE. 1.) then
        write(6,*) '!!!Error in KINIT_4f!!!'
        write(6,*) '  ---> RNN456 =', Rnn_456
        write(6,*) '  ---> It should be >1.0.'
        write(6,*) '  ---> Good-bye!'
        STOP
      endif
*################# For Multi-Jacobian method #################
      if (IRESNS(2) .NE. 40)  GOTO 1999
*>>> QED
      weight_mj(0) = weimj_QED
      num_reso_mj = 0
*>>> Z^0
      if ( .false.
     +     .OR. jgra_sel.EQ. 4
     +     .OR. jgra_sel.EQ. 5
     +     .OR. jgra_sel.EQ.14
     +   ) then
         KF_mj = 23
         num_reso_mj = num_reso_mj + 1
         if (num_reso_mj .GT. max_reso_mj)  GOTO 1100
         call AMG_RESO( KF_mj
     &                 ,mass_mj(num_reso_mj), width_mj(num_reso_mj) )
         weight_mj(num_reso_mj)  = weimj_Z0   !!!
            if ( .false.
     +         .OR. mass_mj(num_reso_mj)  .LE.0
     +         .OR. width_mj(num_reso_mj) .LE.0
     +         .OR. weight_mj(num_reso_mj).LE.0
     +         )  GOTO 1200
         mass2_mj(num_reso_mj) = mass_mj(num_reso_mj)**2
      endif
*>>> Normalization of weights
      if (num_reso_mj .GT. 0) then
        sumsum = 0D0
        do i = 0, num_reso_mj
          sumsum = sumsum + weight_mj(i)
        enddo
        sum = 0D0
        do i = 0, num_reso_mj
          sum = sum + weight_mj(i)
          weisel_mj(i) = sum/sumsum
          weight_mj(i) = weight_mj(i)/sumsum
        enddo
        write(6,*) ' '
        write(6,*) 'MJ>> # of resonance particles =', num_reso_mj
        write(6,*) 'MJ>> Masses  ='
     &                ,(real(mass_mj(i)),i=1,num_reso_mj)
        write(6,*) 'MJ>> Widths  ='
     &                ,(real(width_mj(i)),i=1,num_reso_mj)
        write(6,*) 'MJ>> Weights ='
     &                ,(real(weight_mj(i)),i=0,num_reso_mj)
        write(6,*) 'MJ>> WeiSELs ='
     &                ,(real(weisel_mj(i)),i=0,num_reso_mj)
        write(6,*) ' '
      endif
      GOTO 1999
 1100 continue
      write(6,*) '!!!Error in KINIT_4f!!!'
      write(6,*) '  ---> num_reso_mj(=', num_reso_mj, ') exceeds'
     &            //' max_reso_mj(=', max_reso_mj, ').'
      write(6,*) '  ---> Please inform the author!'
      write(6,*) '  ---> Good-bye!'
      STOP
 1200 continue
      write(6,*) '!!!Error in KINIT_4f!!!'
      write(6,*) '  ---> Check parameters for resonance particles.'
      write(6,*) '  ---> num_reso_mj =', num_reso_mj
      write(6,*) '  ---> Mass   =', mass_mj(num_reso_mj)
      write(6,*) '  ---> Width  =', width_mj(num_reso_mj)
      write(6,*) '  ---> Weight =', weight_mj(num_reso_mj)
      write(6,*) '  ---> Good-bye!'
      STOP
 1999 continue
*#############################################################
      RETURN
      END
