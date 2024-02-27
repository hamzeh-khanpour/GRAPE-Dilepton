      integer function  Ibtest(Ivar,bit)
      implicit NONE
*----------- Argument -----------
      integer   Ivar,bit
*--------------------------------
*-------- Local variables --------
      integer*4  Ivar4,bit4
*---------------------------------
      Ivar4 = Ivar
      bit4  = bit
      if (btest(Ivar4,bit4)) then
         Ibtest = 1
       else
         Ibtest = 0
      endif
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine  Add_bit(Ivar,bit)
      implicit NONE
* --------------- Argument ---------------
      integer   Ivar,bit
* ----------------------------------------
      Ivar = Ivar + 2**bit
      return
      end
*=====================================================================
      subroutine  Add_Graph(flag, num_flag, Igraph)
      implicit NONE
* --------------- Argument ---------------
      integer   num_flag              ! Input
      integer   flag(num_flag)        ! Input/Output
      integer   Igraph                ! Input
* ----------------------------------------
* ---------- Local variables ----------
      integer   Iflag, Igraph_local
* -------------------------------------
      Iflag = int( (Igraph-1)/30 ) + 1
      if (Iflag .GT. num_flag) then
        write(6,*) '!!!Error in Add_Graph!!!'
        write(6,*) ' ---> Iflag    =', Iflag
        write(6,*) ' ---> num_flag =', num_flag
        write(6,*) ' ---> Iflag should be <= num_flag.'
        write(6,*) ' ---> Good-bye!'
        STOP
      endif
      Igraph_local = mod(Igraph-1,30) + 1      ! from 1 to 30
      if (btest(flag(Iflag),Igraph_local-1)) then
        write(6,*) ' '
        write(6,*) '!!!Error in Add_Graph!!!'
        write(6,*) ' ---> Graph(',Igraph,') is already fired.'
        write(6,*) ' ---> This call is ignored.'
        write(6,*) ' '
      endif
      flag(Iflag) = flag(Iflag) + 2**(Igraph_local-1)
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine  Sub_bit(Ivar,bit)
      implicit NONE
* --------------- Argument ---------------
      integer   Ivar,bit
* ----------------------------------------
      Ivar = Ivar - 2**bit
      return
      end
*=====================================================================
      subroutine  Sub_Graph(flag, num_flag, Igraph)
      implicit NONE
* --------------- Argument ---------------
      integer   num_flag              ! Input
      integer   flag(num_flag)        ! Input/Output
      integer   Igraph                ! Input
* ----------------------------------------
* ---------- Local variables ----------
      integer   Iflag, Igraph_local
* -------------------------------------
      Iflag = int( (Igraph-1)/30 ) + 1
      if (Iflag .GT. num_flag) then
        write(6,*) '!!!Error in Sub_Graph!!!'
        write(6,*) ' ---> Iflag    =', Iflag
        write(6,*) ' ---> num_flag =', num_flag
        write(6,*) ' ---> Iflag should be <= num_flag.'
        write(6,*) ' ---> Good-bye!'
        STOP
      endif
      Igraph_local = mod(Igraph-1,30) + 1      ! from 1 to 30
      if (.NOT.btest(flag(Iflag),Igraph_local-1)) then
        write(6,*) ' '
        write(6,*) '!!!Error in Sub_Graph!!!'
        write(6,*) ' ---> Graph(',Igraph,') is not yet fired.'
        write(6,*) ' ---> This call is ignored.'
        write(6,*) ' '
      endif
      flag(Iflag) = flag(Iflag) - 2**(Igraph_local-1)
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      integer function  KFencode(my_code)
      implicit NONE
* --------------- Argument ---------------
      integer  my_code
* ----------------------------------------
      if      (my_code .EQ. 1) then        
         KFencode =   2
       elseif (my_code .EQ. 2) then        
         KFencode =  -2
       elseif (my_code .EQ. 3) then        
         KFencode =   1
       elseif (my_code .EQ. 4) then        
         KFencode =  -1
       elseif (my_code .EQ. 5) then        
         KFencode =   3
       elseif (my_code .EQ. 6) then        
         KFencode =  -3
       elseif (my_code .EQ. 7) then        
         KFencode =   4
       elseif (my_code .EQ. 8) then        
         KFencode =  -4
       elseif (my_code .EQ. 9) then        
         KFencode =   5
       elseif (my_code .EQ.10 .OR. my_code.EQ.0) then        
         KFencode =  -5
       elseif (my_code .EQ.11) then        
         KFencode =   6
       elseif (my_code .EQ.12) then        
         KFencode =  -6
       else
         write(6,*) '!!!Error in KFencode!!!'
         write(6,*) '  ---> unknown my_code(=',my_code,')'
         STOP
      endif
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine  START_GRAPE(istep)
      implicit NONE
* --------- Argument ---------
      integer  istep
* ----------------------------
* ------ Local variables ------
      character  c(30)*64, cside(8)*8,cleft*8,cright*8, cstep*64
      integer    i,j,iclast,icend(2)
* -----------------------------
      c( 1)='########################################################'
      c( 2)='********************************************************'
      c( 3)='++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      c( 4)='   GGGGGG     RRRRRR        A      PPPPPPP     EEEEEE   '
      c( 5)='  G      G   R      R      A A     P      P   E         '
      c( 6)='  G          R       R    A   A    P       P  E         '
      c( 7)='  G  GGGGGG  R  RRRR     A     A   PPPPPPPP   EEEEEEE   '
      c( 8)='  G     G G  R   R       AAAAAAA   P          E         '
      c( 9)='  G     G G  R    R     A       A  P          E         '
      c(10)='   GGGGG  G  R      RR  A       A  P          EEEEEEEE  '
      c(11)='                                                        '
      c(12)='  GRAce-based generator for Proton-Electron collisions  '
      c(13)='                                                        '
      c(14)='               GRAPE-Dilepton_version1.1k               '
      c(15)='                     ^^^^^^^^                           '
      c(16)='                      Mar.27, 2003                      '
      c(17)='  Comments/bug-report to Tetsuo ABE (tabe@post.kek.jp)  '
      c(18)='++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      c(19)='********************************************************'
      c(20)='########################################################'
      iclast = 20
      cside(1) = '######'
       cside(2) = '##**++'
       cside(3) = '++**##'
        cside(4) = '##****'
        cside(5) = '****##'
         cside(6) = '##**++'
         cside(7) = '++**##'
      do 110 j = 1, 64
          if (c(1)(j:j) .EQ. ' ')  then
            icend(1) = j-1
            GOTO 121
          endif
 110  continue
 121  do 120 j = 1, 8
          if (cside(1)(j:j) .EQ. ' ') then
            icend(2) = j-1
            GOTO 210
          endif
 120  continue
 210  write(6,*) ' '
      write(6,*) ' '
      do 200 i = 1, iclast
        if      ((i.eq.1).or.(i.eq.iclast)) then
           cleft  = cside(1)
           cright = cside(1)
         elseif ((i.eq.2).or.(i.eq.(iclast-1))) then
           cleft  = cside(4)
           cright = cside(5)
         elseif ((i.eq.3).or.(i.eq.(iclast-2))) then
           cleft  = cside(6)
           cright = cside(7)
         else
           cleft  = cside(2)
           cright = cside(3)
        endif
        write(6,*) '     ' // cleft(1:icend(2)) // c(i)(1:icend(1))
     &                     // cright(1:icend(2))
 200  continue
      if (istep.eq.1) then
         cstep = '<<<<<<<<<< This is an INTEGRATION step. >>>>>>>>>>'
       elseif (istep.eq.2) then
         cstep = '<<<<<<< This is an EVENT-GENERATION step. >>>>>>>'
       else
         write(6,*) '!!!Error in START_GRAPE!!!'
         write(6,*) ' ---> Invalid istep value(=',istep,' )'
         write(6,*) ' ---> Good-bye!'
         STOP
      endif
      write(6,*) ' '
      write(6,*)  '               ' // cstep
      write(6,*) ' '
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine  make_title_merge(
     &                  process, lpair, merge, Ce           ! Inputs
     &                 ,Ctitle                              ! Output
     &                            )
      implicit NONE
* --------------- Argument ---------------
      integer        process, lpair, merge        ! Input
      character*(*)  Ce                           ! Input
      character*(*)  Ctitle                       ! Output
* ----------------------------------------
* -------- Local variables --------
      character  Clpair(3)*5
* ---------------------------------
      Clpair(1) = 'e+ e-'
      Clpair(2) = 'm+ m-'
      Clpair(3) = 't+ t-'
      if (process .EQ. 3) then
        if     (merge .EQ. 12)       then        ! u+u^
          Ctitle = Ce // 'uu^'         // ' -> '
     &         //  Ce // 'q'           // ' ' // Clpair(lpair)
        elseif (merge .EQ. 34)       then        ! d+d^
          Ctitle = Ce // 'dd^'         // ' -> '
     &         //  Ce // 'q'           // ' ' // Clpair(lpair)
        elseif (merge .EQ. 56)       then        ! s+s^
          Ctitle = Ce // 'ss^'         // ' -> '
     &         //  Ce // 'q'           // ' ' // Clpair(lpair)
        elseif (merge .EQ. 78)       then        ! c+c^
          Ctitle = Ce // 'cc^'         // ' -> '
     &         //  Ce // 'q'           // ' ' // Clpair(lpair)
        elseif (merge .EQ. 90)       then        ! b+b^
          Ctitle = Ce // 'bb^'         // ' -> '
     &         //  Ce // 'q'           // ' ' // Clpair(lpair)
        elseif (merge .EQ. 1234)     then        ! u+u^+d+d^
          Ctitle = Ce // 'uu^dd^'      // ' -> '
     &         //  Ce // 'q'           // ' ' // Clpair(lpair)
        elseif (merge .EQ. 123456)   then        ! u+u^+d+d^+s+s^
          Ctitle = Ce // 'uu^dd^ss^'   // ' -> '
     &         //  Ce // 'q'           // ' ' // Clpair(lpair)
        elseif (merge .EQ. 12345678)   then      ! u+u^+d+d^+s+s^+c+c^
          Ctitle = Ce // 'uu^-cc^'     // ' -> '
     &         //  Ce // 'q'           // ' ' // Clpair(lpair)
        elseif (merge .EQ. 1234567890) then      ! u+u^+d+d^+s+s^+c+c^+b+b^
          Ctitle = Ce // 'uu^-bb^'     // ' -> '
     &         //  Ce // 'q'           // ' ' // Clpair(lpair)
        elseif (merge .EQ. 17) then              ! u+c
          Ctitle = Ce // 'uc'          // ' -> '
     &          // Ce // 'q'           // ' ' // Clpair(lpair)
        elseif (merge .EQ. 28) then              ! u^+c^
          Ctitle = Ce // 'u^c^'        // ' -> '
     &          // Ce // 'q'           // ' ' // Clpair(lpair)
        elseif (merge .EQ. 35) then              ! d+s
          Ctitle = Ce // 'ds'          // ' -> '
     &          // Ce // 'q'           // ' ' // Clpair(lpair)
        elseif (merge .EQ. 46) then              ! d^+s^
          Ctitle = Ce // 'd^s^'        // ' -> '
     &          // Ce // 'q'           // ' ' // Clpair(lpair)
        elseif (merge .EQ. 359) then             ! d+s+b
          Ctitle = Ce // 'dsb'         // ' -> '
     &          // Ce // 'q'           // ' ' // Clpair(lpair)
        elseif (merge .EQ. 460) then             ! d^+s^+b^
          Ctitle = Ce // 'd^s^b^'      // ' -> '
     &          // Ce // 'q'           // ' ' // Clpair(lpair)
        else
           write(6,*) '!!!Error in Make_Title_merge!!!'
           write(6,*) ' ---> Not-supported merge(merge=',merge,')'
           write(6,*) ' ---> Good-bye!'
           STOP
        endif
      else
        write(6,*) '!!!Error in Make_Title_merge!!!'
        write(6,*) ' ---> Not-supported process(process=',process,')'
        write(6,*) ' ---> Good-bye!'
        STOP
      endif
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine  Get_q_merged(
     &                  merge                      ! Input
     &                 ,N_merged, qflv_merged      ! Outputs
     &                            )
      implicit NONE
* --------------- Argument ---------------
      integer   merge
      integer   N_merged, qflv_merged(12)
* ----------------------------------------
* -------- Local variables --------
      integer  i,j, Isum
* ---------------------------------
      if     (merge .LT.   10) then
        N_merged = 0
        RETURN
      elseif (merge .LT.   100) then             ! 2 figures
        N_merged = 2
      elseif (merge .LT.   1000) then            ! 3 figures
        N_merged = 3
      elseif (merge .LT.   10000) then           ! 4 figures
        N_merged = 4
      elseif (merge .LT.   100000) then          ! 5 figures
        N_merged = 5
      elseif (merge .LT.   1000000) then         ! 6 figures
        N_merged = 6
      elseif (merge .LT.   10000000) then        ! 7 figures
        N_merged = 7
      elseif (merge .LT.   100000000) then       ! 8 figures
        N_merged = 8
      elseif (merge .LT.   1000000000) then      ! 9 figures
        N_merged = 9
      else                                       ! 10 figures
        N_merged = 10
      endif
      Isum = 0
      do i = 1, N_merged
        j = N_merged - i + 1
        qflv_merged(j) = mod(merge,10**i) - Isum
        Isum = Isum + qflv_merged(j)
        qflv_merged(j) = qflv_merged(j) / 10**(i-1)
      enddo
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine  Print_x_sec( Nprc, xsec, xsec_err   ! Input
     &                       )
      implicit NONE
* --------------- Argument ---------------
      integer            Nprc
      double precision   xsec(Nprc), xsec_err(Nprc)
* ----------------------------------------
* -------- Local variables --------
      double precision   XSEC_tot, XSEC_err_tot, Rorder
      integer            Iprc, Iorder
* ---------------------------------
      if (Nprc .LT. 1) then
        write(6,*) '!!!Error in Print_x_sec!!!'
        write(6,*) '  ---> Nprc:(# of processes) =', Nprc
        write(6,*) '  ---> Good-bye!'
        STOP
      endif
      XSEC_tot     = 0.D0
      XSEC_err_tot = 0.D0
      do Iprc = 1, Nprc
        if (xsec(Iprc) .LE. 0.D0) then
          write(6,*) '!!!Warning in Print_x_sec!!!'
          write(6,*) '  ---> Cross-section of process:',Iprc
          write(6,*) '       is 0 or negative value:', xsec(Iprc),'.'
        endif
        XSEC_tot     = XSEC_tot + xsec(Iprc)
        XSEC_err_tot = XSEC_err_tot + xsec_err(Iprc)**2
      enddo
      XSEC_err_tot = sqrt(XSEC_err_tot)
      Rorder  =  log10(XSEC_tot)
      Iorder  =  int( Rorder )
       if (Rorder .LT. 0.D0)  Iorder = Iorder - 1
      Rorder  =  10.D0**Iorder
      write(6,*) ' '
      write(6,100)
      write(6,300)
      if     (Iorder .LE. -100) then
        write(6,1230)  XSEC_tot/Rorder, XSEC_err_tot/Rorder, abs(Iorder)
      elseif (Iorder .LE.  -10) then
        write(6,1220)  XSEC_tot/Rorder, XSEC_err_tot/Rorder, abs(Iorder)
      elseif (Iorder .LT.    0) then
        write(6,1210)  XSEC_tot/Rorder, XSEC_err_tot/Rorder, abs(Iorder)
      elseif (Iorder .LT.   10) then
        write(6,1110)  XSEC_tot/Rorder, XSEC_err_tot/Rorder, abs(Iorder)
      elseif (Iorder .LT. 100) then
        write(6,1120)  XSEC_tot/Rorder, XSEC_err_tot/Rorder, abs(Iorder)
      else
        write(6,1130)  XSEC_tot/Rorder, XSEC_err_tot/Rorder, abs(Iorder)
      endif
      write(6,300)
      write(6,200)
      write(6,*) ' '
 100  format( '             '
     &       ,'****************'
     &       ,' << Cross Section >> '
     &       ,'****************'
     &      )
 200  format( '             '
     &       ,'*****************************************************'
     &      )
 300  format( '             '
     &       ,'*                                                   *'
     &      )
 1110 format( '             '
     &       ,'*          ( ', F8.6
     &       ,' +- ',          F8.6
     &       ,' )E+0',          I1
     &       ,' pb          *'
     &      )
 1120 format( '             '
     &       ,'*          ( ', F8.6
     &       ,' +- ',          F8.6
     &       ,' )E+',          I2
     &       ,' pb          *'
     &      )
 1130 format( '             '
     &       ,'*          ( ', F8.6
     &       ,' +- ',          F8.6
     &       ,' )E+',          I3
     &       ,' pb         *'
     &      )
 1210 format( '             '
     &       ,'*          ( ', F8.6
     &       ,' +- ',          F8.6
     &       ,' )E-0',          I1
     &       ,' pb          *'
     &      )
 1220 format( '             '
     &       ,'*          ( ', F8.6
     &       ,' +- ',          F8.6
     &       ,' )E-',          I2
     &       ,' pb          *'
     &      )
 1230 format( '             '
     &       ,'*          ( ', F8.6
     &       ,' +- ',          F8.6
     &       ,' )E-',          I3
     &       ,' pb         *'
     &      )
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function RND_GAU(mean, sigma)
      implicit NONE
* -------- Argument --------
      double precision  mean, sigma
* --------------------------
* -------- Functions --------
      double precision  PYR
       external         PYR
* ---------------------------
* ------ Local variables ------
      logical           lwhich
      double precision  two_pi, rnd1,rnd2, xxx, x1,x2
      parameter( two_pi=2D0*3.14159265359D0 )
      save  lwhich, x2
* -----------------------------
* ---------- DATA ----------
      data  lwhich/.true./
* --------------------------
      if (lwhich) then
        lwhich = .false.
        rnd1 = PYR(0)
        rnd2 = PYR(0)
        xxx = sqrt(abs(-2D0*log(rnd1)))
        x1 = xxx * cos(two_pi*rnd2) * sigma + mean
        x2 = xxx * sin(two_pi*rnd2) * sigma + mean
        RND_GAU = x1
      else
        lwhich = .true.
        RND_GAU = x2
      endif
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function  F_VMp(t)
      implicit NONE
*----------- Argument -------------
      double precision  t      
*----------------------------------
      include './inc/graepia.h'
*-------- Local variables --------
      double precision  tmp
      integer  icnt
      save     icnt
*---------------------------------
*------------- DATA -------------
      data icnt/0/
*--------------------------------
      if (t .GT. 0) then
         if (icnt .LT. 10) then
            write(6,*) '!!!Warning in F_VMp!!!'
            write(6,*) '  ---> t(=', t, ') should be negative.'
            write(6,*) '  ---> t is set to ZERO.'
         endif
         t = 0D0
         icnt = icnt + 1
         if (icnt .EQ. 10) then
            write(6,*) '  ===> Too many warnings.'
            write(6,*) '  ===> Giving up printing any more.'
         endif
      endif
      if (t_slope_VM .LE. 0) then
         tmp = 1D0 - t/0.71D0
         F_VMp  =  1D0  /  ( tmp*tmp )
       else
         F_VMp  =  exp( - dble(t_slope_VM) * abs(t) )
      endif
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function  G_VMp(x,q2)
      implicit NONE
*----------- Argument -------------
      double precision  x,q2
*----------------------------------
*-------- Local variables --------
      integer  icount
      save     icount
      double precision  UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GL
C      double precision  VAL(20)
C      character         PARM(20)*20
*---------------------------------
*------------- DATA -------------
      data icount/0/
*--------------------------------
C      if (icount .EQ. 0) then
C         PARM(1) = 'Init0'   ! Initialization of the PDFLIB common
C         VAL(1) = 0.D0
C         call PDFSET(PARM, VAL)
C         PARM(1) = 'Nptype'
C         PARM(2) = 'Ngroup'
C         PARM(3) = 'Nset'
C         VAL(1)  = 1         ! 1:proton
C         VAL(2)  = 5         ! (5,5):GRV94(LO)
C         VAL(3)  = 5         ! (4,32):CTEQ4L
C         CALL PDFSET(PARM, VAL)
C         icount = icount + 1
C      endif
      CALL STRUCTM(x,sqrt(abs(q2)),UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GL)
      G_VMp = GL      !!! x * (Gluon_density) !!!
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Z = Z * exp(i*phase*pi)
      subroutine Z_Rotate(Z, phase)
      implicit NONE
*-------- Argument --------
      complex*16        Z
      double precision  phase
*--------------------------
*------ Local variables ------
      double precision  pi, tmp
      parameter ( pi = 3.14159265359 )
      complex*16  ctmp
*-----------------------------
      ctmp = (0D0,0D0)
      if      (    phase  .EQ. 0.0D0) then
         ctmp = (1D0,0D0)
       elseif (abs(phase) .EQ. 0.5D0) then
          tmp = sign(1D0, phase)
         ctmp = cmplx(0D0,tmp)
       elseif (abs(phase) .EQ. 1.0D0) then
         ctmp = (-1D0, 0D0)
       elseif (    phase  .EQ. 1.5D0) then
         ctmp = (0D0,-1D0)
       elseif (    phase  .EQ. 2.0D0) then
         ctmp = (1D0,0D0)
       else
          tmp = phase * pi
         ctmp = cmplx(cos(tmp),sin(tmp))
      endif
      Z = Z * ctmp
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine  AMG_RESO( KF_code, am, ag )
      IMPLICIT REAL*8(A-H,O-Z)
*-------- Argument --------
      integer  KF_code                ! Input
      double precision  am, ag        ! Outputs
*--------------------------
      include 'inclm.h'
*------ Local variables ------
*-----------------------------
      if     (KF_code .EQ.  443) then
         am = amjp1s
         ag = agjp1s
      elseif (KF_code .EQ.  444) then
         am = amjp2s
         ag = agjp2s
      elseif (KF_code .EQ.  445) then
         am = amjp37
         ag = agjp37
      elseif (KF_code .EQ.  446) then
         am = amjp40
         ag = agjp40
      elseif (KF_code .EQ.  447) then
         am = amjp41
         ag = agjp41
      elseif (KF_code .EQ.  448) then
         am = amjp44
         ag = agjp44
      elseif (KF_code .EQ.  553) then
         am = amyy1s
         ag = agyy1s
      elseif (KF_code .EQ.  554) then
         am = amyy2s
         ag = agyy2s
      elseif (KF_code .EQ.  555) then
         am = amyy3s
         ag = agyy3s
      elseif (KF_code .EQ.  556) then
         am = amyy4s
         ag = agyy4s
      elseif (KF_code .EQ.  557) then
         am = amyy10
         ag = agyy10
      elseif (KF_code .EQ.  558) then
         am = amyy11
         ag = agyy11
      elseif (KF_code .EQ.   23) then
         am = amz
         ag = agz
      elseif (KF_code .EQ.   25) then
         am = amh
         ag = agh
      else
         write(6,*) '!!!Error in AMG_RESO!!!'
         write(6,*) '  ---> Unknown KF_code(=', KF_code, ')'
         write(6,*) '  ---> Good-bye!'
         STOP
      endif
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine  AGee_VM( KF_code, ag_ee )
      implicit NONE
*-------- Argument --------
      integer           KF_code      ! Input
      double precision  ag_ee        ! Output
*--------------------------
*------ Local variables ------
*-----------------------------
      if     (KF_code .EQ.  443) then
         ag_ee = 5.26D-6
      elseif (KF_code .EQ.  444) then
         ag_ee = 2.14D-6
      elseif (KF_code .EQ.  445) then
         ag_ee = 0.26D-6
      elseif (KF_code .EQ.  446) then
         ag_ee = 0.75D-6
      elseif (KF_code .EQ.  447) then
         ag_ee = 0.77D-6
      elseif (KF_code .EQ.  448) then
         ag_ee = 0.47D-6
      elseif (KF_code .EQ.  553) then
         ag_ee = 1.32D-6
      elseif (KF_code .EQ.  554) then
         ag_ee = 0.520D-6
      elseif (KF_code .EQ.  555) then
         ag_ee = 0.476D-6
      elseif (KF_code .EQ.  556) then
         ag_ee = 0.248D-6
      elseif (KF_code .EQ.  557) then
         ag_ee = 0.31D-6
      elseif (KF_code .EQ.  558) then
         ag_ee = 0.130D-6
      else
         write(6,*) '!!!Error in AGee_VM!!!'
         write(6,*) '  ---> Unknown KF_code(=', KF_code, ')'
         write(6,*) '  ---> Good-bye!'
         STOP
      endif
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      integer function  KF_from_SOPHIA(S_code_in)
      implicit NONE
*------------ Argument ------------
      integer  S_code_in
*----------------------------------
*------- Local variables -------
      integer     S_code_max
      parameter ( S_code_max = 49 )
      integer  KF_code(S_code_max)
*-------------------------------
*---------- DATA ----------
      data  KF_code/
* S_code:    1     2     3     4     5     6     7     8     9    10
*          gam    e+    e-    mu+   mu-  pi0    pi+   pi-   k+    k-
     &      22,  -11,   11,  -13,   13,  111,  211, -211,  321, -321
* S_code:   11    12    13    14    15    16    17    18    19    20
*          k0l   k0s     p     n   nue  nueb   num  numb  pbar  nbar
     &    ,130,  310, 2212, 2112,   12,  -12,   14,  -14,-2212,-2112
* S_code:   21    22    23    24    25    26    27    28    29    30
*           k0   k0b   eta   etap  rho+  rho-  rho0  k*+   k*-   k*0
     &    ,311, -311,  221,  331,  213, -213,  113,  323, -323,  313
* S_code:   31    32    33    34    35    36    37    38    39    40
*         k*0b  omeg   phi   SIG+  SIG0  SIG-  XI0   XI-   LAM  DELT++
     &   ,-313,  223,  333, 3222, 3212, 3112, 3322, 3312, 3122, 2224
* S_code:   41    42    43    44    45    46    47    48    49    50
*         DELT+ DELT0 DELT- SIG*+ SIG*0 SIG*- XI*0  XI*-  OME*-
     &   ,2214, 2114, 1114, 3224, 3214, 3114, 3324, 3314, 3334
     &             /
*--------------------------
      if (abs(S_code_in).LT.1 .OR. abs(S_code_in).GT.S_code_max) then
        write(6,*) '!!!Error in KF_from_SOPHIA!!!'
        write(6,*) '  ---> Invalid S_code_in(=', S_code_in,')'
        write(6,*) '  ---> S_code_in should be in [1,',S_code_max,']'
        write(6,*) '  ---> Good-bye!'
        STOP
      endif
      KF_from_SOPHIA = KF_code(abs(S_code_in)) *sign(1,S_code_in)
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine TAswap(A,B)
      implicit NONE
*------------ Argument ------------
      double precision  A,B
*----------------------------------
*------- Local variables -------
      double precision  Dtmp
*-------------------------------
      Dtmp = A
      A = B
      B = Dtmp
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function R_LT(x,Q2)
      implicit NONE
*---------- Argument ----------
      double precision  x,Q2      ! Input
*------------------------------
*------ Local variables ------
      double precision  Q2_loc
*-----------------------------
      Q2_loc = abs(Q2)
      if     (Q2_loc .LT. 0.05) then
        R_LT = 0.00D0
      elseif (Q2_loc .LT. 0.1) then
        R_LT = 0.05D0
      elseif (Q2_loc .LT. 0.2) then
        R_LT = 0.10D0
      elseif (Q2_loc .LT. 0.3) then
        R_LT = 0.15D0
      elseif (Q2_loc .LT. 2.0) then
        R_LT = 0.20D0
      elseif (Q2_loc .LT. 10) then
        R_LT = 0.18D0
      elseif (Q2_loc .LT. 50) then
        R_LT = 0.1D0
      else
        R_LT = 0.0D0
      endif
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
