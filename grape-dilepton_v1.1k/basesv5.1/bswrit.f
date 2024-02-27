*   This program generates an output file with RZ format.
*       (written by T.Abe  on Oct.18 in 1998)
      subroutine BSWRIT( LUN )
      implicit NONE
* -------- Argument --------
      integer  LUN
* --------------------------
* ----- Local variables -----
      character*68  NT_NAME, Htitle
     &             ,NT_TITLE_i4, NT_TITLE_r4, NT_TITLE_r8, NT_TITLE_asc
     &             ,NT_TITLE_reslt
      integer       NTID, LREC, IERR, icycle, i,j,k,n, NXBIN
     &             ,Ii4,Ir4   ! counter for integer, real*4
* ---------------------------
* --------------- HBOOK common ---------------
      INTEGER    NWRD_HLIMIT
      PARAMETER (NWRD_HLIMIT = 100000)
      REAL*4       PP
      COMMON /PAWC/PP(NWRD_HLIMIT)
      integer       IQUEST(100)
      common /QUEST/IQUEST

      integer         max_i4,  max_r4
      character*(*)  cmax_i4, cmax_r4
       parameter ( max_i4 = 2047,   max_r4 = 1023 )   !!!
       parameter (cmax_i4 ='2047', cmax_r4 ='1023')   !!!
      character*32  char1, char2
      integer*4    num_i4, num_r4
      integer*4         i4(max_i4)
      real*4            r4(max_r4)
      double precision  r8
      double precision  estimate, sigma
       common /BN_CHAR/   char1, char2
       common /BN_INTE/   num_i4, i4
       common /BN_REAL4/  num_r4, r4
       common /BN_REAL8/  r8
       common /BN_RESLT/  estimate, sigma

      character*(*)  CHARACTR, INTEGER4, REAL4, REAL8, RESULT
       parameter (CHARACTR = 'char1:c*32, char2:c*32')
       parameter (INTEGER4 = 'num_i4[0,'//cmax_i4//']:i*4'
     &                    //',i4(num_i4):i*4')
       parameter (REAL4 = 'num_r4[0,'//cmax_r4//']:i*4, r4(num_r4):r*4')
       parameter (REAL8 = 'r8:r*8')
       parameter (RESULT= 'estimate:r*8, sigma:r*8')
* --------------------------------------------
      integer     MXDIM,NDMX,LENG,ITM,NHS,NSC
       PARAMETER (MXDIM = 50, NDMX = 50, LENG  = 279936 )
       PARAMETER (ITM  = 50 )
       PARAMETER (NHS = 50, NSC = 50)
*----------
      integer           NDIM,NWILD,IG,NCALL
      double precision  XL,XU
       COMMON /BASE1/ XL(MXDIM),XU(MXDIM),NDIM,NWILD
     &               ,IG(MXDIM),NCALL
*----------
      integer           IT
      double precision  SCALLS,WGT,TI,TSI,TACC
       COMMON /BASE3/ SCALLS,WGT,TI,TSI,TACC,IT
*----------
      integer           ND,NG,NPG,MA
      double precision  XI,DX,DXD,DXP
       COMMON /BASE4/ XI(NDMX,MXDIM),DX(MXDIM),DXD(LENG),DXP(LENG),
     &                ND,NG,NPG,MA(MXDIM)
*----------
      integer           ITRAT
      REAL*4            TIME, EFF, WRONG, TRSLT, TSTD, PCNT
      double precision  RESLT,ACSTD
       COMMON /BASE5/ ITRAT(ITM,0:1),TIME(ITM,0:2),EFF(ITM,0:1),
     &                WRONG(ITM,0:1),RESLT(ITM,0:1),ACSTD(ITM,0:1),
     &                TRSLT(ITM,0:1),TSTD(ITM,0:1),PCNT(ITM,0:1)
*----------
      integer           IA1,IC1,M1,IX1,IA2,IC2,M2,IX2,IA3,IC3,M3,IX3
      real              RDM,RM1,RM2
       COMMON /RANDM/ RDM(31),RM1,RM2,IA1,IC1,M1,IX1,
     &                IA2,IC2,M2,IX2,
     &                IA3,IC3,M3,IX3
*----------
      integer           NW
      INTEGER*4         XHASH,DHASH,NHIST,MAPL,IFBASE,NSCAT,MAPD
       COMMON/PLOTH/ XHASH(NHS+1,13),DHASH(NSC+1,14),IFBASE(NHS),
     &               NHIST, MAPL(4,NHS),
     &               NSCAT, MAPD(4,NSC),
     &               NW
*----------
      integer*4       IBUF
       COMMON /PLOTB/ IBUF( 281*NHS + 2527*NSC )
      real*4          BUFF( 281*NHS + 2527*NSC )
      EQUIVALENCE (IBUF(1),BUFF(1))
*----------
      integer           ITG,ITF
      real*4            STIME
      double precision  AVGI,SD,CHI2A
      COMMON /BSRSLT/AVGI,SD,CHI2A,STIME,ITG,ITF
*----------

* ------------ Initialization of HBOOK ------------
      NTID = 0                     !!! For Ntuple ID
      NT_NAME =  'bases.rz'        !!! Name of Ntuple-file
      NT_TITLE_i4    = 'BASES data(integer*4)'    !!! Name of Ntuple-title
      NT_TITLE_r4    = 'BASES data(real*4)'       !!! Name of Ntuple-title
      NT_TITLE_r8    = 'BASES data(real*8)'       !!! Name of Ntuple-title
      NT_TITLE_asc   = 'BASES data(character)'    !!! Name of Ntuple-title
      NT_TITLE_reslt = 'BASES data(result)'       !!! Name of Ntuple-title
      LREC = 1024                  !!! record length
*                  ! maximum safe value of record length : 8191 (words)

      call hlimit(NWRD_HLIMIT)

      IQUEST(10) = 1000000         !!! max. # of records(in blocks)

      call hbset('BSIZE', LREC, IERR)
      if (IERR.NE.0) then
         write(6,*) '!!!ERROR from HBSET in BSWRIT!!!'
         write(6,*) '  ---> IERR =', IERR
         STOP
      endif

* ------------ Opening a Ntuple-file ------------
      close(LUN)
      call hropen(LUN, 'bn', NT_NAME, 'NQ', LREC, IERR)
      if (IERR .NE. 0) then
         write(6,*) '!!!ERROR from HROPEN in BSWRIT!!!'
         write(6,*) '  ---> IERR =', IERR
         STOP
      endif

      call hbnt(NTID+1,   NT_TITLE_i4,   ' ')           ! Booking a CWN
      call hbnt(NTID+2,   NT_TITLE_r4,   ' ')           ! Booking a CWN
      call hbnt(NTID+3,   NT_TITLE_r8,   ' ')           ! Booking a CWN
      call hbnt(NTID+4,   NT_TITLE_asc,  ' ')           ! Booking a CWN
      call hbnt(NTID+5,   NT_TITLE_reslt,' ')           ! Booking a CWN

      call hbname(NTID+1, 'INTEGER4',   num_i4,  INTEGER4)
      call hbname(NTID+2,    'REAL4',   num_r4,  REAL4)
      call hbname(NTID+3,    'REAL8',       r8,  REAL8)
      call hbnamc(NTID+4, 'CHARACTR',    char1,  CHARACTR)
*               ^
      call hbname(NTID+5,   'RESULT', estimate,  RESULT)

* ---------------- Filling Ntuple ----------------
      write(6,*) ' '
      write(6,*) 'Making bases.rz...'
* (BASE1)
      num_i4 = MXDIM + 3
      num_r4  = 0
C      num_r8  = MXDIM * 2
        i4(1) = NDIM
        i4(2) = NWILD
        i4(3) = NCALL
         Ii4 = 3
        do 110 i = 1, MXDIM
          r8 = XL(i)
          call hfnt(NTID+3)
           Ii4 = Ii4 + 1
           i4(Ii4) = IG(i)
 110    continue
        do 120 i = 1, MXDIM
          r8 = XU(i)
          call hfnt(NTID+3)
 120    continue
      call hfnt(NTID+1)
C      call hfnt(NTID+3)
      write(6,*) '   ---> BASE1 : finished'

* (BASE3)
      num_i4 = 1
      num_r4  = 0
C      num_r8  = 5
        i4(1) = IT
        call hfnt(NTID+1)
        r8 = SCALLS
        call hfnt(NTID+3)
        r8 = WGT
        call hfnt(NTID+3)
        r8 = TI
        call hfnt(NTID+3)
        r8 = TSI
        call hfnt(NTID+3)
        r8 = TACC
        call hfnt(NTID+3)
C      call hfnt(NTID+1)
C      call hfnt(NTID+3)
      write(6,*) '   ---> BASE3 : finished'

* (BASE4)
C      num_i4 = 4
      num_r4  = 0
C      num_r8  = NDMX + 1
        do 410 i = 1, MXDIM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if (i.eq.1) then
             num_i4 = 4
             i4(1) = ND
             i4(2) = NG
             i4(3) = NPG
             Ii4 = 3
           else
             num_i4 = 1
             Ii4 = 0
          endif
          Ii4 = Ii4 + 1
           i4(Ii4) = MA(i)
          do 420 j = 1, NDMX ! ----------
            r8 = XI(j,i)
            call hfnt(NTID+3)
 420      continue ! ---------------------
          r8 = DX(i)
          call hfnt(NTID+3)
          call hfnt(NTID+1)
C          call hfnt(NTID+3)
 410    continue !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      num_i4 = 0
      num_r4  = 0
C      num_r8  = 2
        do 430 i = 1, LENG
          r8 = DXD(i)
          call hfnt(NTID+3)
          r8 = DXP(i)
          call hfnt(NTID+3)
 430    continue
      write(6,*) '   ---> BASE4 : finished'

* (BASE5)
      num_i4 =  2
      num_r4  = 13
C      num_r8  =  4
        do 510 i = 1, ITM
          Ii4 = 0
          Ir4  = 0
C          Ir8  = 0
          do 520 j = 0, 1
            Ii4 = Ii4 + 1
             i4(Ii4) = ITRAT(i,j)

            Ir4  = Ir4 + 1
             r4(Ir4)  = TIME(i,j)
            Ir4  = Ir4 + 1
             r4(Ir4)  = EFF(i,j)
            Ir4  = Ir4 + 1
             r4(Ir4)  = WRONG(i,j)
            Ir4  = Ir4 + 1
             r4(Ir4)  = TRSLT(i,j)
            Ir4  = Ir4 + 1
             r4(Ir4)  = TSTD(i,j)
            Ir4  = Ir4 + 1
             r4(Ir4)  = PCNT(i,j)

            r8 = RESLT(i,j)
            call hfnt(NTID+3)
            r8 = ACSTD(i,j)
            call hfnt(NTID+3)
 520      continue
          Ir4  = Ir4 + 1
           r4(Ir4) = TIME(i,2)
          call hfnt(NTID+1)
          call hfnt(NTID+2)
C          call hfnt(NTID+3)
 510    continue
      write(6,*) '   ---> BASE5 : finished'

* (RANDM)
      num_i4 = 12
      num_r4  = 33
C      num_r8  =  0
        i4( 1) = IA1
        i4( 2) = IC1
        i4( 3) = M1
        i4( 4) = IX1
        i4( 5) = IA2
        i4( 6) = IC2
        i4( 7) = M2
        i4( 8) = IX2
        i4( 9) = IA3
        i4(10) = IC3
        i4(11) = M3
        i4(12) = IX3
        do 610 i = 1, 31
          r4(i) = RDM(i)
 610    continue
        r4(32) = RM1
        r4(33) = RM2
      call hfnt(NTID+1)
      call hfnt(NTID+2)
      write(6,*) '   ---> RANDM : finished'

* (PLOTH)
      num_i4 = 18*(NHS+NSC)+30
      num_r4  = 0
C      num_r8  = 0
        Ii4 = 0
        do 710 i = 1, 13
          do 710 j = 1, NHS+1
            Ii4 = Ii4 + 1
             i4(Ii4) = XHASH(j,i)
 710    continue
        do 720 i = 1, 14
          do 720 j = 1, NSC+1
            Ii4 = Ii4 + 1
             i4(Ii4) = DHASH(j,i)
 720    continue
        do 730 i = 1, NHS
          Ii4 = Ii4 + 1
           i4(Ii4) = IFBASE(i)
 730    continue
        Ii4 = Ii4 + 1
 740     i4(Ii4) = NHIST
        do 750 i = 1, NHS
          do 750 j = 1, 4
            Ii4 = Ii4 + 1
             i4(Ii4) = MAPL(j,i)
 750    continue
        Ii4 = Ii4 + 1
 760     i4(Ii4) = NSCAT
        do 770 i = 1, NSC
          do 770 j = 1, 4
            Ii4 = Ii4 + 1
             i4(Ii4) = MAPD(j,i)
 770    continue
        if (NW .EQ. 0)  NW = 281
        Ii4 = Ii4 + 1
 780     i4(Ii4) = NW
      call hfnt(NTID+1)
      write(6,*) '   ---> PLOTH : finished'

* (PLOTB)
      if (NSCAT.gt.0) then
        write(6,*) ' At present, 2-D histograms are not supported.'
        write(6,*) ' Sorry.'
        STOP
      endif
      k = 0      ! Pointer of BUFFER
      do 810 i = 1, max(1,NHIST)
        NXBIN = IBUF(281*(i-1)+3)      ! # of bins
        num_i4  = NXBIN + 1 + 1
        num_r4  = (281 - num_i4) - 16
C        num_r8  = 0
        do 820 j = 1, num_i4-1
          k = k + 1
           i4(j) = IBUF(k)
 820    continue
        do 830 j = 1, num_r4
          k = k + 1
           r4(j) = BUFF(k)
 830    continue
        k = k + 1
         i4(num_i4) = IBUF(k)   ! = -1
        write(Htitle,840) (BUFF(n),n=k+1,k+16)
 840    format(16A4)
        k = k + 16
         char1 = Htitle(1:32)
         char2 = Htitle(33:64)

        call hfnt(NTID+1)
        call hfnt(NTID+2)
        call hfnt(NTID+4)
 810  continue
      write(6,*) '   ---> PLOTB : finished'

* (BSRSLT)
      estimate = AVGI
      sigma    = SD
      call hfnt(NTID+5)
      write(6,*) '   ---> BSRSLT: finished'


* ------------ Closing Ntuple-file ------------
 9999 call hldir('//bn', 'T')       ! Listing the directories of tree
                                    !  starting from the TOP directory(//bn).
      call hcdir('//bn', ' ')       ! Changing the CWD to TOP(//bn).
      call hrout(NTID+1,icycle,' ') ! Writing the histogram(Ntuple) from
                                    !  the CWD in memory onto the current
                                    !  directory on the direct access file.
      call hrout(NTID+2,icycle,' ')
      call hrout(NTID+3,icycle,' ')
      call hrout(NTID+4,icycle,' ')
      call hrout(NTID+5,icycle,' ')
      call hrend('bn')              ! Closing the direct access file.


      return
      end
