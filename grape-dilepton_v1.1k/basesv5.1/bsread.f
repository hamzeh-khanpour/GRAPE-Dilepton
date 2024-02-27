*   This program reads an input file with RZ format.
*       (written by T.Abe  on Oct.18 in 1998)
      subroutine BSREAD( LUN )
      implicit NONE
* -------- Argument --------
      integer  LUN
* --------------------------
* ----- Local variables -----
      character*68  NT_NAME, Htitle
      integer       NTID, LREC, IERR, icycle, iofset, i,j,k,n
     &             ,Ii4,Ir4   ! counter for integer, real*4
     &             ,Ievt_i4, Ievt_r4, Ievt_r8   ! counter for events
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
       parameter (cmax_i4 ='2047', cmax_r4 ='1023')  !!!
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
      LREC = 1024                  !!! record length
*                  ! maximum safe value of record length : 8191 (words)
      call hlimit(NWRD_HLIMIT)

* ------------ Opening a Ntuple-file ------------
      close(LUN)
      call hropen(LUN, 'bn', NT_NAME, ' ', LREC, IERR)
      if (IERR .NE. 0) then
         write(6,*) '!!!ERROR from HROPEN in BSREAD!!!'
         write(6,*) '  ---> IERR =', IERR
C         write(6,*) ' ---> Check whether input file exists or not.'
         STOP
      endif
      iofset = 0
      icycle = 999   ! to get the correct quantity of data in ex. subsequent
      call hrin(NTID+1, icycle, iofset)                  ! calls to HGNTB
      call hrin(NTID+2, icycle, iofset)
      call hrin(NTID+3, icycle, iofset)
      call hrin(NTID+4, icycle, iofset)
      call hrin(NTID+5, icycle, iofset)
* Clearing all the addresses stored by HBOOK for the specified Ntuple
C      call hbname(NTID+1, ' ', 0, '$CLEAR')
C      call hbname(NTID+2, ' ', 0, '$CLEAR')
C      call hbname(NTID+3, ' ', 0, '$CLEAR')
C      call hbnamc(NTID+4, ' ', 0, '$CLEAR')
* Definition of addresses where data of Ntuple are stored. 
      call hbname(NTID+1, 'INTEGER4',   num_i4,  '$SET')
      call hbname(NTID+2,    'REAL4',   num_r4,  '$SET')
      call hbname(NTID+3,    'REAL8',       r8,  '$SET')
      call hbnamc(NTID+4, 'CHARACTR',    char1,  '$SET')
*               ^
      call hbname(NTID+5,   'RESULT', estimate,  '$SET')
* ---------------- Reading Ntuple ----------------
      write(6,*) ' '
      write(6,*) 'Loading bases.rz...'
      Ievt_i4 = 0
      Ievt_r4 = 0
      Ievt_r8 = 0
* (BASE1)
      Ievt_i4 = Ievt_i4 + 1
      call hgnt(NTID+1, Ievt_i4, IERR)
C      Ievt_r8 = Ievt_r8 + 1
C      call hgnt(NTID+3, Ievt_r8, IERR)
        NDIM  = i4(1)
        NWILD = i4(2)
        NCALL = i4(3)
         Ii4 = 3
        do 110 i = 1, MXDIM
          Ievt_r8 = Ievt_r8 + 1
          call hgnt(NTID+3, Ievt_r8, IERR)
          XL(i) = r8
           Ii4 = Ii4 + 1
           IG(i) = i4(Ii4)
 110    continue
        do 120 i = 1, MXDIM
          Ievt_r8 = Ievt_r8 + 1
          call hgnt(NTID+3, Ievt_r8, IERR)
          XU(i) = r8
 120    continue
      write(6,*) '   ---> BASE1 : finished'

* (BASE3)
      Ievt_i4 = Ievt_i4 + 1
      call hgnt(NTID+1, Ievt_i4, IERR)
       IT     = i4(1)
      Ievt_r8 = Ievt_r8 + 1
      call hgnt(NTID+3, Ievt_r8, IERR)
       SCALLS = r8
      Ievt_r8 = Ievt_r8 + 1
      call hgnt(NTID+3, Ievt_r8, IERR)
       WGT    = r8
      Ievt_r8 = Ievt_r8 + 1
      call hgnt(NTID+3, Ievt_r8, IERR)
       TI     = r8
      Ievt_r8 = Ievt_r8 + 1
      call hgnt(NTID+3, Ievt_r8, IERR)
       TSI    = r8
      Ievt_r8 = Ievt_r8 + 1
      call hgnt(NTID+3, Ievt_r8, IERR)
       TACC   = r8
      write(6,*) '   ---> BASE3 : finished'

* (BASE4)
        do 410 i = 1, MXDIM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          Ievt_i4 = Ievt_i4 + 1
          call hgnt(NTID+1, Ievt_i4, IERR)
C          Ievt_r8 = Ievt_r8 + 1
C          call hgnt(NTID+3, Ievt_r8, IERR)
          if (i.eq.1) then
             ND  = i4(1)
             NG  = i4(2)
             NPG = i4(3)
             Ii4 = 3
           else
             Ii4 = 0
          endif
          Ii4 = Ii4 + 1
           MA(i) = i4(Ii4)
          do 420 j = 1, NDMX ! ----------
            Ievt_r8 = Ievt_r8 + 1
            call hgnt(NTID+3, Ievt_r8, IERR)
             XI(j,i) = r8
 420      continue ! ---------------------
          Ievt_r8 = Ievt_r8 + 1
          call hgnt(NTID+3, Ievt_r8, IERR)
           DX(i) = r8
 410    continue !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do 430 i = 1, LENG
          Ievt_r8 = Ievt_r8 + 1
          call hgnt(NTID+3, Ievt_r8, IERR)
           DXD(i) = r8
          Ievt_r8 = Ievt_r8 + 1
          call hgnt(NTID+3, Ievt_r8, IERR)
           DXP(i) = r8
 430    continue
      write(6,*) '   ---> BASE4 : finished'

* (BASE5)
      do 510 i = 1, ITM
        Ievt_i4 = Ievt_i4 + 1
        call hgnt(NTID+1, Ievt_i4, IERR)
        Ievt_r4 = Ievt_r4 + 1
        call hgnt(NTID+2, Ievt_r4, IERR)
C        Ievt_r8 = Ievt_r8 + 1
C        call hgnt(NTID+3, Ievt_r8, IERR)
          Ii4 = 0
          Ir4  = 0
C          Ir8  = 0
          do 520 j = 0, 1
            Ii4 = Ii4 + 1
             ITRAT(i,j) = i4(Ii4)

            Ir4  = Ir4 + 1
             TIME(i,j)  = r4(Ir4)
            Ir4  = Ir4 + 1
             EFF(i,j)   = r4(Ir4)
            Ir4  = Ir4 + 1
             WRONG(i,j) = r4(Ir4)
            Ir4  = Ir4 + 1
             TRSLT(i,j) = r4(Ir4)
            Ir4  = Ir4 + 1
             TSTD(i,j)  = r4(Ir4)
            Ir4  = Ir4 + 1
             PCNT(i,j)  = r4(Ir4)

            Ievt_r8 = Ievt_r8 + 1
            call hgnt(NTID+3, Ievt_r8, IERR)
             RESLT(i,j) = r8
            Ievt_r8 = Ievt_r8 + 1
            call hgnt(NTID+3, Ievt_r8, IERR)
             ACSTD(i,j) = r8
 520      continue
          Ir4  = Ir4 + 1
           TIME(i,2) = r4(Ir4)
 510  continue
      write(6,*) '   ---> BASE5 : finished'

* (RANDM)
      Ievt_i4 = Ievt_i4 + 1
      call hgnt(NTID+1, Ievt_i4, IERR)
      Ievt_r4 = Ievt_r4 + 1
      call hgnt(NTID+2, Ievt_r4, IERR)
        IA1 = i4( 1)
        IC1 = i4( 2)
        M1  = i4( 3)
        IX1 = i4( 4)
        IA2 = i4( 5)
        IC2 = i4( 6)
        M2  = i4( 7)
        IX2 = i4( 8)
        IA3 = i4( 9)
        IC3 = i4(10)
        M3  = i4(11)
        IX3 = i4(12)
        do 610 i = 1, 31
          RDM(i) = r4(i)
 610    continue
        RM1 = r4(32)
        RM2 = r4(33)
      write(6,*) '   ---> RANDM : finished'

* (PLOTH)
      Ievt_i4 = Ievt_i4 + 1
      call hgnt(NTID+1, Ievt_i4, IERR)
        Ii4 = 0
        do 710 i = 1, 13
          do 710 j = 1, NHS+1
            Ii4 = Ii4 + 1
             XHASH(j,i) = i4(Ii4)
 710    continue
        do 720 i = 1, 14
          do 720 j = 1, NSC+1
            Ii4 = Ii4 + 1
             DHASH(j,i) = i4(Ii4)
 720    continue
        do 730 i = 1, NHS
          Ii4 = Ii4 + 1
           IFBASE(i) = i4(Ii4)
 730    continue
        Ii4 = Ii4 + 1
 740     NHIST = i4(Ii4)
        do 750 i = 1, NHS
          do 750 j = 1, 4
            Ii4 = Ii4 + 1
             MAPL(j,i) = i4(Ii4)
 750    continue
        Ii4 = Ii4 + 1
 760     NSCAT = i4(Ii4)
        do 770 i = 1, NSC
          do 770 j = 1, 4
            Ii4 = Ii4 + 1
             MAPD(j,i) = i4(Ii4)
 770    continue
        Ii4 = Ii4 + 1
 780     NW = i4(Ii4)
      write(6,*) '   ---> PLOTH : finished'

* (PLOTB)
      if (NSCAT.gt.0) then
        write(6,*) ' At present, 2-D histograms are not supported.'
        write(6,*) ' Sorry.'
        STOP
      endif
      k = 0      ! Pointer of BUFFER
      do 810 i = 1, max(1,NHIST)
        Ievt_i4 = Ievt_i4 + 1
        call hgnt(NTID+1, Ievt_i4, IERR)
        Ievt_r4 = Ievt_r4 + 1
        call hgnt(NTID+2, Ievt_r4, IERR)
        call hgnt(NTID+4,       i, IERR)
          do 820 j = 1, num_i4-1
            k = k + 1
             IBUF(k) = i4(j)
 820      continue
          do 830 j = 1, num_r4
            k = k + 1
             BUFF(k) = r4(j)
 830      continue
          k = k + 1
           IBUF(k) = i4(num_i4)   ! = -1
          Htitle = char1//char2
           read(Htitle,840) (BUFF(n),n=k+1,k+16)
 840       format(16A4)
          k = k + 16
 810  continue
      write(6,*) '   ---> PLOTB : finished'

* (BSRSLT)
      call hgnt(NTID+5, 1, IERR)
      AVGI = estimate
      SD   = sigma
      write(6,*) '   ---> BSRSLT: finished'


      call hrend('bn')              ! Closing the direct access file.


      return
      end
