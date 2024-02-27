*     File et_ettt/aet_etttmptbl.f : Sat Mar 18 19:45:07 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aet_etttmptbl
      implicit real*8(a-h,o-z)
      include 'inclet_ettt1.h'
      include 'incl2.h'
      include 'inclk.h'
      include 'inclet_etttp.h'
*-----------------------------------------------------------------------
      jgraph = 0

* External lines
      call smextf(lincom,amuq(3),pfet_ettt2,ptet_ettt2r,cfet_ettt2r)
      exet_ettt2r(1) = lprtcl
      call smextf(loutgo,amlp(1),pfet_ettt4,ptet_ettt4m,cfet_ettt4m)
      exet_ettt4m(1) = lantip
      call smextf(loutgo,amuq(3),pfet_ettt9,ptet_ettt9r,cfet_ettt9r)
      exet_ettt9r(1) = lprtcl
      call smextf(lincom,amlp(1),pfet_ettt17,ptet_ettt17m,cfet_ettt17m)
      exet_ettt17m(1) = lantip
      call smextf(loutgo,amlp(3),pfet_ettt33,ptet_ettt33o,cfet_ettt33o)
      exet_ettt33o(1) = lprtcl
      call smextf(lincom,amlp(3),pfet_ettt62,ptet_ettt62o,cfet_ettt62o)
      exet_ettt62o(1) = lantip


* Buffer clear for amplitude
      do 200 ih = 0 , let_etttag-1
         agc(ih) = (0.0d0,0.0d0)
  200 continue

      call aet_etttmpt0

      return
      end
************************************************************************
      subroutine aet_etttmpt0
      implicit real*8(a-h,o-z)

      include 'inclet_ettt1.h'
      include 'incl2.h'
      include 'inclk.h'
*-----------------------------------------------------------------------
      if(jselg(1) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 1
        call aet_etttg1
      endif

      if(jselg(2) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 2
        call aet_etttg2
      endif

      if(jselg(3) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 3
        call aet_etttg3
      endif

      if(jselg(4) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 4
        call aet_etttg4
      endif

      if(jselg(5) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 5
        call aet_etttg5
      endif

      if(jselg(6) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 6
        call aet_etttg6
      endif

      if(jselg(7) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 7
        call aet_etttg7
      endif

      if(jselg(8) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 8
        call aet_etttg8
      endif

      if(jselg(9) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 9
        call aet_etttg9
      endif

      if(jselg(10) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 10
        call aet_etttg10
      endif

      if(jselg(11) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 11
        call aet_etttg11
      endif

      if(jselg(12) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 12
        call aet_etttg12
      endif

      if(jhiggs .ne. 0) then
      if(jselg(13) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 13
        call aet_etttg13
      endif
      endif

      if(jhiggs .ne. 0) then
      if(jselg(14) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 14
        call aet_etttg14
      endif
      endif

      if(jselg(15) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 15
        call aet_etttg15
      endif

      if(jselg(16) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 16
        call aet_etttg16
      endif

      if(jselg(17) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 17
        call aet_etttg17
      endif

      if(jselg(18) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 18
        call aet_etttg18
      endif

      if(jhiggs .ne. 0) then
      if(jselg(19) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 19
        call aet_etttg19
      endif
      endif

      if(jhiggs .ne. 0) then
      if(jselg(20) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 20
        call aet_etttg20
      endif
      endif

      if(jhiggs .ne. 0) then
      if(jselg(21) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 21
        call aet_etttg21
      endif
      endif

      if(jhiggs .ne. 0) then
      if(jselg(22) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 22
        call aet_etttg22
      endif
      endif

      if(jhiggs .ne. 0) then
      if(igauzb .ne. 0) then
      if(jselg(23) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 23
        call aet_etttg23
      endif
      endif
      endif

      if(jselg(24) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 24
        call aet_etttg24
      endif

      if(jselg(25) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 25
        call aet_etttg25
      endif

      if(jhiggs .ne. 0) then
      if(jselg(26) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 26
        call aet_etttg26
      endif
      endif

      if(jselg(27) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 27
        call aet_etttg27
      endif

      if(jselg(28) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 28
        call aet_etttg28
      endif

      if(jhiggs .ne. 0) then
      if(jselg(29) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 29
        call aet_etttg29
      endif
      endif

      if(jselg(30) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 30
        call aet_etttg30
      endif

      if(jselg(31) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 31
        call aet_etttg31
      endif

      if(jselg(32) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 32
        call aet_etttg32
      endif

      if(jselg(33) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 33
        call aet_etttg33
      endif

      if(jhiggs .ne. 0) then
      if(jselg(34) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 34
        call aet_etttg34
      endif
      endif

      if(jhiggs .ne. 0) then
      if(jselg(35) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 35
        call aet_etttg35
      endif
      endif

      return
      end
