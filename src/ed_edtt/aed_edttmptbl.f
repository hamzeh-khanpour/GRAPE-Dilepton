*     File ed_edtt/aed_edttmptbl.f : Sat Mar 18 19:45:00 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aed_edttmptbl
      implicit real*8(a-h,o-z)
      include 'incled_edtt1.h'
      include 'incl2.h'
      include 'inclk.h'
      include 'incled_edttp.h'
*-----------------------------------------------------------------------
      jgraph = 0

* External lines
      call smextf(lincom,amdq(1),pfed_edtt2,pted_edtt2s,cfed_edtt2s)
      exed_edtt2s(1) = lprtcl
      call smextf(loutgo,amlp(1),pfed_edtt4,pted_edtt4m,cfed_edtt4m)
      exed_edtt4m(1) = lantip
      call smextf(loutgo,amdq(1),pfed_edtt9,pted_edtt9s,cfed_edtt9s)
      exed_edtt9s(1) = lprtcl
      call smextf(lincom,amlp(1),pfed_edtt17,pted_edtt17m,cfed_edtt17m)
      exed_edtt17m(1) = lantip
      call smextf(loutgo,amlp(3),pfed_edtt33,pted_edtt33o,cfed_edtt33o)
      exed_edtt33o(1) = lprtcl
      call smextf(lincom,amlp(3),pfed_edtt62,pted_edtt62o,cfed_edtt62o)
      exed_edtt62o(1) = lantip


* Buffer clear for amplitude
      do 200 ih = 0 , led_edttag-1
         agc(ih) = (0.0d0,0.0d0)
  200 continue

      call aed_edttmpt0

      return
      end
************************************************************************
      subroutine aed_edttmpt0
      implicit real*8(a-h,o-z)

      include 'incled_edtt1.h'
      include 'incl2.h'
      include 'inclk.h'
*-----------------------------------------------------------------------
      if(jselg(1) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 1
        call aed_edttg1
      endif

      if(jselg(2) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 2
        call aed_edttg2
      endif

      if(jselg(3) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 3
        call aed_edttg3
      endif

      if(jselg(4) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 4
        call aed_edttg4
      endif

      if(jselg(5) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 5
        call aed_edttg5
      endif

      if(jselg(6) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 6
        call aed_edttg6
      endif

      if(jselg(7) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 7
        call aed_edttg7
      endif

      if(jselg(8) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 8
        call aed_edttg8
      endif

      if(jselg(9) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 9
        call aed_edttg9
      endif

      if(jselg(10) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 10
        call aed_edttg10
      endif

      if(jselg(11) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 11
        call aed_edttg11
      endif

      if(jselg(12) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 12
        call aed_edttg12
      endif

      if(jselg(13) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 13
        call aed_edttg13
      endif

      if(jselg(14) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 14
        call aed_edttg14
      endif

      if(jselg(15) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 15
        call aed_edttg15
      endif

      if(jselg(16) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 16
        call aed_edttg16
      endif

      if(jhiggs .ne. 0) then
      if(jselg(17) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 17
        call aed_edttg17
      endif
      endif

      if(jselg(18) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 18
        call aed_edttg18
      endif

      if(jselg(19) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 19
        call aed_edttg19
      endif

      if(jselg(20) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 20
        call aed_edttg20
      endif

      if(jselg(21) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 21
        call aed_edttg21
      endif

      if(jselg(22) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 22
        call aed_edttg22
      endif

      if(jselg(23) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 23
        call aed_edttg23
      endif

      if(jselg(24) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 24
        call aed_edttg24
      endif

      if(jselg(25) .ne. 0) then
        jgraph = jgraph + 1
        igraph(jgraph) = 25
        call aed_edttg25
      endif

      return
      end