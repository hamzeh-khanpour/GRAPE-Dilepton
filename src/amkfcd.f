*     File amkfcd.f : Sat Mar 18 19:44:57 2000
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
************************************************************************
      subroutine amkfcd
      implicit real*8(a-h,o-z)
      include 'inclm.h'
      include 'inclc.h'
*-----------------------------------------------------------------------
      do 10 i = 1, 2212
          amprtl(i) = -1.0d30
          agprtl(i) = -1.0d30
   10 continue
* mass as a function of KF-code
      amprtl(24) =    amw
      agprtl(24) =    agw
      amprtl(23) =    amz
      agprtl(23) =    agz
      amprtl(22) =    ama
      agprtl(22) =      0
      amprtl(21) =    amg
      agprtl(21) =      0
*** no KF-code for g4aux
*     amprtl(??) =    amg
*     agprtl(??) =      0
      amprtl(25) =    amh
      agprtl(25) =    agh
*** no KF-code for chi-plus
*     amprtl(??) =    amx
*     agprtl(??) =    agx
      amprtl(18) =    amy
      agprtl(18) =    agy
      amprtl(12) = amnu(1)
      agprtl(12) =      0
      amprtl(14) = amnu(2)
      agprtl(14) =      0
      amprtl(16) = amnu(3)
      agprtl(16) =      0
      amprtl(11) = amlp(1)
      agprtl(11) =      0
      amprtl(13) = amlp(2)
      agprtl(13) =      0
      amprtl(15) = amlp(3)
      agprtl(15) =      0
      amprtl( 2) = amuq(1)
      agprtl( 2) = aguq(1)
      amprtl( 4) = amuq(2)
      agprtl( 4) = aguq(2)
      amprtl( 6) = amuq(3)
      agprtl( 6) = aguq(3)
      amprtl( 1) = amdq(1)
      agprtl( 1) = agdq(1)
      amprtl( 3) = amdq(2)
      agprtl( 3) = agdq(2)
      amprtl( 5) = amdq(3)
      agprtl( 5) = agdq(3)
      amprtl(2212) =    amp
      agprtl(2212) =      0
      amprtl(2112) =    amn
      agprtl(2112) =      0
      amprtl(2212) =    amp
      agprtl(2212) =      0
*** no KF-code for c-plus
*     amprtl(??) =   amcp
*     agprtl(??) =   agcp
*** no KF-code for c-minus
*     amprtl(??) =   amcm
*     agprtl(??) =   agcm
*** no KF-code for c-z
*     amprtl(??) =   amcz
*     agprtl(??) =   agcz
*** no KF-code for c-a
*     amprtl(??) =   amca
*     agprtl(??) =      0
*** no KF-code for c-g
*     amprtl(??) =   amcg
*     agprtl(??) =      0
      amprtl(443) = amjp1s
      agprtl(443) = agjp1s
      amprtl(444) = amjp2s
      agprtl(444) = agjp2s
      amprtl(445) = amjp37
      agprtl(445) = agjp37
      amprtl(446) = amjp40
      agprtl(446) = agjp40
      amprtl(447) = amjp41
      agprtl(447) = agjp41
      amprtl(448) = amjp44
      agprtl(448) = agjp44
      amprtl(553) = amyy1s
      agprtl(553) = agyy1s
      amprtl(554) = amyy2s
      agprtl(554) = agyy2s
      amprtl(555) = amyy3s
      agprtl(555) = agyy3s
      amprtl(556) = amyy4s
      agprtl(556) = agyy4s
      amprtl(557) = amyy10
      agprtl(557) = agyy10
      amprtl(558) = amyy11
      agprtl(558) = agyy11
      return
      end
