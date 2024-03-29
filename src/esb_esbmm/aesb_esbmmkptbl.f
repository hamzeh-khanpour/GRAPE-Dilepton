*     File esb_esbmm/aesb_esbmmkptbl.f : Sat Mar 18 19:45:02 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aesb_esbmmkptbl
      implicit real*8(a-h,o-z)

      include 'inclesb_esbmm1.h'
      include 'inclk.h'
      include 'inclesb_esbmmp.h'
*-----------------------------------------------------------------------
      do 100 i = 1, 4
        pfesb_esbmm2(i) =  + peesb_esbmm1(i)
        pfesb_esbmm4(i) =  + peesb_esbmm2(i)
        pfesb_esbmm6(i) =  + peesb_esbmm1(i) + peesb_esbmm2(i)
        pfesb_esbmm8(i) =  - peesb_esbmm3(i)
        pfesb_esbmm10(i) =  + peesb_esbmm1(i) - peesb_esbmm3(i)
        pfesb_esbmm12(i) =  + peesb_esbmm2(i) - peesb_esbmm3(i)
        pfesb_esbmm14(i) =  + peesb_esbmm1(i) + peesb_esbmm2(i) - peesb_esbmm3(i)
        pfesb_esbmm16(i) =  - peesb_esbmm4(i)
        pfesb_esbmm18(i) =  + peesb_esbmm1(i) - peesb_esbmm4(i)
        pfesb_esbmm20(i) =  + peesb_esbmm2(i) - peesb_esbmm4(i)
        pfesb_esbmm22(i) =  + peesb_esbmm1(i) + peesb_esbmm2(i) - peesb_esbmm4(i)
        pfesb_esbmm24(i) =  - peesb_esbmm3(i) - peesb_esbmm4(i)
        pfesb_esbmm26(i) =  + peesb_esbmm1(i) - peesb_esbmm3(i) - peesb_esbmm4(i)
        pfesb_esbmm28(i) =  + peesb_esbmm2(i) - peesb_esbmm3(i) - peesb_esbmm4(i)
        pfesb_esbmm30(i) =  + peesb_esbmm5(i) + peesb_esbmm6(i)
        pfesb_esbmm32(i) =  - peesb_esbmm5(i)
        pfesb_esbmm34(i) =  + peesb_esbmm1(i) - peesb_esbmm5(i)
        pfesb_esbmm36(i) =  + peesb_esbmm2(i) - peesb_esbmm5(i)
        pfesb_esbmm38(i) =  + peesb_esbmm1(i) + peesb_esbmm2(i) - peesb_esbmm5(i)
        pfesb_esbmm40(i) =  - peesb_esbmm3(i) - peesb_esbmm5(i)
        pfesb_esbmm42(i) =  + peesb_esbmm1(i) - peesb_esbmm3(i) - peesb_esbmm5(i)
        pfesb_esbmm44(i) =  + peesb_esbmm2(i) - peesb_esbmm3(i) - peesb_esbmm5(i)
        pfesb_esbmm46(i) =  + peesb_esbmm4(i) + peesb_esbmm6(i)
        pfesb_esbmm48(i) =  - peesb_esbmm4(i) - peesb_esbmm5(i)
        pfesb_esbmm50(i) =  + peesb_esbmm1(i) - peesb_esbmm4(i) - peesb_esbmm5(i)
        pfesb_esbmm52(i) =  + peesb_esbmm2(i) - peesb_esbmm4(i) - peesb_esbmm5(i)
        pfesb_esbmm54(i) =  + peesb_esbmm3(i) + peesb_esbmm6(i)
        pfesb_esbmm56(i) =  - peesb_esbmm3(i) - peesb_esbmm4(i) - peesb_esbmm5(i)
        pfesb_esbmm58(i) =  - peesb_esbmm2(i) + peesb_esbmm6(i)
        pfesb_esbmm60(i) =  - peesb_esbmm1(i) + peesb_esbmm6(i)
        pfesb_esbmm62(i) =  + peesb_esbmm6(i)

        pfesb_esbmm3(i) = - pfesb_esbmm2(i)
        pfesb_esbmm5(i) = - pfesb_esbmm4(i)
        pfesb_esbmm7(i) = - pfesb_esbmm6(i)
        pfesb_esbmm9(i) = - pfesb_esbmm8(i)
        pfesb_esbmm11(i) = - pfesb_esbmm10(i)
        pfesb_esbmm13(i) = - pfesb_esbmm12(i)
        pfesb_esbmm15(i) = - pfesb_esbmm14(i)
        pfesb_esbmm17(i) = - pfesb_esbmm16(i)
        pfesb_esbmm19(i) = - pfesb_esbmm18(i)
        pfesb_esbmm21(i) = - pfesb_esbmm20(i)
        pfesb_esbmm23(i) = - pfesb_esbmm22(i)
        pfesb_esbmm25(i) = - pfesb_esbmm24(i)
        pfesb_esbmm27(i) = - pfesb_esbmm26(i)
        pfesb_esbmm29(i) = - pfesb_esbmm28(i)
        pfesb_esbmm31(i) = - pfesb_esbmm30(i)
        pfesb_esbmm33(i) = - pfesb_esbmm32(i)
        pfesb_esbmm35(i) = - pfesb_esbmm34(i)
        pfesb_esbmm37(i) = - pfesb_esbmm36(i)
        pfesb_esbmm39(i) = - pfesb_esbmm38(i)
        pfesb_esbmm41(i) = - pfesb_esbmm40(i)
        pfesb_esbmm43(i) = - pfesb_esbmm42(i)
        pfesb_esbmm45(i) = - pfesb_esbmm44(i)
        pfesb_esbmm47(i) = - pfesb_esbmm46(i)
        pfesb_esbmm49(i) = - pfesb_esbmm48(i)
        pfesb_esbmm51(i) = - pfesb_esbmm50(i)
        pfesb_esbmm53(i) = - pfesb_esbmm52(i)
        pfesb_esbmm55(i) = - pfesb_esbmm54(i)
        pfesb_esbmm57(i) = - pfesb_esbmm56(i)
        pfesb_esbmm59(i) = - pfesb_esbmm58(i)
        pfesb_esbmm61(i) = - pfesb_esbmm60(i)
        pfesb_esbmm63(i) = - pfesb_esbmm62(i)
  100 continue

      vnesb_esbmm2 =  + amass2(1)
      vnesb_esbmm4 =  + amass2(2)
      vnesb_esbmm6 =  + 2.0d0*prod(1,2) + amass2(1) + amass2(2)
      vnesb_esbmm8 =  + amass2(3)
      vnesb_esbmm10 =  - 2.0d0*prod(1,3) + amass2(1) + amass2(3)
      vnesb_esbmm12 =  - 2.0d0*prod(2,3) + amass2(2) + amass2(3)
      vnesb_esbmm14 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,3) - 2.0d0*
     &      prod(2,3) + amass2(1) + amass2(2) + amass2(3)
      vnesb_esbmm16 =  + amass2(4)
      vnesb_esbmm18 =  - 2.0d0*prod(1,4) + amass2(1) + amass2(4)
      vnesb_esbmm20 =  - 2.0d0*prod(2,4) + amass2(2) + amass2(4)
      vnesb_esbmm22 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,4) - 2.0d0*
     &      prod(2,4) + amass2(1) + amass2(2) + amass2(4)
      vnesb_esbmm24 =  + 2.0d0*prod(3,4) + amass2(3) + amass2(4)
      vnesb_esbmm26 =  - 2.0d0*prod(1,3) - 2.0d0*prod(1,4) + 2.0d0*
     &      prod(3,4) + amass2(1) + amass2(3) + amass2(4)
      vnesb_esbmm28 =  - 2.0d0*prod(2,3) - 2.0d0*prod(2,4) + 2.0d0*
     &      prod(3,4) + amass2(2) + amass2(3) + amass2(4)
      vnesb_esbmm30 =  + 2.0d0*prod(5,6) + amass2(5) + amass2(6)
      vnesb_esbmm32 =  + amass2(5)
      vnesb_esbmm34 =  - 2.0d0*prod(1,5) + amass2(1) + amass2(5)
      vnesb_esbmm36 =  - 2.0d0*prod(2,5) + amass2(2) + amass2(5)
      vnesb_esbmm38 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,5) - 2.0d0*
     &      prod(2,5) + amass2(1) + amass2(2) + amass2(5)
      vnesb_esbmm40 =  + 2.0d0*prod(3,5) + amass2(3) + amass2(5)
      vnesb_esbmm42 =  - 2.0d0*prod(1,3) - 2.0d0*prod(1,5) + 2.0d0*
     &      prod(3,5) + amass2(1) + amass2(3) + amass2(5)
      vnesb_esbmm44 =  - 2.0d0*prod(2,3) - 2.0d0*prod(2,5) + 2.0d0*
     &      prod(3,5) + amass2(2) + amass2(3) + amass2(5)
      vnesb_esbmm46 =  + 2.0d0*prod(4,6) + amass2(4) + amass2(6)
      vnesb_esbmm48 =  + 2.0d0*prod(4,5) + amass2(4) + amass2(5)
      vnesb_esbmm50 =  - 2.0d0*prod(1,4) - 2.0d0*prod(1,5) + 2.0d0*
     &      prod(4,5) + amass2(1) + amass2(4) + amass2(5)
      vnesb_esbmm52 =  - 2.0d0*prod(2,4) - 2.0d0*prod(2,5) + 2.0d0*
     &      prod(4,5) + amass2(2) + amass2(4) + amass2(5)
      vnesb_esbmm54 =  + 2.0d0*prod(3,6) + amass2(3) + amass2(6)
      vnesb_esbmm56 =  + 2.0d0*prod(3,4) + 2.0d0*prod(3,5) + 2.0d0*
     &      prod(4,5) + amass2(3) + amass2(4) + amass2(5)
      vnesb_esbmm58 =  - 2.0d0*prod(2,6) + amass2(2) + amass2(6)
      vnesb_esbmm60 =  - 2.0d0*prod(1,6) + amass2(1) + amass2(6)
      vnesb_esbmm62 =  + amass2(6)
      vnesb_esbmm3 = vnesb_esbmm2
      vnesb_esbmm5 = vnesb_esbmm4
      vnesb_esbmm7 = vnesb_esbmm6
      vnesb_esbmm9 = vnesb_esbmm8
      vnesb_esbmm11 = vnesb_esbmm10
      vnesb_esbmm13 = vnesb_esbmm12
      vnesb_esbmm15 = vnesb_esbmm14
      vnesb_esbmm17 = vnesb_esbmm16
      vnesb_esbmm19 = vnesb_esbmm18
      vnesb_esbmm21 = vnesb_esbmm20
      vnesb_esbmm23 = vnesb_esbmm22
      vnesb_esbmm25 = vnesb_esbmm24
      vnesb_esbmm27 = vnesb_esbmm26
      vnesb_esbmm29 = vnesb_esbmm28
      vnesb_esbmm31 = vnesb_esbmm30
      vnesb_esbmm33 = vnesb_esbmm32
      vnesb_esbmm35 = vnesb_esbmm34
      vnesb_esbmm37 = vnesb_esbmm36
      vnesb_esbmm39 = vnesb_esbmm38
      vnesb_esbmm41 = vnesb_esbmm40
      vnesb_esbmm43 = vnesb_esbmm42
      vnesb_esbmm45 = vnesb_esbmm44
      vnesb_esbmm47 = vnesb_esbmm46
      vnesb_esbmm49 = vnesb_esbmm48
      vnesb_esbmm51 = vnesb_esbmm50
      vnesb_esbmm53 = vnesb_esbmm52
      vnesb_esbmm55 = vnesb_esbmm54
      vnesb_esbmm57 = vnesb_esbmm56
      vnesb_esbmm59 = vnesb_esbmm58
      vnesb_esbmm61 = vnesb_esbmm60
      vnesb_esbmm63 = vnesb_esbmm62
      return
      end
