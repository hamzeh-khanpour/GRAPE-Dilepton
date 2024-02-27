*     File ecb_ecbee/aecb_ecbeekptbl.f : Sat Mar 18 19:45:04 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aecb_ecbeekptbl
      implicit real*8(a-h,o-z)

      include 'inclecb_ecbee1.h'
      include 'inclk.h'
      include 'inclecb_ecbeep.h'
*-----------------------------------------------------------------------
      do 100 i = 1, 4
        pfecb_ecbee2(i) =  + peecb_ecbee1(i)
        pfecb_ecbee4(i) =  + peecb_ecbee2(i)
        pfecb_ecbee6(i) =  + peecb_ecbee1(i) + peecb_ecbee2(i)
        pfecb_ecbee8(i) =  - peecb_ecbee3(i)
        pfecb_ecbee10(i) =  + peecb_ecbee1(i) - peecb_ecbee3(i)
        pfecb_ecbee12(i) =  + peecb_ecbee2(i) - peecb_ecbee3(i)
        pfecb_ecbee14(i) =  + peecb_ecbee1(i) + peecb_ecbee2(i) - peecb_ecbee3(i)
        pfecb_ecbee16(i) =  - peecb_ecbee4(i)
        pfecb_ecbee18(i) =  + peecb_ecbee1(i) - peecb_ecbee4(i)
        pfecb_ecbee20(i) =  + peecb_ecbee2(i) - peecb_ecbee4(i)
        pfecb_ecbee22(i) =  + peecb_ecbee1(i) + peecb_ecbee2(i) - peecb_ecbee4(i)
        pfecb_ecbee24(i) =  - peecb_ecbee3(i) - peecb_ecbee4(i)
        pfecb_ecbee26(i) =  + peecb_ecbee1(i) - peecb_ecbee3(i) - peecb_ecbee4(i)
        pfecb_ecbee28(i) =  + peecb_ecbee2(i) - peecb_ecbee3(i) - peecb_ecbee4(i)
        pfecb_ecbee30(i) =  + peecb_ecbee5(i) + peecb_ecbee6(i)
        pfecb_ecbee32(i) =  - peecb_ecbee5(i)
        pfecb_ecbee34(i) =  + peecb_ecbee1(i) - peecb_ecbee5(i)
        pfecb_ecbee36(i) =  + peecb_ecbee2(i) - peecb_ecbee5(i)
        pfecb_ecbee38(i) =  + peecb_ecbee1(i) + peecb_ecbee2(i) - peecb_ecbee5(i)
        pfecb_ecbee40(i) =  - peecb_ecbee3(i) - peecb_ecbee5(i)
        pfecb_ecbee42(i) =  + peecb_ecbee1(i) - peecb_ecbee3(i) - peecb_ecbee5(i)
        pfecb_ecbee44(i) =  + peecb_ecbee2(i) - peecb_ecbee3(i) - peecb_ecbee5(i)
        pfecb_ecbee46(i) =  + peecb_ecbee4(i) + peecb_ecbee6(i)
        pfecb_ecbee48(i) =  - peecb_ecbee4(i) - peecb_ecbee5(i)
        pfecb_ecbee50(i) =  + peecb_ecbee1(i) - peecb_ecbee4(i) - peecb_ecbee5(i)
        pfecb_ecbee52(i) =  + peecb_ecbee2(i) - peecb_ecbee4(i) - peecb_ecbee5(i)
        pfecb_ecbee54(i) =  + peecb_ecbee3(i) + peecb_ecbee6(i)
        pfecb_ecbee56(i) =  - peecb_ecbee3(i) - peecb_ecbee4(i) - peecb_ecbee5(i)
        pfecb_ecbee58(i) =  - peecb_ecbee2(i) + peecb_ecbee6(i)
        pfecb_ecbee60(i) =  - peecb_ecbee1(i) + peecb_ecbee6(i)
        pfecb_ecbee62(i) =  + peecb_ecbee6(i)

        pfecb_ecbee3(i) = - pfecb_ecbee2(i)
        pfecb_ecbee5(i) = - pfecb_ecbee4(i)
        pfecb_ecbee7(i) = - pfecb_ecbee6(i)
        pfecb_ecbee9(i) = - pfecb_ecbee8(i)
        pfecb_ecbee11(i) = - pfecb_ecbee10(i)
        pfecb_ecbee13(i) = - pfecb_ecbee12(i)
        pfecb_ecbee15(i) = - pfecb_ecbee14(i)
        pfecb_ecbee17(i) = - pfecb_ecbee16(i)
        pfecb_ecbee19(i) = - pfecb_ecbee18(i)
        pfecb_ecbee21(i) = - pfecb_ecbee20(i)
        pfecb_ecbee23(i) = - pfecb_ecbee22(i)
        pfecb_ecbee25(i) = - pfecb_ecbee24(i)
        pfecb_ecbee27(i) = - pfecb_ecbee26(i)
        pfecb_ecbee29(i) = - pfecb_ecbee28(i)
        pfecb_ecbee31(i) = - pfecb_ecbee30(i)
        pfecb_ecbee33(i) = - pfecb_ecbee32(i)
        pfecb_ecbee35(i) = - pfecb_ecbee34(i)
        pfecb_ecbee37(i) = - pfecb_ecbee36(i)
        pfecb_ecbee39(i) = - pfecb_ecbee38(i)
        pfecb_ecbee41(i) = - pfecb_ecbee40(i)
        pfecb_ecbee43(i) = - pfecb_ecbee42(i)
        pfecb_ecbee45(i) = - pfecb_ecbee44(i)
        pfecb_ecbee47(i) = - pfecb_ecbee46(i)
        pfecb_ecbee49(i) = - pfecb_ecbee48(i)
        pfecb_ecbee51(i) = - pfecb_ecbee50(i)
        pfecb_ecbee53(i) = - pfecb_ecbee52(i)
        pfecb_ecbee55(i) = - pfecb_ecbee54(i)
        pfecb_ecbee57(i) = - pfecb_ecbee56(i)
        pfecb_ecbee59(i) = - pfecb_ecbee58(i)
        pfecb_ecbee61(i) = - pfecb_ecbee60(i)
        pfecb_ecbee63(i) = - pfecb_ecbee62(i)
  100 continue

      vnecb_ecbee2 =  + amass2(1)
      vnecb_ecbee4 =  + amass2(2)
      vnecb_ecbee6 =  + 2.0d0*prod(1,2) + amass2(1) + amass2(2)
      vnecb_ecbee8 =  + amass2(3)
      vnecb_ecbee10 =  - 2.0d0*prod(1,3) + amass2(1) + amass2(3)
      vnecb_ecbee12 =  - 2.0d0*prod(2,3) + amass2(2) + amass2(3)
      vnecb_ecbee14 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,3) - 2.0d0*
     &      prod(2,3) + amass2(1) + amass2(2) + amass2(3)
      vnecb_ecbee16 =  + amass2(4)
      vnecb_ecbee18 =  - 2.0d0*prod(1,4) + amass2(1) + amass2(4)
      vnecb_ecbee20 =  - 2.0d0*prod(2,4) + amass2(2) + amass2(4)
      vnecb_ecbee22 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,4) - 2.0d0*
     &      prod(2,4) + amass2(1) + amass2(2) + amass2(4)
      vnecb_ecbee24 =  + 2.0d0*prod(3,4) + amass2(3) + amass2(4)
      vnecb_ecbee26 =  - 2.0d0*prod(1,3) - 2.0d0*prod(1,4) + 2.0d0*
     &      prod(3,4) + amass2(1) + amass2(3) + amass2(4)
      vnecb_ecbee28 =  - 2.0d0*prod(2,3) - 2.0d0*prod(2,4) + 2.0d0*
     &      prod(3,4) + amass2(2) + amass2(3) + amass2(4)
      vnecb_ecbee30 =  + 2.0d0*prod(5,6) + amass2(5) + amass2(6)
      vnecb_ecbee32 =  + amass2(5)
      vnecb_ecbee34 =  - 2.0d0*prod(1,5) + amass2(1) + amass2(5)
      vnecb_ecbee36 =  - 2.0d0*prod(2,5) + amass2(2) + amass2(5)
      vnecb_ecbee38 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,5) - 2.0d0*
     &      prod(2,5) + amass2(1) + amass2(2) + amass2(5)
      vnecb_ecbee40 =  + 2.0d0*prod(3,5) + amass2(3) + amass2(5)
      vnecb_ecbee42 =  - 2.0d0*prod(1,3) - 2.0d0*prod(1,5) + 2.0d0*
     &      prod(3,5) + amass2(1) + amass2(3) + amass2(5)
      vnecb_ecbee44 =  - 2.0d0*prod(2,3) - 2.0d0*prod(2,5) + 2.0d0*
     &      prod(3,5) + amass2(2) + amass2(3) + amass2(5)
      vnecb_ecbee46 =  + 2.0d0*prod(4,6) + amass2(4) + amass2(6)
      vnecb_ecbee48 =  + 2.0d0*prod(4,5) + amass2(4) + amass2(5)
      vnecb_ecbee50 =  - 2.0d0*prod(1,4) - 2.0d0*prod(1,5) + 2.0d0*
     &      prod(4,5) + amass2(1) + amass2(4) + amass2(5)
      vnecb_ecbee52 =  - 2.0d0*prod(2,4) - 2.0d0*prod(2,5) + 2.0d0*
     &      prod(4,5) + amass2(2) + amass2(4) + amass2(5)
      vnecb_ecbee54 =  + 2.0d0*prod(3,6) + amass2(3) + amass2(6)
      vnecb_ecbee56 =  + 2.0d0*prod(3,4) + 2.0d0*prod(3,5) + 2.0d0*
     &      prod(4,5) + amass2(3) + amass2(4) + amass2(5)
      vnecb_ecbee58 =  - 2.0d0*prod(2,6) + amass2(2) + amass2(6)
      vnecb_ecbee60 =  - 2.0d0*prod(1,6) + amass2(1) + amass2(6)
      vnecb_ecbee62 =  + amass2(6)
      vnecb_ecbee3 = vnecb_ecbee2
      vnecb_ecbee5 = vnecb_ecbee4
      vnecb_ecbee7 = vnecb_ecbee6
      vnecb_ecbee9 = vnecb_ecbee8
      vnecb_ecbee11 = vnecb_ecbee10
      vnecb_ecbee13 = vnecb_ecbee12
      vnecb_ecbee15 = vnecb_ecbee14
      vnecb_ecbee17 = vnecb_ecbee16
      vnecb_ecbee19 = vnecb_ecbee18
      vnecb_ecbee21 = vnecb_ecbee20
      vnecb_ecbee23 = vnecb_ecbee22
      vnecb_ecbee25 = vnecb_ecbee24
      vnecb_ecbee27 = vnecb_ecbee26
      vnecb_ecbee29 = vnecb_ecbee28
      vnecb_ecbee31 = vnecb_ecbee30
      vnecb_ecbee33 = vnecb_ecbee32
      vnecb_ecbee35 = vnecb_ecbee34
      vnecb_ecbee37 = vnecb_ecbee36
      vnecb_ecbee39 = vnecb_ecbee38
      vnecb_ecbee41 = vnecb_ecbee40
      vnecb_ecbee43 = vnecb_ecbee42
      vnecb_ecbee45 = vnecb_ecbee44
      vnecb_ecbee47 = vnecb_ecbee46
      vnecb_ecbee49 = vnecb_ecbee48
      vnecb_ecbee51 = vnecb_ecbee50
      vnecb_ecbee53 = vnecb_ecbee52
      vnecb_ecbee55 = vnecb_ecbee54
      vnecb_ecbee57 = vnecb_ecbee56
      vnecb_ecbee59 = vnecb_ecbee58
      vnecb_ecbee61 = vnecb_ecbee60
      vnecb_ecbee63 = vnecb_ecbee62
      return
      end