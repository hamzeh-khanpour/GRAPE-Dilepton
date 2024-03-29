*     File ep_epee/aep_epeekptbl.f : Sat Mar 18 19:44:57 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aep_epeekptbl
      implicit real*8(a-h,o-z)

      include 'inclep_epee1.h'
      include 'inclk.h'
      include 'inclep_epeep.h'
*-----------------------------------------------------------------------
      do 100 i = 1, 4
        pfep_epee2(i) =  + peep_epee1(i)
        pfep_epee4(i) =  + peep_epee2(i)
        pfep_epee6(i) =  + peep_epee1(i) + peep_epee2(i)
        pfep_epee8(i) =  - peep_epee3(i)
        pfep_epee10(i) =  + peep_epee1(i) - peep_epee3(i)
        pfep_epee12(i) =  + peep_epee2(i) - peep_epee3(i)
        pfep_epee14(i) =  + peep_epee1(i) + peep_epee2(i) - peep_epee3(i)
        pfep_epee16(i) =  - peep_epee4(i)
        pfep_epee18(i) =  + peep_epee1(i) - peep_epee4(i)
        pfep_epee20(i) =  + peep_epee2(i) - peep_epee4(i)
        pfep_epee22(i) =  + peep_epee1(i) + peep_epee2(i) - peep_epee4(i)
        pfep_epee24(i) =  - peep_epee3(i) - peep_epee4(i)
        pfep_epee26(i) =  + peep_epee1(i) - peep_epee3(i) - peep_epee4(i)
        pfep_epee28(i) =  + peep_epee2(i) - peep_epee3(i) - peep_epee4(i)
        pfep_epee30(i) =  + peep_epee5(i) + peep_epee6(i)
        pfep_epee32(i) =  - peep_epee5(i)
        pfep_epee34(i) =  + peep_epee1(i) - peep_epee5(i)
        pfep_epee36(i) =  + peep_epee2(i) - peep_epee5(i)
        pfep_epee38(i) =  + peep_epee1(i) + peep_epee2(i) - peep_epee5(i)
        pfep_epee40(i) =  - peep_epee3(i) - peep_epee5(i)
        pfep_epee42(i) =  + peep_epee1(i) - peep_epee3(i) - peep_epee5(i)
        pfep_epee44(i) =  + peep_epee2(i) - peep_epee3(i) - peep_epee5(i)
        pfep_epee46(i) =  + peep_epee4(i) + peep_epee6(i)
        pfep_epee48(i) =  - peep_epee4(i) - peep_epee5(i)
        pfep_epee50(i) =  + peep_epee1(i) - peep_epee4(i) - peep_epee5(i)
        pfep_epee52(i) =  + peep_epee2(i) - peep_epee4(i) - peep_epee5(i)
        pfep_epee54(i) =  + peep_epee3(i) + peep_epee6(i)
        pfep_epee56(i) =  - peep_epee3(i) - peep_epee4(i) - peep_epee5(i)
        pfep_epee58(i) =  - peep_epee2(i) + peep_epee6(i)
        pfep_epee60(i) =  - peep_epee1(i) + peep_epee6(i)
        pfep_epee62(i) =  + peep_epee6(i)

        pfep_epee3(i) = - pfep_epee2(i)
        pfep_epee5(i) = - pfep_epee4(i)
        pfep_epee7(i) = - pfep_epee6(i)
        pfep_epee9(i) = - pfep_epee8(i)
        pfep_epee11(i) = - pfep_epee10(i)
        pfep_epee13(i) = - pfep_epee12(i)
        pfep_epee15(i) = - pfep_epee14(i)
        pfep_epee17(i) = - pfep_epee16(i)
        pfep_epee19(i) = - pfep_epee18(i)
        pfep_epee21(i) = - pfep_epee20(i)
        pfep_epee23(i) = - pfep_epee22(i)
        pfep_epee25(i) = - pfep_epee24(i)
        pfep_epee27(i) = - pfep_epee26(i)
        pfep_epee29(i) = - pfep_epee28(i)
        pfep_epee31(i) = - pfep_epee30(i)
        pfep_epee33(i) = - pfep_epee32(i)
        pfep_epee35(i) = - pfep_epee34(i)
        pfep_epee37(i) = - pfep_epee36(i)
        pfep_epee39(i) = - pfep_epee38(i)
        pfep_epee41(i) = - pfep_epee40(i)
        pfep_epee43(i) = - pfep_epee42(i)
        pfep_epee45(i) = - pfep_epee44(i)
        pfep_epee47(i) = - pfep_epee46(i)
        pfep_epee49(i) = - pfep_epee48(i)
        pfep_epee51(i) = - pfep_epee50(i)
        pfep_epee53(i) = - pfep_epee52(i)
        pfep_epee55(i) = - pfep_epee54(i)
        pfep_epee57(i) = - pfep_epee56(i)
        pfep_epee59(i) = - pfep_epee58(i)
        pfep_epee61(i) = - pfep_epee60(i)
        pfep_epee63(i) = - pfep_epee62(i)
  100 continue

      vnep_epee2 =  + amass2(1)
      vnep_epee4 =  + amass2(2)
      vnep_epee6 =  + 2.0d0*prod(1,2) + amass2(1) + amass2(2)
      vnep_epee8 =  + amass2(3)
      vnep_epee10 =  - 2.0d0*prod(1,3) + amass2(1) + amass2(3)
      vnep_epee12 =  - 2.0d0*prod(2,3) + amass2(2) + amass2(3)
      vnep_epee14 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,3) - 2.0d0*
     &      prod(2,3) + amass2(1) + amass2(2) + amass2(3)
      vnep_epee16 =  + amass2(4)
      vnep_epee18 =  - 2.0d0*prod(1,4) + amass2(1) + amass2(4)
      vnep_epee20 =  - 2.0d0*prod(2,4) + amass2(2) + amass2(4)
      vnep_epee22 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,4) - 2.0d0*
     &      prod(2,4) + amass2(1) + amass2(2) + amass2(4)
      vnep_epee24 =  + 2.0d0*prod(3,4) + amass2(3) + amass2(4)
      vnep_epee26 =  - 2.0d0*prod(1,3) - 2.0d0*prod(1,4) + 2.0d0*
     &      prod(3,4) + amass2(1) + amass2(3) + amass2(4)
      vnep_epee28 =  - 2.0d0*prod(2,3) - 2.0d0*prod(2,4) + 2.0d0*
     &      prod(3,4) + amass2(2) + amass2(3) + amass2(4)
      vnep_epee30 =  + 2.0d0*prod(5,6) + amass2(5) + amass2(6)
      vnep_epee32 =  + amass2(5)
      vnep_epee34 =  - 2.0d0*prod(1,5) + amass2(1) + amass2(5)
      vnep_epee36 =  - 2.0d0*prod(2,5) + amass2(2) + amass2(5)
      vnep_epee38 =  + 2.0d0*prod(1,2) - 2.0d0*prod(1,5) - 2.0d0*
     &      prod(2,5) + amass2(1) + amass2(2) + amass2(5)
      vnep_epee40 =  + 2.0d0*prod(3,5) + amass2(3) + amass2(5)
      vnep_epee42 =  - 2.0d0*prod(1,3) - 2.0d0*prod(1,5) + 2.0d0*
     &      prod(3,5) + amass2(1) + amass2(3) + amass2(5)
      vnep_epee44 =  - 2.0d0*prod(2,3) - 2.0d0*prod(2,5) + 2.0d0*
     &      prod(3,5) + amass2(2) + amass2(3) + amass2(5)
      vnep_epee46 =  + 2.0d0*prod(4,6) + amass2(4) + amass2(6)
      vnep_epee48 =  + 2.0d0*prod(4,5) + amass2(4) + amass2(5)
      vnep_epee50 =  - 2.0d0*prod(1,4) - 2.0d0*prod(1,5) + 2.0d0*
     &      prod(4,5) + amass2(1) + amass2(4) + amass2(5)
      vnep_epee52 =  - 2.0d0*prod(2,4) - 2.0d0*prod(2,5) + 2.0d0*
     &      prod(4,5) + amass2(2) + amass2(4) + amass2(5)
      vnep_epee54 =  + 2.0d0*prod(3,6) + amass2(3) + amass2(6)
      vnep_epee56 =  + 2.0d0*prod(3,4) + 2.0d0*prod(3,5) + 2.0d0*
     &      prod(4,5) + amass2(3) + amass2(4) + amass2(5)
      vnep_epee58 =  - 2.0d0*prod(2,6) + amass2(2) + amass2(6)
      vnep_epee60 =  - 2.0d0*prod(1,6) + amass2(1) + amass2(6)
      vnep_epee62 =  + amass2(6)
      vnep_epee3 = vnep_epee2
      vnep_epee5 = vnep_epee4
      vnep_epee7 = vnep_epee6
      vnep_epee9 = vnep_epee8
      vnep_epee11 = vnep_epee10
      vnep_epee13 = vnep_epee12
      vnep_epee15 = vnep_epee14
      vnep_epee17 = vnep_epee16
      vnep_epee19 = vnep_epee18
      vnep_epee21 = vnep_epee20
      vnep_epee23 = vnep_epee22
      vnep_epee25 = vnep_epee24
      vnep_epee27 = vnep_epee26
      vnep_epee29 = vnep_epee28
      vnep_epee31 = vnep_epee30
      vnep_epee33 = vnep_epee32
      vnep_epee35 = vnep_epee34
      vnep_epee37 = vnep_epee36
      vnep_epee39 = vnep_epee38
      vnep_epee41 = vnep_epee40
      vnep_epee43 = vnep_epee42
      vnep_epee45 = vnep_epee44
      vnep_epee47 = vnep_epee46
      vnep_epee49 = vnep_epee48
      vnep_epee51 = vnep_epee50
      vnep_epee53 = vnep_epee52
      vnep_epee55 = vnep_epee54
      vnep_epee57 = vnep_epee56
      vnep_epee59 = vnep_epee58
      vnep_epee61 = vnep_epee60
      vnep_epee63 = vnep_epee62
      return
      end
