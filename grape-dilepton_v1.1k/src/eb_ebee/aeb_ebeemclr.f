*     File eb_ebee/aeb_ebeemclr.f : Sat Mar 18 19:45:05 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aeb_ebeemclr
      implicit real*8(a-h,o-z)

      include 'incleb_ebee1.h'
      include 'inclk.h'

           do 300 igr = 0, neb_ebeegraph
             ansp(igr) = 0.0d0
  300      continue
           do 310 igr = 1, neb_ebeegraph
             ancp(igr) = 0.0d0
  310      continue
           fkcall = 0
           nkcall = 0

      return
      end