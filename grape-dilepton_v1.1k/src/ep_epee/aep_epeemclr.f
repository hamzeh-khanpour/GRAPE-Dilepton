*     File ep_epee/aep_epeemclr.f : Sat Mar 18 19:44:57 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aep_epeemclr
      implicit real*8(a-h,o-z)

      include 'inclep_epee1.h'
      include 'inclk.h'

           do 300 igr = 0, nep_epeegraph
             ansp(igr) = 0.0d0
  300      continue
           do 310 igr = 1, nep_epeegraph
             ancp(igr) = 0.0d0
  310      continue
           fkcall = 0
           nkcall = 0

      return
      end
