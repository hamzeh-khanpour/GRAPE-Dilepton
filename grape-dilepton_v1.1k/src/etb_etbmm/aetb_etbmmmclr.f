*     File etb_etbmm/aetb_etbmmmclr.f : Sat Mar 18 19:45:07 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine aetb_etbmmmclr
      implicit real*8(a-h,o-z)

      include 'incletb_etbmm1.h'
      include 'inclk.h'

           do 300 igr = 0, netb_etbmmgraph
             ansp(igr) = 0.0d0
  300      continue
           do 310 igr = 1, netb_etbmmgraph
             ancp(igr) = 0.0d0
  310      continue
           fkcall = 0
           nkcall = 0

      return
      end
