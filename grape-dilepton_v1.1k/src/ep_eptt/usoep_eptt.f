*     File ep_eptt/usoep_eptt.f : Sat Mar 18 19:44:58 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      subroutine usoep_eptt
      implicit real*8(a-h,o-z)

      include 'inclep_eptt1.h'
      include 'inclk.h'
*-----------------------------------------------------------------------
      if(jgraph .lt. 1) then
        write(*,*) '*** No graphs ***'
        return
      endif

      write(*,620)
  620 format('1'/' integrated value of square of each graph',/
     & ' graph  ',8x,'absolute',10x,'relative')

      rn = 1.0d0/(fkcall + nkcall)
      do 100 ig = 1, jgraph
        write(*,*) igraph(ig), ' :',  ansp(ig)*rn, ansp(ig)/ansp(0)
  100 continue
      write(*,*) ' total :', ansp(0)*rn
  600 format(1x,i5,a, e18.8, e18.8)
  610 format(a, e18.8)

      write(*,*)
      write(*,*)
      write(*,*) 'WARNING: Squared Amplitude:'
      write(*,*) '       It should be multiplied by the statistical ',
     &                  'factor 1/N! for N' 
      write(*,*) '       identical final state particle to get the ',
     &                  'total X-section.'

      return
      end
