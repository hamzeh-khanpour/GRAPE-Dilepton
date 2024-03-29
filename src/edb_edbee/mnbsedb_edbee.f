*     File edb_edbee/mnbsedb_edbee.f : Sat Mar 18 19:45:01 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
      subroutine mnbsedb_edbee
************************************************************************
* main program for BASES 5.1
      implicit real*8(a-h,o-z)

      common / knmflg / knmgvs, knmhst, knmsph, knmisr
      external func

*-----------------------------------------------------------------------

      write(*,'(20x,a/)') 'grace 2.1(5)'
      write(*,'(10x,a/)')
     .'(c)Copyright 1990-1998 Minami-Tateya Group (Japan)'
C     write(*,'(1x/)')
*-----------------------------------------------------------------------
*               initialization Kinematics Flags
*     knmgvs  : (=1)sqrt(S) is given in kinit.f
*     knmhst  : (=1)Fill historgram 
*     knmsph  : reserved for soft photon check 
*     knmisr  : reserved for initial state radiation
      knmgvs = 1
      knmhst = 1
      knmsph = 0
      knmisr = 0
************************************************************************
*               initialization of BASES/SPRING 5.1
************************************************************************
*=======================================================================
*          initialization of bases by calling bsinit
*=======================================================================
*         -------------
           call bsinit
*         -------------
*=======================================================================
*      initialization of parameters
*          for kinematics and matrix elements
*      initialization of histograms
*=======================================================================
      lu     = 11

      open(lu,file='bases.result',status='unknown',form='formatted')

*         ------------------
           call userin( lu )
*         ------------------
************************************************************************
*              numerical integration by BASES 5.1
************************************************************************

      call bases( func, estim, sigma, ctime, it1, it2 )

      call bsinfo( lu )

      call bhplot( lu )

      close( lu )

************************************************************************
*             event generation by SPRING 5.1
************************************************************************
*=======================================================================
*     save the probability information to the file
*=======================================================================

      lun = 23
      

      call bswrit( lun )
      close ( lun )

      stop
      end
