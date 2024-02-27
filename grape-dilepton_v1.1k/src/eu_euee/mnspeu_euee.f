*     File eu_euee/mnspeu_euee.f : Sat Mar 18 19:44:59 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
      subroutine mnspeu_euee
************************************************************************
* main program for SPRING 5.1
      implicit real*8(a-h,o-z)

      include 'inclk.h'
      common / knmflg / knmgvs, knmhst, knmsph, knmisr
      external func
      parameter (nextn  =6)
      common /sp4vec/ vec(4,nextn)

*-----------------------------------------------------------------------

      write(*,'(20x,a/)') 'grace 2.1(5)'
      write(*,'(10x,a/)')
     .'(c)Copyright 1990-1998 Minami-Tateya Group (Japan)'
      write(*,'(5x,a)') 'Initial particles:'
      write(*,'(25x,a)') 'u'
      write(*,'(25x,a)') 'positron'
      write(*,'(5x,a)') 'Final particles:'
      write(*,'(25x,a)') 'u'
      write(*,'(25x,a)') 'positron'
      write(*,'(25x,a)') 'electron'
      write(*,'(25x,a)') 'positron'
      write(*,'(1x/)')
*-----------------------------------------------------------------------
*               initialization Kinematics Flags
*     knmgvs  : (=1)sqrt(S) is given in kinit.f
*     knmhst  : (=1)Fill historgram for BASES
*     knmsph  : reserved for soft photon check 
*     knmisr  : reserved for initial state radiation
      knmgvs = 1
      knmhst = 0
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
*          read the probability information from the file
*=======================================================================
           lun = 23
           open(lun,file='bases.data',status='old',form='unformatted')
*         -------------
           call bsread( lun )
*         -------------
           close( lun )

*=======================================================================
*      initialization of parameters
*          for kinematics and matrix elements
*      initialization of histograms
*=======================================================================
      lu     = 12

      open(lu,file='spring.result',status='unknown',form='formatted')

*         ------------------
           call userin( lu )
*         ------------------

*=======================================================================
*     initialization of additional histograms for spring
*=======================================================================

*=======================================================================
*     event generation
*=======================================================================

      mxtry  = 50

      mxevnt = 10000

*     write(6,*)' number of events ? '

*     read(5,*) mxevnt

      do 100 nevnt = 1, mxevnt

         call spring( func, mxtry )


*-----------------------------------------------------------------------
*        compute the four vectors of generated event
*          from the kinematical variables
*-----------------------------------------------------------------------

*        do 90 k = 1 , nextn
*           write(6,*) (vec(j,k),j=1,4)
*  90    continue
  100 continue


      call spinfo( lu )

      call shplot( lu )

      close( lu )

      stop
      end

