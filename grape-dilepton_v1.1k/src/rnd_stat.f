*******************************************************
*      Reading/Dumping random number parameters       *
*                 written by T.Abe                    *
*                on May. 04 in 1999                   *
*******************************************************
*********************************************************************
*  -- LRND : Main switch                                            *
*  -- Parameters are read   from a file named 'rndstat.dat.prev'.   *
*       (If it doesn't exist, default parameters are used.)         *
*  -- Parameters are dumped into a file named 'rndstat.dat'.        *
*********************************************************************
      subroutine RND_Init(LUN, LRND)
      implicit NONE
* ------------ Argument ------------
      integer   LUN
      logical   LRND
* ----------------------------------
* ------------- DRN common --------------
      COMMON/RANDM/RDM(31),RM1,RM2,IA1,IC1,M1,IX1,
     .                             IA2,IC2,M2,IX2,
     .                             IA3,IC3,M3,IX3
      real         RDM    ,RM1,RM2
      integer                      IA1,IC1,M1,IX1,
     &                             IA2,IC2,M2,IX2,
     &                             IA3,IC3,M3,IX3
* ---------------------------------------
* -------- Local variables --------
      integer  i
* ---------------------------------
      if (.not.LRND)  RETURN
      write(6,*) '========== RND_Init started =========='
      open(LUN, file='rndstat.dat.prev', status='OLD', err=1000
     &            ,form='unformatted')
      write(6,*) 'RND_Init>> rndstat.dat.prev was FOUND and opened.'
 100  read(LUN, err=166) (RDM(i),i=1,31)
     &                          ,RM1,RM2,IA1,IC1,M1,IX1,
     &                             IA2,IC2,M2,IX2,
     &                             IA3,IC3,M3,IX3
      GOTO 300
 166  write(6,*) '!!!Error in reading DRN_stuff!!!'
      write(6,*) '  ---> Good-bye!'
      STOP
 300  write(6,*) 'RND_Init>> rndstat.dat.prev was read successfully.'
      close(LUN)
      GOTO 9999
 1000 continue
      write(6,*) 'RND_Init>> rndstat.dat.prev was NOT found.'
      write(6,*) 'RND_Init>> Random numbers initialized with DEFAULTs'
      GOTO 9999
 9999 continue
      write(6,*) '========== RND_Init  ended  =========='
      return
      end
*########################################################################
      subroutine RND_Dump(LUN, LRND)
      implicit NONE
* ------------ Argument ------------
      integer   LUN
      logical   LRND
* ----------------------------------
* ------------- DRN common --------------
      COMMON/RANDM/RDM(31),RM1,RM2,IA1,IC1,M1,IX1,
     .                             IA2,IC2,M2,IX2,
     .                             IA3,IC3,M3,IX3
      real         RDM    ,RM1,RM2
      integer                      IA1,IC1,M1,IX1,
     &                             IA2,IC2,M2,IX2,
     &                             IA3,IC3,M3,IX3
* ---------------------------------------
* -------- Local variables --------
      integer  i
* ---------------------------------
      if (.not.LRND)  RETURN
      write(6,*) '========== RND_Dump started =========='
      open(LUN, file='rndstat.dat', status='unknown', err=1000
     &            ,form='unformatted')
      write(6,*) 'RND_Dump>> rndstat.dat was opened.'
 100  write(LUN, err=166) (RDM(i),i=1,31)
     &                          ,RM1,RM2,IA1,IC1,M1,IX1,
     &                             IA2,IC2,M2,IX2,
     &                             IA3,IC3,M3,IX3
      GOTO 300
 166  write(6,*) '!!!Error in writing DRN_stuff!!!'
      write(6,*) '  ---> Good-bye!'
      STOP
 300  write(6,*) 'RND_Dump>> Dumping was finished successfully.'
      close(LUN)
      GOTO 9999
 1000 continue
      write(6,*) 'RND_Dump>> Opening rndstat.dat was failed.'
      write(6,*) 'RND_Dump>> Random number parameters were not dumped.'
      GOTO 9999
 9999 continue
      write(6,*) '========== RND_Dump  ended  =========='
      return
      end
*########################################################################
