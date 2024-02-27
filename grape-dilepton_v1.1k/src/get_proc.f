*************************************************
*      Getting process ID: /amjprc/jproc        *
*               written by T.Abe                *
*              on Sep. 04 in 1998               *
*************************************************
* Input : process, lpair, qflv
* Output: /amjprc/jproc
      subroutine  Get_Proc
      implicit NONE
      include './inc/graepia.h'
      integer         jproc      
       COMMON /amjprc/jproc
* -------- Local variables --------
* ---------------------------------
      jproc = 0
      if ((lpair.LT.1).or.(lpair.GT.3)) then
        write(6,*) '!!!Error!!! : unknown lepton flavor(=',lpair,')'
        write(6,*) ' ---> Good-bye!'
        STOP
      endif
      if      (process .EQ. 1) then
         if (lpair.EQ.1)  jproc = 1
         if (lpair.EQ.2)  jproc = 2
         if (lpair.EQ.3)  jproc = 3
         qflv = 0
         write(6,*) '---> Elastic process'
       elseif (process .EQ. 2) then
         if (lpair.EQ.1)  jproc = 1 + num_prc_ela
         if (lpair.EQ.2)  jproc = 2 + num_prc_ela
         if (lpair.EQ.3)  jproc = 3 + num_prc_ela
         qflv = 0
         write(6,*) '---> Quasi-elastic process'
       elseif (process .EQ. 3) then
         jproc = (qflv-1)*3 + lpair + num_prc_ela_qela
         qflv = int( (jproc-1-num_prc_ela_qela)/3 ) + 1   
         write(6,*)   '---> DIS process'
        if      (qflv.EQ. 1) then
           write(6,*) '       (Scattering of e and u quark) '
         elseif (qflv.EQ. 2) then
           write(6,*) '       (Scattering of e and u-bar quark) '
         elseif (qflv.EQ. 3) then
           write(6,*) '       (Scattering of e and d quark) '
         elseif (qflv.EQ. 4) then
           write(6,*) '       (Scattering of e and d-bar quark) '
         elseif (qflv.LE. 5) then
           write(6,*) '       (Scattering of e and s quark) '
         elseif (qflv.LE. 6) then
           write(6,*) '       (Scattering of e and s-bar quark) '
         elseif (qflv.LE. 7) then
           write(6,*) '       (Scattering of e and c quark) '
         elseif (qflv.LE. 8) then
           write(6,*) '       (Scattering of e and c-bar quark) '
         elseif (qflv.LE. 9) then
           write(6,*) '       (Scattering of e and b quark) '
         elseif (qflv.LE.10) then
           write(6,*) '       (Scattering of e and b-bar quark) '
         elseif (qflv.LE.11) then
           write(6,*) '       (Scattering of e and t quark) '
         elseif (qflv.LE.12) then
           write(6,*) '       (Scattering of e and t-bar quark) '
         else
           write(6,*) '!!!Error!!! : unknown quark flavor(=',qflv,')'
           write(6,*) ' ---> Good-bye!'
           STOP
        endif
       elseif (process .EQ. 4) then
         jproc = (qflv-1)*3 + lpair + num_prc_4f
         qflv = int( (jproc-1-num_prc_4f)/3 ) + 1   
         write(6,*)   '---> Resolved process'
        if      (qflv.EQ. 1) then
           write(6,*) '       (Scattering of u and u-bar quark) '
         elseif (qflv.EQ. 2) then
           write(6,*) '       (Scattering of u-bar and u quark) '
         elseif (qflv.EQ. 3) then
           write(6,*) '       (Scattering of d and d-bar quark) '
         elseif (qflv.EQ. 4) then
           write(6,*) '       (Scattering of d-bar and d quark) '
         elseif (qflv.LE. 5) then
           write(6,*) '       (Scattering of s and s-bar quark) '
         elseif (qflv.LE. 6) then
           write(6,*) '       (Scattering of s-bar and s quark) '
         elseif (qflv.LE. 7) then
           write(6,*) '       (Scattering of c and c-bar quark) '
         elseif (qflv.LE. 8) then
           write(6,*) '       (Scattering of c-bar and c quark) '
         elseif (qflv.LE. 9) then
           write(6,*) '       (Scattering of b and b-bar quark) '
         elseif (qflv.LE.10) then
           write(6,*) '       (Scattering of b-bar and b quark) '
         elseif (qflv.LE.11) then
           write(6,*) '       (Scattering of t and t-bar quark) '
         elseif (qflv.LE.12) then
           write(6,*) '       (Scattering of t-bar and t quark) '
         else
           write(6,*) '!!!Error!!! : unknown quark flavor(=',qflv,')'
           write(6,*) ' ---> Good-bye!'
           STOP
        endif
       elseif (process .EQ. 5) then
         jproc = (qflv-1)*3 + lpair + num_prc_4f_qq
         qflv = int( (jproc-1-num_prc_4f_qq)/3 ) + 1   
         write(6,*)   '---> Direct process'
        if      (qflv.EQ. 1) then
           write(6,*) '       (Scattering of u and photon quark) '
         elseif (qflv.EQ. 2) then
           write(6,*) '       (Scattering of u-bar and photon quark) '
         elseif (qflv.EQ. 3) then
           write(6,*) '       (Scattering of d and photon quark) '
         elseif (qflv.EQ. 4) then
           write(6,*) '       (Scattering of d-bar and photon quark) '
         elseif (qflv.LE. 5) then
           write(6,*) '       (Scattering of s and photon quark) '
         elseif (qflv.LE. 6) then
           write(6,*) '       (Scattering of s-bar and photon quark) '
         elseif (qflv.LE. 7) then
           write(6,*) '       (Scattering of c and photon quark) '
         elseif (qflv.LE. 8) then
           write(6,*) '       (Scattering of c-bar and photon quark) '
         elseif (qflv.LE. 9) then
           write(6,*) '       (Scattering of b and photon quark) '
         elseif (qflv.LE.10) then
           write(6,*) '       (Scattering of b-bar and photon quark) '
         elseif (qflv.LE.11) then
           write(6,*) '       (Scattering of t and photon quark) '
         elseif (qflv.LE.12) then
           write(6,*) '       (Scattering of t-bar and photon quark) '
         else
           write(6,*) '!!!Error!!! : unknown quark flavor(=',qflv,')'
           write(6,*) ' ---> Good-bye!'
           STOP
        endif
       else
         write(6,*) '!!!Error!!! : unknown process(=',process,')'
         write(6,*) ' ---> Good-bye!'
         STOP
      endif
      if      (lpair.EQ.1) then
         write(6,*) '---> Electron-pair production '
       elseif (lpair.EQ.2) then
         write(6,*) '---> Muon-pair production '
       elseif (lpair.EQ.3) then
         write(6,*) '---> Tau-pair production '
       else
           write(6,*) '!!!Error!!! : unknown lepton flavor(=',lpair,')'
           write(6,*) ' ---> Good-bye!'
           STOP
      endif
      if (jproc.EQ.0) then
        write(6,*) '!!!Error!!! : jproc was not set correctly.'
        write(6,*) ' ---> Good-bye!'
        STOP
      endif
      return
      end
