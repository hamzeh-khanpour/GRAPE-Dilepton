*******************************************************
*      Getting flags for Feynman Graph Selection      *
*                  written by T.Abe                   *
*                 on Aug. 09 in 1998                  *
*******************************************************
* Input  : jgra_sel(<>0), process, lpair, qflv *
* Output : jgra_flag (--> used in gfinit)      *
************************************************
* Input  : jgra_sel(=0), jgra_flag             *
* Output : jgra_flag (--> used in gfinit)      *
************************************************
      subroutine  Get_Graph_flag
      implicit NONE
      include './inc/graepia.h'
      integer         jproc      
       COMMON /amjprc/jproc
      integer  i, sum, N
      if (jgra_sel .LT. 0) then
         jgra_sel = - jgra_sel
         write(6,*)
     &         '!!!!!! USER Specification for Graph Selection !!!!!!'
         GOTO 222
       else
         do i = 1, num_gra_flg
           jgra_flag(i) = 0
         enddo
      endif
      if ((jgra_sel.LE.0).or.(jgra_sel.GT.14)) then
         write(6,*) '!!!Error in Get_Graph_flag!!!'
         write(6,*) '  ---> Unknown graph selection'
         write(6,*) '        ( jgra_sel =',jgra_sel,' )'
         write(6,*) '  ---> Good-bye!'
         STOP
      endif
* ---------- For ep --> ep(X)ee process ----------
      if (      (process.EQ.1 .OR. process.EQ.2)
     +    .and. (lpair.EQ.1) ) then
        if      (jgra_sel .EQ. 0) then
           continue
         elseif (jgra_sel .EQ. 1) then
           call Add_Graph(jgra_flag, num_gra_flg,   11)   ! BH(direct)
           call Add_Graph(jgra_flag, num_gra_flg,   15)   ! BH(direct)
         elseif (jgra_sel .EQ. 2) then
           call Add_Graph(jgra_flag, num_gra_flg,   11)   ! BH(direct)
           call Add_Graph(jgra_flag, num_gra_flg,   15)   ! BH(direct)
           call Add_Graph(jgra_flag, num_gra_flg,    7)   ! BH(exchanged)
           call Add_Graph(jgra_flag, num_gra_flg,    9)   ! BH(exchanged)
         elseif (jgra_sel .EQ. 3) then
           call Add_Graph(jgra_flag, num_gra_flg,   11)   ! BH(direct)
           call Add_Graph(jgra_flag, num_gra_flg,   15)   ! BH(direct)
           call Add_Graph(jgra_flag, num_gra_flg,    7)   ! BH(exchanged)
           call Add_Graph(jgra_flag, num_gra_flg,    9)   ! BH(exchanged)
           call Add_Graph(jgra_flag, num_gra_flg,    1)   ! CO(direct)
           call Add_Graph(jgra_flag, num_gra_flg,    5)   ! CO(direct)
           call Add_Graph(jgra_flag, num_gra_flg,    3)   ! CO(exchanged)
           call Add_Graph(jgra_flag, num_gra_flg,   13)   ! CO(exchanged)
         elseif ((jgra_sel .EQ. 4).or.(jgra_sel .EQ. 5)) then
           do i=1,16
             call Add_Graph(jgra_flag, num_gra_flg,  i)
           enddo
         elseif (jgra_sel .EQ. 13) then
           call Add_Graph(jgra_flag, num_gra_flg,    1)   ! CO(direct)
           call Add_Graph(jgra_flag, num_gra_flg,    5)   ! CO(direct)
           call Add_Graph(jgra_flag, num_gra_flg,    3)   ! CO(exchanged)
           call Add_Graph(jgra_flag, num_gra_flg,   13)   ! CO(exchanged)
         elseif (jgra_sel .EQ. 14) then
           call Add_Graph(jgra_flag, num_gra_flg,    2)   ! Z0(dir)
           call Add_Graph(jgra_flag, num_gra_flg,    6)   ! Z0(dir)
           call Add_Graph(jgra_flag, num_gra_flg,    4)   ! Z0(exch)
           call Add_Graph(jgra_flag, num_gra_flg,   14)   ! Z0(exch)
         else
           GOTO 9100   ! Unknown jgra_sel
        endif
        GOTO 222
      endif
* ---------- For ep --> ep(X)ll process ----------
      if (      (process.EQ.1 .OR. process.EQ.2)
     +    .and. (lpair.EQ.2 .OR. lpair.EQ.3) ) then
        if      (jgra_sel .EQ. 0) then
           continue
         elseif ((jgra_sel .EQ. 1).or.(jgra_sel .EQ. 2)) then
           call Add_Graph(jgra_flag, num_gra_flg,    5)   ! BH
           call Add_Graph(jgra_flag, num_gra_flg,    7)   ! BH
         elseif (jgra_sel .EQ. 3) then
           call Add_Graph(jgra_flag, num_gra_flg,    5)   ! BH
           call Add_Graph(jgra_flag, num_gra_flg,    7)   ! BH
           call Add_Graph(jgra_flag, num_gra_flg,    1)   ! CO
           call Add_Graph(jgra_flag, num_gra_flg,    3)   ! CO
         elseif ((jgra_sel .EQ. 4).or.(jgra_sel .EQ. 5)) then
           do i=1,8
             call Add_Graph(jgra_flag, num_gra_flg,  i)
           enddo
         elseif (jgra_sel .EQ. 13) then
           call Add_Graph(jgra_flag, num_gra_flg,    1)   ! CO
           call Add_Graph(jgra_flag, num_gra_flg,    3)   ! CO
         elseif (jgra_sel .EQ. 14) then
           call Add_Graph(jgra_flag, num_gra_flg,    2)   ! Z0
           call Add_Graph(jgra_flag, num_gra_flg,    4)   ! Z0
         else
           GOTO 9100   ! Unknown jgra_sel
        endif
        GOTO 222
      endif
* ---------- For eq --> eqee process (light quarks) ----------
      if ( (process.EQ.3) .and. (lpair.EQ.1) .and. (qflv.LE.6) ) then
        if      (jgra_sel .EQ. 1) then
           call Add_Graph(jgra_flag, num_gra_flg,   21)   ! BH(direct)
           call Add_Graph(jgra_flag, num_gra_flg,   29)   ! BH(direct)
         elseif (jgra_sel .EQ. 2) then
           call Add_Graph(jgra_flag, num_gra_flg,   21)   ! BH(direct)
           call Add_Graph(jgra_flag, num_gra_flg,   29)   ! BH(direct)
           call Add_Graph(jgra_flag, num_gra_flg,   13)   ! BH(exchanged)
           call Add_Graph(jgra_flag, num_gra_flg,   17)   ! BH(exchanged)
         elseif (jgra_sel .EQ. 3) then
           call Add_Graph(jgra_flag, num_gra_flg,   21)   ! BH(direct)
           call Add_Graph(jgra_flag, num_gra_flg,   29)   ! BH(direct)
           call Add_Graph(jgra_flag, num_gra_flg,   13)   ! BH(exchanged)
           call Add_Graph(jgra_flag, num_gra_flg,   17)   ! BH(exchanged)
           call Add_Graph(jgra_flag, num_gra_flg,    1)   ! CO(direct)
           call Add_Graph(jgra_flag, num_gra_flg,    9)   ! CO(direct)
           call Add_Graph(jgra_flag, num_gra_flg,   33)   ! CO(direct)
           call Add_Graph(jgra_flag, num_gra_flg,   45)   ! CO(direct)
           call Add_Graph(jgra_flag, num_gra_flg,    5)   ! CO(exchanged)
           call Add_Graph(jgra_flag, num_gra_flg,   25)   ! CO(exchanged)
           call Add_Graph(jgra_flag, num_gra_flg,   37)   ! CO(exchanged)
           call Add_Graph(jgra_flag, num_gra_flg,   41)   ! CO(exchanged)
         elseif ((jgra_sel .EQ. 4).or.(jgra_sel .EQ. 5)) then
           do i=1,48
             call Add_Graph(jgra_flag, num_gra_flg, i)
           enddo
         elseif (jgra_sel .EQ. 13) then
           call Add_Graph(jgra_flag, num_gra_flg,    1)   ! CO(direct)
           call Add_Graph(jgra_flag, num_gra_flg,    9)   ! CO(direct)
           call Add_Graph(jgra_flag, num_gra_flg,   33)   ! CO(direct)
           call Add_Graph(jgra_flag, num_gra_flg,   45)   ! CO(direct)
           call Add_Graph(jgra_flag, num_gra_flg,    5)   ! CO(exchanged)
           call Add_Graph(jgra_flag, num_gra_flg,   25)   ! CO(exchanged)
           call Add_Graph(jgra_flag, num_gra_flg,   37)   ! CO(exchanged)
           call Add_Graph(jgra_flag, num_gra_flg,   41)   ! CO(exchanged)
         elseif (jgra_sel .EQ. 14) then
           call Add_Graph(jgra_flag, num_gra_flg,    2)   ! Z0(direct)
           call Add_Graph(jgra_flag, num_gra_flg,   10)   ! Z0(direct)
           call Add_Graph(jgra_flag, num_gra_flg,   34)   ! Z0(direct)
           call Add_Graph(jgra_flag, num_gra_flg,   47)   ! Z0(direct)
           call Add_Graph(jgra_flag, num_gra_flg,    6)   ! Z0(exchanged)
           call Add_Graph(jgra_flag, num_gra_flg,   26)   ! Z0(exchanged)
           call Add_Graph(jgra_flag, num_gra_flg,   38)   ! Z0(exchanged)
           call Add_Graph(jgra_flag, num_gra_flg,   43)   ! Z0(exchanged)
         else
           GOTO 9100   ! Unknown jgra_sel
        endif
        GOTO 222
      endif
* ---------- For eq --> eqll process (light quarks) ----------
      if (   (process.EQ.3)
     +     .and. (lpair.EQ.2 .OR. lpair.EQ.3)
     +     .and. (qflv.LE.6) ) then
        if     ((jgra_sel .EQ. 1).or.(jgra_sel .EQ. 2)) then
           call Add_Graph(jgra_flag, num_gra_flg,   9)   ! BH
           call Add_Graph(jgra_flag, num_gra_flg,  13)   ! BH
         elseif (jgra_sel .EQ. 3) then
           call Add_Graph(jgra_flag, num_gra_flg,   9)   ! BH
           call Add_Graph(jgra_flag, num_gra_flg,  13)   ! BH
           call Add_Graph(jgra_flag, num_gra_flg,   1)   ! CO(e-line)
           call Add_Graph(jgra_flag, num_gra_flg,   5)   ! CO(e-line)
           call Add_Graph(jgra_flag, num_gra_flg,  18)   ! CO(q-line)
           call Add_Graph(jgra_flag, num_gra_flg,  22)   ! CO(q-line)
         elseif (jgra_sel .EQ. 4) then
           do i=1,25
             call Add_Graph(jgra_flag, num_gra_flg, i)
           enddo
           call Sub_Graph(jgra_flag, num_gra_flg,  17)   
         elseif (jgra_sel .EQ. 5) then
           do i=1,25
             call Add_Graph(jgra_flag, num_gra_flg, i)
           enddo
         elseif (jgra_sel .EQ. 13) then
           call Add_Graph(jgra_flag, num_gra_flg,   1)   ! CO(e-line)
           call Add_Graph(jgra_flag, num_gra_flg,   5)   ! CO(e-line)
           call Add_Graph(jgra_flag, num_gra_flg,  18)   ! CO(q-line)
           call Add_Graph(jgra_flag, num_gra_flg,  22)   ! CO(q-line)
         elseif (jgra_sel .EQ. 14) then
           call Add_Graph(jgra_flag, num_gra_flg,   2)   ! Z0(e-line)
           call Add_Graph(jgra_flag, num_gra_flg,   6)   ! Z0(e-line)
           call Add_Graph(jgra_flag, num_gra_flg,  19)   ! Z0(q-line)
           call Add_Graph(jgra_flag, num_gra_flg,  24)   ! Z0(q-line)
         else
           GOTO 9100   ! Unknown jgra_sel
        endif
        GOTO 222
      endif
* ------------------------------------------------------------
* ---------- For eq --> eqee process (heavy quarks) ----------
      if ( (process.EQ.3) .and. (lpair.EQ.1) .and. (qflv.GT.6) ) then
        if      (jgra_sel .EQ. 1) then
           call Add_Graph(jgra_flag, num_gra_flg,  21)   ! BH(direct)
           call Add_Graph(jgra_flag, num_gra_flg,  29)   ! BH(direct)
         elseif (jgra_sel .EQ. 2) then
           call Add_Graph(jgra_flag, num_gra_flg,  21)   ! BH(direct)
           call Add_Graph(jgra_flag, num_gra_flg,  29)   ! BH(direct)
           call Add_Graph(jgra_flag, num_gra_flg,  13)   ! BH(exchanged)
           call Add_Graph(jgra_flag, num_gra_flg,  17)   ! BH(exchanged)
         elseif (jgra_sel .EQ. 3) then
           call Add_Graph(jgra_flag, num_gra_flg,  21)   ! BH(direct)
           call Add_Graph(jgra_flag, num_gra_flg,  29)   ! BH(direct)
           call Add_Graph(jgra_flag, num_gra_flg,  13)   ! BH(exchanged)
           call Add_Graph(jgra_flag, num_gra_flg,  17)   ! BH(exchanged)
           call Add_Graph(jgra_flag, num_gra_flg,   1)   ! CO(direct)
           call Add_Graph(jgra_flag, num_gra_flg,   9)   ! CO(direct)
           call Add_Graph(jgra_flag, num_gra_flg,  35)   ! CO(direct)
           call Add_Graph(jgra_flag, num_gra_flg,  47)   ! CO(direct)
           call Add_Graph(jgra_flag, num_gra_flg,   5)   ! CO(exchanged)
           call Add_Graph(jgra_flag, num_gra_flg,  25)   ! CO(exchanged)
           call Add_Graph(jgra_flag, num_gra_flg,  39)   ! CO(exchanged)
           call Add_Graph(jgra_flag, num_gra_flg,  43)   ! CO(exchanged)
         elseif (jgra_sel .EQ. 4) then
           do i=1,50
             call Add_Graph(jgra_flag, num_gra_flg, i)
           enddo
           call Sub_Graph(jgra_flag, num_gra_flg,  33)   
           call Sub_Graph(jgra_flag, num_gra_flg,  34)   
         elseif (jgra_sel .EQ. 13) then
           call Add_Graph(jgra_flag, num_gra_flg,   1)   ! CO(direct)
           call Add_Graph(jgra_flag, num_gra_flg,   9)   ! CO(direct)
           call Add_Graph(jgra_flag, num_gra_flg,  35)   ! CO(direct)
           call Add_Graph(jgra_flag, num_gra_flg,  47)   ! CO(direct)
           call Add_Graph(jgra_flag, num_gra_flg,   5)   ! CO(exchanged)
           call Add_Graph(jgra_flag, num_gra_flg,  25)   ! CO(exchanged)
           call Add_Graph(jgra_flag, num_gra_flg,  39)   ! CO(exchanged)
           call Add_Graph(jgra_flag, num_gra_flg,  43)   ! CO(exchanged)
         elseif (jgra_sel .EQ. 14) then
           call Add_Graph(jgra_flag, num_gra_flg,   2)   ! Z0(direct)
           call Add_Graph(jgra_flag, num_gra_flg,  10)   ! Z0(direct)
           call Add_Graph(jgra_flag, num_gra_flg,  36)   ! Z0(direct)
           call Add_Graph(jgra_flag, num_gra_flg,  49)   ! Z0(direct)
           call Add_Graph(jgra_flag, num_gra_flg,   6)   ! Z0(exchanged)
           call Add_Graph(jgra_flag, num_gra_flg,  26)   ! Z0(exchanged)
           call Add_Graph(jgra_flag, num_gra_flg,  40)   ! Z0(exchanged)
           call Add_Graph(jgra_flag, num_gra_flg,  45)   ! Z0(exchanged)
         else
           GOTO 9100   ! Unknown jgra_sel
        endif
        GOTO 222
      endif
* ---------- For eq --> eqll process (heavy quarks) ----------
      if (   (process.EQ.3)
     +     .and. (lpair.EQ.2 .OR. lpair.EQ.3)
     +     .and. (qflv.GT.6) ) then
        if     ((jgra_sel .EQ. 1).or.(jgra_sel .EQ. 2)) then
           call Add_Graph(jgra_flag, num_gra_flg,   9)   ! BH
           call Add_Graph(jgra_flag, num_gra_flg,  15)   ! BH
         elseif (jgra_sel .EQ. 3) then
           call Add_Graph(jgra_flag, num_gra_flg,   9)   ! BH
           call Add_Graph(jgra_flag, num_gra_flg,  15)   ! BH
           call Add_Graph(jgra_flag, num_gra_flg,   1)   ! CO(e-line)
           call Add_Graph(jgra_flag, num_gra_flg,   5)   ! CO(e-line)
           call Add_Graph(jgra_flag, num_gra_flg,  24)   ! CO(q-line)
           call Add_Graph(jgra_flag, num_gra_flg,  30)   ! CO(q-line)
         elseif (jgra_sel .EQ. 4) then
           do i=1,35
             call Add_Graph(jgra_flag, num_gra_flg, i)
           enddo
           call Sub_Graph(jgra_flag, num_gra_flg,  13)   
           call Sub_Graph(jgra_flag, num_gra_flg,  14)   
           call Sub_Graph(jgra_flag, num_gra_flg,  19)   
           call Sub_Graph(jgra_flag, num_gra_flg,  20)   
           call Sub_Graph(jgra_flag, num_gra_flg,  21)   
           call Sub_Graph(jgra_flag, num_gra_flg,  22)   
           call Sub_Graph(jgra_flag, num_gra_flg,  23)   
           call Sub_Graph(jgra_flag, num_gra_flg,  26)   
           call Sub_Graph(jgra_flag, num_gra_flg,  29)   
           call Sub_Graph(jgra_flag, num_gra_flg,  34)   
           call Sub_Graph(jgra_flag, num_gra_flg,  35)   
         elseif (jgra_sel .EQ. 13) then
           call Add_Graph(jgra_flag, num_gra_flg,   1)   ! CO(e-line)
           call Add_Graph(jgra_flag, num_gra_flg,   5)   ! CO(e-line)
           call Add_Graph(jgra_flag, num_gra_flg,  24)   ! CO(q-line)
           call Add_Graph(jgra_flag, num_gra_flg,  30)   ! CO(q-line)
         elseif (jgra_sel .EQ. 14) then
           call Add_Graph(jgra_flag, num_gra_flg,   2)   ! Z0(e-line)
           call Add_Graph(jgra_flag, num_gra_flg,   6)   ! Z0(e-line)
           call Add_Graph(jgra_flag, num_gra_flg,  25)   ! Z0(q-line)
           call Add_Graph(jgra_flag, num_gra_flg,  32)   ! Z0(q-line)
         else
           GOTO 9100   ! Unknown jgra_sel
        endif
        GOTO 222
      endif
      write(6,*) '!!!Error in Get_Graph_flag!!!'
      write(6,*) '  ---> jgra_flag =', (jgra_flag(i),i=1,num_gra_flg)
      write(6,*) '  ---> jgra_flag is not specified.'
      write(6,*) '  ---> Good-bye!'
      STOP
 222  continue
      if (print_flag .EQ. 1) then
        write(6,*) ' '
        write(6,1100)     ' ---> jgra_sel  =', jgra_sel
        N = min(num_gra_flg, 5)
        write(6,1500)     ' ---> jgra_flag =', (jgra_flag(i),i=1,N)
        if (num_gra_flg .GT. 5) then
           N = min(num_gra_flg, 10)
           write(6,1500)  '                 ', (jgra_flag(i),i=6,N)
        endif
        if (num_gra_flg .GT.10) then
           N = min(num_gra_flg, 15)
           write(6,1500)  '                 ', (jgra_flag(i),i=11,N)
        endif
        if (num_gra_flg .GT.15) then
           N = min(num_gra_flg, 20)
           write(6,1500)  '                 ', (jgra_flag(i),i=16,N)
        endif
        if (num_gra_flg .GT.20) then
           N = num_gra_flg   ! min(num_gra_flg, 25)
           write(6,1500)  '                 ', (jgra_flag(i),i=21,N)
        endif
        write(6,*) ' '
C        write(6,*) ' '
      endif
1100  format(A17,  (I10,1X))
1500  format(A17, 5(I10,1X))
*>>> Checking whether jgra_sel is set correctly or not.
      sum  = 0
      do i=1,num_gra_flg
        sum  = sum  + jgra_flag(i)
        if (jgra_flag(i) .LT. 0 ) then
          write(6,*) '!!!Error in Get_Graph_flag!!!'
          write(6,*) '  ---> jgra_flag is not set correctly.'
          write(6,*) '  ---> Good-bye!'
          STOP
        endif
      enddo
      if (sum .EQ. 0) then
        write(6,*) '!!!Error in Get_Graph_flag!!!'
        write(6,*) '  ---> jgra_flag is not set correctly.'
        write(6,*) '  ---> Good-bye!'
        STOP
      endif
      RETURN
 9100 continue
      write(6,*) '!!!Error in Get_Graph_Flag!!!'
      write(6,*) '  ---> Unknown jgra_sel(=',jgra_sel,') for EW'
      write(6,*) '  ---> Good-bye!'
      STOP
 9300 continue
 9400 continue
 9500 continue
      return
      end
