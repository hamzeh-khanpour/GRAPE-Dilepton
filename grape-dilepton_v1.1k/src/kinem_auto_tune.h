*(Modified by T.Abe on May 23, 2002)

*(on isym34)
      if (isym_34 .EQ. 1122) then
        if (  (process .LE. 2)
     +    .OR.((jgra_sel .GT. 0).and.(jgra_sel .LE. 2))
     +     ) then
          isym_34 = 0
          ii34    = 3
        else
          isym_34 = 1
        endif
      endif
*(on u-pol cut)
      if (u_cut(1) .LT. 0.) then
        if ((process .EQ. 3).and.(jgra_sel .GE. 3)) then
           u_cut(1) = 25.
           u_cut(2) = 1.E20
         else
           u_cut(1) =  0.
           u_cut(2) = 1.E20
        endif
      endif
*(On eps_p)
      if (Neps_p.LT.0) then
        if (process .LE. 2) then
          Neps_p = 12
        else
          Neps_p = 11
        endif
      endif
*(On eps_p)
      if (NN_ISR.LT.0) then
        if (process .LE. 2) then
          NN_ISR = 56
        else
          NN_ISR = 55
        endif
      endif
*(on reso56)
      if (Ireso56 .EQ. 1122) then
        if     (jgra_sel .LE. 3) then
          Ireso56  = -1
        elseif (jgra_sel .EQ. 4) then
          Ireso56  =  3
          thresh56 = 0.4
        elseif (jgra_sel .EQ. 14) then
          Ireso56  =  3
          thresh56 = 0.0
        else
          Ireso56  = -1
        endif
      endif
*(Max. of Q2p)
      if (.true.
     +     .AND. (process .EQ. 2)       ! Quasi-elastic process
     +     .AND. (W_cut(1) .LT. 2.0)    ! Using Brasse et. al param
     +     .AND. (Q2p_cut(2) .GT. 100)  ! Too large
     +   ) then
C       Q2p_cut(2) = 10000
        write(6,*) ' '
        write(6,*) '!!! Warning !!!'
        write(6,*) '  ---> (Max. cut for Q2p) =', Q2p_cut(2)
        write(6,*) '       is too large for the Brasse et. al param.'
        write(6,*) '  ---> It is set to be 100 GeV^2.'
        Q2p_cut(2) = 100
      endif
