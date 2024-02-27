*(warnings_cards.h)
      if (process .EQ. 2) then      
      endif
      if (process .EQ. 3) then      
        if ((u_cut(1) .LT. 25.).and.(jgra_sel .GE. 3)) then
          write(6,*) ' '
          write(6,*) '!!!Warning!!!'
          write(6,*) '  ---> Min. cut for u =', u_cut(1), ' GeV^2'
          write(6,*) '  ---> Are you sure?'
          write(6,*) ' '
        endif
        if (W_cut(1) .LT. 5.) then
          write(6,*) ' '
          write(6,*) '!!!Warning!!!'
          write(6,*) '  ---> Min. cut for mass of hadronic system ='
     &                                               ,W_cut(1), ' GeV'
          write(6,*) '  ---> Are you sure?'
          write(6,*) ' '
        endif
        if ((MassELL_cut(1) .LT. 1.).and.(jgra_sel .GE. 3)) then
          write(6,*) ' '
          write(6,*) '!!!Warning!!!'
          write(6,*) '  ---> Min. cut for mass of e,l^+,l^- system ='
     &                                         ,MassELL_cut(1), ' GeV'
          write(6,*) '  ---> Are you sure?'
          write(6,*) ' '
        endif
        if ((MassQLL_cut(1) .LT. 5.).and.(jgra_sel .GE. 3)) then
          write(6,*) ' '
          write(6,*) '!!!Warning!!!'
          write(6,*) '  ---> Min. cut for mass of q,l^+,l^- system ='
     &                                         ,MassQLL_cut(1), ' GeV'
          write(6,*) '  ---> Are you sure?'
          write(6,*) ' '
        endif
      endif
