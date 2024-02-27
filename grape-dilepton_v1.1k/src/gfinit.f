*     File gfinit.f : Sat Mar 18 19:45:07 2000
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
************************************************************************
      subroutine gfinit
      implicit real*8(a-h,o-z)
      common /amjprc/jproc
      include './inc/graepia.h'
      if (print_flag .EQ. 1) then
        write(6,*) ' '
        write(6,*) '>>> Graph selection'
        write(6,*) 'jselg ='
      endif
        if(jproc .eq. 1) then
            call gep_epeefinit
        endif
        if(jproc .eq. 2) then
            call gep_epmmfinit
        endif
        if(jproc .eq. 3) then
            call gep_epttfinit
        endif
        if(jproc .eq. 4) then
            call gep_eXeefinit
        endif
        if(jproc .eq. 5) then
            call gep_eXmmfinit
        endif
        if(jproc .eq. 6) then
            call gep_eXttfinit
        endif
        if(jproc .eq. 7) then
            call geu_eueefinit
        endif
        if(jproc .eq. 8) then
            call geu_eummfinit
        endif
        if(jproc .eq. 9) then
            call geu_euttfinit
        endif
        if(jproc .eq. 10) then
            call geub_eubeefinit
        endif
        if(jproc .eq. 11) then
            call geub_eubmmfinit
        endif
        if(jproc .eq. 12) then
            call geub_eubttfinit
        endif
        if(jproc .eq. 13) then
            call ged_edeefinit
        endif
        if(jproc .eq. 14) then
            call ged_edmmfinit
        endif
        if(jproc .eq. 15) then
            call ged_edttfinit
        endif
        if(jproc .eq. 16) then
            call gedb_edbeefinit
        endif
        if(jproc .eq. 17) then
            call gedb_edbmmfinit
        endif
        if(jproc .eq. 18) then
            call gedb_edbttfinit
        endif
        if(jproc .eq. 19) then
            call ges_eseefinit
        endif
        if(jproc .eq. 20) then
            call ges_esmmfinit
        endif
        if(jproc .eq. 21) then
            call ges_esttfinit
        endif
        if(jproc .eq. 22) then
            call gesb_esbeefinit
        endif
        if(jproc .eq. 23) then
            call gesb_esbmmfinit
        endif
        if(jproc .eq. 24) then
            call gesb_esbttfinit
        endif
        if(jproc .eq. 25) then
            call gec_eceefinit
        endif
        if(jproc .eq. 26) then
            call gec_ecmmfinit
        endif
        if(jproc .eq. 27) then
            call gec_ecttfinit
        endif
        if(jproc .eq. 28) then
            call gecb_ecbeefinit
        endif
        if(jproc .eq. 29) then
            call gecb_ecbmmfinit
        endif
        if(jproc .eq. 30) then
            call gecb_ecbttfinit
        endif
        if(jproc .eq. 31) then
            call geb_ebeefinit
        endif
        if(jproc .eq. 32) then
            call geb_ebmmfinit
        endif
        if(jproc .eq. 33) then
            call geb_ebttfinit
        endif
        if(jproc .eq. 34) then
            call gebb_ebbeefinit
        endif
        if(jproc .eq. 35) then
            call gebb_ebbmmfinit
        endif
        if(jproc .eq. 36) then
            call gebb_ebbttfinit
        endif
        if(jproc .eq. 37) then
            call get_eteefinit
        endif
        if(jproc .eq. 38) then
            call get_etmmfinit
        endif
        if(jproc .eq. 39) then
            call get_etttfinit
        endif
        if(jproc .eq. 40) then
            call getb_etbeefinit
        endif
        if(jproc .eq. 41) then
            call getb_etbmmfinit
        endif
        if(jproc .eq. 42) then
            call getb_etbttfinit
        endif
        return
        end
