*     File spdetc.f : Sat Mar 18 19:44:57 2000
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
************************************************************************
      subroutine spdetc
      implicit real*8(a-h,o-z)
      common /amjprc/jproc
        if(jproc .eq. 1) then
            call spep_epeedetc
        endif
        if(jproc .eq. 2) then
            call spep_epmmdetc
        endif
        if(jproc .eq. 3) then
            call spep_epttdetc
        endif
        if(jproc .eq. 4) then
            call spep_eXeedetc
        endif
        if(jproc .eq. 5) then
            call spep_eXmmdetc
        endif
        if(jproc .eq. 6) then
            call spep_eXttdetc
        endif
        if(jproc .eq. 7) then
            call speu_eueedetc
        endif
        if(jproc .eq. 8) then
            call speu_eummdetc
        endif
        if(jproc .eq. 9) then
            call speu_euttdetc
        endif
        if(jproc .eq. 10) then
            call speub_eubeedetc
        endif
        if(jproc .eq. 11) then
            call speub_eubmmdetc
        endif
        if(jproc .eq. 12) then
            call speub_eubttdetc
        endif
        if(jproc .eq. 13) then
            call sped_edeedetc
        endif
        if(jproc .eq. 14) then
            call sped_edmmdetc
        endif
        if(jproc .eq. 15) then
            call sped_edttdetc
        endif
        if(jproc .eq. 16) then
            call spedb_edbeedetc
        endif
        if(jproc .eq. 17) then
            call spedb_edbmmdetc
        endif
        if(jproc .eq. 18) then
            call spedb_edbttdetc
        endif
        if(jproc .eq. 19) then
            call spes_eseedetc
        endif
        if(jproc .eq. 20) then
            call spes_esmmdetc
        endif
        if(jproc .eq. 21) then
            call spes_esttdetc
        endif
        if(jproc .eq. 22) then
            call spesb_esbeedetc
        endif
        if(jproc .eq. 23) then
            call spesb_esbmmdetc
        endif
        if(jproc .eq. 24) then
            call spesb_esbttdetc
        endif
        if(jproc .eq. 25) then
            call spec_eceedetc
        endif
        if(jproc .eq. 26) then
            call spec_ecmmdetc
        endif
        if(jproc .eq. 27) then
            call spec_ecttdetc
        endif
        if(jproc .eq. 28) then
            call specb_ecbeedetc
        endif
        if(jproc .eq. 29) then
            call specb_ecbmmdetc
        endif
        if(jproc .eq. 30) then
            call specb_ecbttdetc
        endif
        if(jproc .eq. 31) then
            call speb_ebeedetc
        endif
        if(jproc .eq. 32) then
            call speb_ebmmdetc
        endif
        if(jproc .eq. 33) then
            call speb_ebttdetc
        endif
        if(jproc .eq. 34) then
            call spebb_ebbeedetc
        endif
        if(jproc .eq. 35) then
            call spebb_ebbmmdetc
        endif
        if(jproc .eq. 36) then
            call spebb_ebbttdetc
        endif
        if(jproc .eq. 37) then
            call spet_eteedetc
        endif
        if(jproc .eq. 38) then
            call spet_etmmdetc
        endif
        if(jproc .eq. 39) then
            call spet_etttdetc
        endif
        if(jproc .eq. 40) then
            call spetb_etbeedetc
        endif
        if(jproc .eq. 41) then
            call spetb_etbmmdetc
        endif
        if(jproc .eq. 42) then
            call spetb_etbttdetc
        endif
        return
        end
