*nt_inc.h
      integer*4           nt_jproc, nt_merge, nt_Nextn, nt_Ngen
       common /NT_PRC_i4/ nt_jproc, nt_merge, nt_Nextn, nt_Ngen
      character*(*) NT_PRC_i4
       parameter (  NT_PRC_i4 ='jproc:i*4,merge:i*4,nextn:i*4,ngen:i*4')
      double precision    nt_P1_lab, nt_P2_lab, nt_xsec(2)
       common /NT_PRC_r8/ nt_P1_lab, nt_P2_lab, nt_xsec
      character*(*) NT_PRC_r8
       parameter (  NT_PRC_r8 = 'P1_lab:r*8, P2_lab:r*8, xsec(2):r*8' )
      character   cmaxisr*1
       parameter (cmaxisr ='1')
      integer*4           nt_Nisr
       common /NT_EVT_i4/ nt_Nisr
      character*(*) NT_EVT_i4
       parameter (  NT_EVT_i4 = 'nisr[0,'//cmaxisr//']:i' )
C      integer     nextn
C      character  cnextn*1
C        parameter ( nextn = 6 )   !!!
C        parameter (cnextn ='6')   !!!
      double precision  vec_px, vec_py, vec_pz, vec_e, vec_m
       common /NT_VEC/  vec_px, vec_py, vec_pz, vec_e, vec_m
      character*(*) NT_VEC
      parameter (   NT_VEC   = 'vec_px:r*8'
     &                      //',vec_py:r*8'
     &                      //',vec_pz:r*8'
     &                      //',vec_e:r*8'
     &                      //',vec_m:r*8'
     &          )
      integer*4        vec_kf
       common /NT_KF/  vec_kf
      character*(*) NT_KF
       parameter (  NT_KF = 'vec_kf[-6000,6000]:i' )
      double precision  nt_vec_isr(4)
       common /NT_ISR/  nt_vec_isr
      character*(*)  NT_ISR
      parameter (    NT_ISR = 'vec_isr(4):r*8' )
      integer    N_max
      character cN_max*2
        parameter ( N_max = 63 )
        parameter (cN_max ='63')
      integer*4                Npy, sta(N_max), kf(N_max), mot(N_max)
      real*4                   px(N_max),py(N_max),pz(N_max),pe(N_max)
     &                        ,pm(N_max)
       common /NT_PYTHIA/ Npy, px, py, pz, pe, pm, kf, sta, mot
      character*(*) PYTHIA
      parameter (   PYTHIA =
     &                       'npy[0,'//cN_max//']:i'
     &                    //',px(npy):r*4'
     &                    //',py(npy):r*4'
     &                    //',pz(npy):r*4'
     &                    //',pe(npy):r*4'
     &                    //',pm(npy):r*4'
     &                    //',kf(npy)[-6000,6000]:i'
     &                    //',sta(npy)[0,41]:i'
     &                    //',mot(npy)[0,'//cN_max//']:i'
     &          )
      integer         NTID, LREC, IQ10
      character*(*)   NT_NAME, NT_TITLE
       parameter(NTID=1, NT_NAME='grp.rz', NT_TITLE = 'GRAPE Ntuple'
     &          ,LREC=4096, IQ10=1000000)
