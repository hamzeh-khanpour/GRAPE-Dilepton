*     File ec_ectt/inclec_ectt1.h : Sat Mar 18 19:45:04 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*

      include 'inclc.h'
      include 'inclg.h'

      parameter (loutgo =  2, lincom =  1)
      parameter (lantip = -1, lprtcl =  1)
      parameter (lscalr =  1)
      parameter (lepexa =  2, lepexv =  3)
      parameter (lepina =  4, lepinv =  4)
      parameter (lextrn =  2, lintrn =  4)

* table of amplitudes
      parameter (nec_ecttgraph =35)
      parameter (nec_ecttextn  =6)


* number of color base
      parameter (ncbase =1)

* number of all helicity states
      parameter (lec_ecttag    =64)

      common /amprck/kmngr,  kmnext, kmlag
      common /amprcs/smpref, smproc
      character*80   smpref, smproc

      parameter (nec_ecttgrpsq = nec_ecttgraph*nec_ecttgraph)
      common /aec_ecttmslct/jselg(nec_ecttgraph),jgraph,jgluon,jhiggs

* Color string information
      common /aec_ecttmcsti/ kmcbas, kmcbmx, kmcstr(2), icinfo(6), icolst
      common /aec_ecttmgrph/agcwrk(0:lec_ecttag-1),agc(0:lec_ecttag-1),
     &              aprop,ancp(nec_ecttgraph),ansp(0:nec_ecttgraph)
     &             ,cfmtx
      common /aec_ecttmgrpi/igraph(nec_ecttgraph)
      complex*16 agc, agcwrk, aprop

* Momenta of external particles
      common /aec_ecttmextr/peec_ectt1(4),peec_ectt2(4),peec_ectt3(4),
     &               peec_ectt4(4),peec_ectt5(4),peec_ectt6(4),
     &               prod(nec_ecttextn, nec_ecttextn)

      common /aec_ecttmextp/ptec_ectt2q,exec_ectt2q,cfec_ectt2q,
     &               ptec_ectt4m,exec_ectt4m,cfec_ectt4m,ptec_ectt9q,
     &               exec_ectt9q,cfec_ectt9q,ptec_ectt17m,exec_ectt17m,
     &               cfec_ectt17m,ptec_ectt33o,exec_ectt33o,
     &               cfec_ectt33o,ptec_ectt62o,exec_ectt62o,
     &               cfec_ectt62o

      real*8     ptec_ectt2q(4,2), exec_ectt2q(1)
      complex*16 cfec_ectt2q(2,2)
      real*8     ptec_ectt4m(4,2), exec_ectt4m(1)
      complex*16 cfec_ectt4m(2,2)
      real*8     ptec_ectt9q(4,2), exec_ectt9q(1)
      complex*16 cfec_ectt9q(2,2)
      real*8     ptec_ectt17m(4,2), exec_ectt17m(1)
      complex*16 cfec_ectt17m(2,2)
      real*8     ptec_ectt33o(4,2), exec_ectt33o(1)
      complex*16 cfec_ectt33o(2,2)
      real*8     ptec_ectt62o(4,2), exec_ectt62o(1)
      complex*16 cfec_ectt62o(2,2)

* Normalization
      common /aec_ecttmdbgg/fknorm,fkcall,nkcall

* Calculated table of amplitudes
      common /aec_ecttmatbl/av, lt, indexg
      complex*16 av(0:lec_ecttag-1)
      integer    lt(0:nec_ecttextn), indexg(nec_ecttextn)

* Spin average
      common/aec_ecttmspin/aspin,aident,jhs(nec_ecttextn),
     &               jhe(nec_ecttextn),jcpol(nec_ecttextn),
     &               kperm(nec_ecttextn)
