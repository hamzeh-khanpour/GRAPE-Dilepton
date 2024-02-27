*     File ep_eXtt/inclep_eXtt1.h : Sat Mar 18 19:45:45 2000
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
      parameter (nep_eXttgraph =8)
      parameter (nep_eXttextn  =6)


* number of all helicity states
      parameter (lep_eXttag    =64)

      common /amprck/kmngr,  kmnext, kmlag
      common /amprcs/smpref, smproc
      character*80   smpref, smproc

      parameter (nep_eXttgrpsq = nep_eXttgraph*nep_eXttgraph)
      common /aep_eXttmslct/jselg(nep_eXttgraph),jgraph,jgluon,jhiggs

* Color string information
      common /aep_eXttmcsti/ kmcbas, kmcbmx, icinfo(6), icolst
      common /aep_eXttmgrph/agcwrk(0:lep_eXttag-1),agc(0:lep_eXttag-1),
     &              aprop,ancp(nep_eXttgraph),ansp(0:nep_eXttgraph)
     &             ,cfmtx
      common /aep_eXttmgrpi/igraph(nep_eXttgraph)
      complex*16 agc, agcwrk, aprop

* Momenta of external particles
      common /aep_eXttmextr/peep_eXtt1(4),peep_eXtt2(4),peep_eXtt3(4),
     &               peep_eXtt4(4),peep_eXtt5(4),peep_eXtt6(4),
     &               prod(nep_eXttextn, nep_eXttextn)

      common /aep_eXttmextp/ptep_eXtt2v,exep_eXtt2v,cfep_eXtt2v,
     &               ptep_eXtt4m,exep_eXtt4m,cfep_eXtt4m,ptep_eXtt9x,
     &               exep_eXtt9x,cfep_eXtt9x,ptep_eXtt17m,exep_eXtt17m,
     &               cfep_eXtt17m,ptep_eXtt33o,exep_eXtt33o,
     &               cfep_eXtt33o,ptep_eXtt62o,exep_eXtt62o,
     &               cfep_eXtt62o

      real*8     ptep_eXtt2v(4,2), exep_eXtt2v(1)
      complex*16 cfep_eXtt2v(2,2)
      real*8     ptep_eXtt4m(4,2), exep_eXtt4m(1)
      complex*16 cfep_eXtt4m(2,2)
      real*8     ptep_eXtt9x(4,2), exep_eXtt9x(1)
      complex*16 cfep_eXtt9x(2,2)
      real*8     ptep_eXtt17m(4,2), exep_eXtt17m(1)
      complex*16 cfep_eXtt17m(2,2)
      real*8     ptep_eXtt33o(4,2), exep_eXtt33o(1)
      complex*16 cfep_eXtt33o(2,2)
      real*8     ptep_eXtt62o(4,2), exep_eXtt62o(1)
      complex*16 cfep_eXtt62o(2,2)

* Normalization
      common /aep_eXttmdbgg/fknorm,fkcall,nkcall

* Calculated table of amplitudes
      common /aep_eXttmatbl/av, lt, indexg
      complex*16 av(0:lep_eXttag-1)
      integer    lt(0:nep_eXttextn), indexg(nep_eXttextn)

* Spin average
      common/aep_eXttmspin/aspin,aident,jhs(nep_eXttextn),
     &               jhe(nep_eXttextn),jcpol(nep_eXttextn),
     &               kperm(nep_eXttextn)
