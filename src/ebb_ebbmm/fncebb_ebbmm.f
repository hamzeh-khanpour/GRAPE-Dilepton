*     File ebb_ebbmm/fncebb_ebbmm.f : Sat Mar 18 19:45:06 2000
*
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
*
************************************************************************
      function fncebb_ebbmm(x)
      implicit real*8(a-h,o-z)
      real*8   fncebb_ebbmm

      include 'inclebb_ebbmm1.h'
      include 'inclk.h'

      parameter ( mxdim = 50 )
      common / loop0 / loop
      common / bparm1 / xl(mxdim),xu(mxdim),ndim,nwild,
     &                 ig(mxdim),ncall
      common / bparm2 / acc1,acc2,itmx1,itmx2

      real*8   x(mxdim) 

      integer         JFLAG,IBASES 
       common /BASE0/ JFLAG,IBASES

      real*8     ansm, ans

      real*8     xx(mxdim),p(4,nebb_ebbmmextn),pp(nebb_ebbmmextn,nebb_ebbmmextn)
      common /sp4vec/ vec(4,nebb_ebbmmextn)
*=======================================================================
*          initialization
*=======================================================================

      ansum = 0.0d0

      do  5 i = 1, ndim
         xx(i) = x(i)
    5 continue

      nreg  = 1
      dft   = 0.d0
*=======================================================================
*          kinematics
*=======================================================================

      kreg = 0 
      do 1000 ireg = 1 , mxreg

         if(ireg .gt. nreg) go to 1000

        call kinem_4f(nebb_ebbmmextn, xx, p, pp, yacob,nreg,ireg,jump)

*-----------------------------------------------------------------------
*          reset the temporal buffer for the region 1
*-----------------------------------------------------------------------
      if(ireg .eq. 1) then
          dft = 0.d0
          do 180 k = 1, nebb_ebbmmextn
          do 180 j = 1, 4
             vec(j,k) = 0.d0
  180     continue
      endif

         if(jump .ne. 0) go to 1000

*-----------------------------------------------------------------------
*          for user's cut
*-----------------------------------------------------------------------

*--> call usrcut(jump)
*--> if(jump .ne. 0) goto 1000

*-----------------------------------------------------------------------
*          four momenta of external particles
*-----------------------------------------------------------------------

      do 20 i = 1, 4
*     1: initial b-bar mass=amdq(3)
      peebb_ebbmm1(i) = p(i, kperm(1))

*     2: initial positron mass=amlp(1)
      peebb_ebbmm2(i) = p(i, kperm(2))

*     3: final b-bar mass=amdq(3)
      peebb_ebbmm3(i) = p(i, kperm(3))

*     4: final positron mass=amlp(1)
      peebb_ebbmm4(i) = p(i, kperm(4))

*     5: final muon mass=amlp(2)
      peebb_ebbmm5(i) = p(i, kperm(5))

*     6: final anti-muon mass=amlp(2)
      peebb_ebbmm6(i) = p(i, kperm(6))

   20 continue

*-----------------------------------------------------------------------
*          inner products of momenta of external particles
*-----------------------------------------------------------------------

      do 30 j = 1, nebb_ebbmmextn
      do 30 i = 1, nebb_ebbmmextn
        prod(i, j) = pp(kperm(i), kperm(j))
   30 continue

*=======================================================================
*          amplitude calculation
*=======================================================================
        call aebb_ebbmmkptbl
*       ===================
        call kinemp_4f( ireg )
*       ===================
        call aebb_ebbmmmptbl
*       ===================
        call aebb_ebbmmmpsum(ansm)
*       ===================

        fknorm = yacob*aspin
        ans    = ansm*fknorm/aident
        ansum  = ansum + ans 

        if (IBASES.EQ.0) then 
          call spstr4(kreg, nebb_ebbmmextn, lebb_ebbmmag, ans, p, agc) 
        endif

*-----------------------------------------------------------------------
*          save four momenta and probabilities of the region 1
*-----------------------------------------------------------------------
      if(ireg .eq. 1) then
          dft = ans
          do 420 k = 1, nebb_ebbmmextn
          do 420 j = 1, 4
             vec(j,k) = p(j,k)
  420     continue
      endif

*=======================================================================
*          fill histograms and scatter plots
*=======================================================================
      call kfill_4f(nebb_ebbmmextn, ndim, x, p, pp, ans*aident)
*=======================================================================
*          update summary table
*=======================================================================

        call bsgetw(wgt)
        ansp(0) = ansp(0) + wgt*ans
        do 60 igr = 1, jgraph
          ansp(igr)=ansp(igr) + wgt*yacob*aspin*ancp(igr)
   60   continue
        nkcall = nkcall + 1
        if(nkcall.gt.10000) then
          nkcall = nkcall - 10000
          fkcall = fkcall + 10000
        endif

        kreg = kreg + 1 
 1000 continue 

      if (IBASES.EQ.0) then 
        call spput4(kreg, nebb_ebbmmextn, lebb_ebbmmag, p, agc) 
      endif

      fncebb_ebbmm = ansum

*-----------------------------------------------------------------------
*          put the final 4 vectors into the arrays vec()
*-----------------------------------------------------------------------
*     if(func .gt. 0.d0) then
*         if(dft/func .lt. drn(dum)) then
*             do 850 k = 1, nebb_ebbmmextn
*             do 850 j = 1, 4
*                vec(j,k) = p(j,k)
* 850         continue
*         endif
*     endif

      return
      end