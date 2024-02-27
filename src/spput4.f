      subroutine spput4(nreg, nextn, lag, p, agc)
      implicit real*8 (a-h,o-z)
      integer  nreg, nextn, lag
      real*8   p(4,nextn)
      complex*16  agc(0:lag-1)
      parameter (MAXREG = 6)
      parameter (mxextn = 10)
      parameter (mxlag  = 6561)   
      complex*16  grcevt
      common /grc4sp/ answrk(0:MAXREG-1)
     &               ,wrkp(4, mxextn,   0:MAXREG-1)
     &               ,grcevt(0:mxlag-1, 0:MAXREG-1)
      dimension  cratio(0:MAXREG-1)
      common /sp4vec/ vec(4,mxextn)
      real*8    drn
      external  drn
*-----------------------------------------------------------------------
* ------- selection ------------------
      ireg = 0
      if( nreg .gt. 1 ) then
        allsum = 0.0d0
        do 1000 i = 0, nreg-1
           allsum = allsum + answrk(i)
 1000   continue
        tmpsum = 0.0d0
        do 2000 i = 0, nreg-1
           tmpsum = tmpsum + answrk(i)
           cratio(i) = tmpsum/allsum
 2000   continue
        cran = drn(idummy)
        do 3000 i = 0, nreg-1
           if( cratio(i) .gt. cran ) then
               ireg = i
               goto 4000
           endif
 3000   continue
 4000   continue
      endif
      do 5000 i = 1 , nextn
      do 5000 j = 1 , 4
         p(j,i) = wrkp(j,i, ireg)
 5000 continue
      do 6000 i = 1 , nextn
      do 6000 j = 1 , 4
             vec(j,i) = p(j,i)
 6000     continue
      do 7000 ih = 0, lag-1
         agc(ih) = grcevt(ih, ireg)
 7000 continue
      return
      end
