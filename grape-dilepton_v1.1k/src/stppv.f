*     File stppv.f : Sat Mar 18 19:44:57 2000
*     Generated by GRACE (Fortran source code generator)
*         Ver.2.1.5(2) 1998/05/11
*     (c)copyright 1990-1998 Minami-Tateya Group, Japan
*-----------------------------------------------------------------------
      subroutine stppv(
     &        kf1,lh1,am1,pf1,idir1,vn1,ex1v,cf1v,pt1v,
     &        kf2,lh2,am2,pf2,idir2,vn2,ex2v,cf2v,pt2v,
     &        kf3,lh3,am3,pf3,idir3,vn3,ex3v, eq,
     &        capp, lt,av)
      implicit real*8(a-h,o-z)
      real*8 pf1(4),pf2(4),pf3(4)
      real*8 ex1v(1),ex2v(1),ex3v
      complex*16 cf1v(2,lh1),cf2v(2,lh2)
      real*8 pt1v(4,lh1),pt2v(4,lh2),eq(4,lh3)
      complex*16 capp(2)
      parameter( ltsize=20, lasize=1024 )
      integer lt(0:3), lts(0:3)
      complex*16 avv(lasize),avs(lasize),avx(lasize),av(lh1*lh2*lh3)
c     complex*16 av1(0:lasize)
*-----------------------------------------------------------------------
      call smffv(lh1,lh2,lh3,ex1v,ex2v,am1,am2,capp,cf1v,cf2v,
     &           pt1v,pt2v,eq,lt,avv)
      call smffs(lh1,lh2,ex1v,ex2v,am1,am2,capp,cf1v,cf2v,
     &           pt1v,pt2v,lts,avs)
      do 1002 ih3  = 1 , lh3
         xxx = eq(4,ih3)*(pf1(4)+pf2(4)) - eq(1,ih3)*(pf1(1)+pf2(1))
     &       - eq(2,ih3)*(pf1(2)+pf2(2)) - eq(3,ih3)*(pf1(3)+pf2(3))
      do 1001 ihs = 1 , lh1*lh2
         ih = (lh1*lh2)*(ih3-1) + ihs
         avx(ih) = avs(ihs)*xxx
 1001 continue
 1002 continue
      ge  = 1.d0/(1.0d0-vn3/0.71d0)**2
      pmu = 2.79284739d0       ! Magnetic moment
      amm = pmu - 1.d0         ! Anomalous magnetic moment
      gm  = pmu*ge
      xu  = 1.0d0 - vn3/(2.0d0*am1)**2
      do 2001 ih = 1 , lh1*lh2*lh3
         av(ih) = pmu * ge * avv(ih)
     &          - amm/xu*ge /(2.0d0*am1)* avx(ih)
 2001 continue
      return
      end
