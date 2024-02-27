      subroutine ktolab(p1,p2,p3,p4, p5,p6)
*
      implicit real*8(a-h,o-z)
      real*8 p1(4),p2(4),p3(4),p4(4),p5(4),p6(4)
      real*8 pk2l(4),pl2k(4),rk2l(3,3),rl2k(3,3),pt1(4)
*
      do 1 i=1,4
       pk2l(i)= p2(i)
       pl2k(i)=-p2(i)
1     continue
       pl2k(4)= p2(4)
*
      call pboost(p1,pl2k, pt1)
      call rotmtx(pt1, rl2k)
      call minvr2(rl2k, rk2l)
*
      call mvmult(rk2l,p3, p5)
      call pboost(p5,pk2l, p5)
      call mvmult(rk2l,p4, p6)
      call pboost(p6,pk2l, p6)
*
      return
      end
