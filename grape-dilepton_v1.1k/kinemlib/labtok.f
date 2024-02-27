      subroutine labtok(p1,p2,p3,p4, p5,p6)
*
      implicit real*8(a-h,o-z)
      real*8 p1(4),p2(4),p3(4),p4(4),p5(4),p6(4)
      real*8 pl2k(4),pt1(4),rl2k(3,3)
*
      do 1 i=1,4
       pl2k(i)=-p2(i)
1     continue
       pl2k(4)= p2(4)
*
      call pboost(p1,pl2k, pt1)
      call rotmtx(pt1, rl2k)
*
      call pboost(p3,pl2k, p5)
      call mvmult(rl2k,p5, p5)
      call pboost(p4,pl2k, p6)
      call mvmult(rl2k,p6, p6)
*
      return
      end
