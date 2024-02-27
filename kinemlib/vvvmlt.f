       SUBROUTINE VVVMLT(P1,P2,P3)
       REAL*8 P1(0:3),P2(0:3),P3(0:3)
       P3(0)=0
       P3(1)=P1(2)*P2(3)-P1(3)*P2(2)
       P3(2)=P1(3)*P2(1)-P1(1)*P2(3)
       P3(3)=P1(1)*P2(2)-P1(2)*P2(1)
       RETURN
       END
