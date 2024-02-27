       SUBROUTINE MVMULT(R,P1,P2)
       REAL*8 R(3,3),P1(4),P2(4),PTMP(4)
C
       PTMP(1)=P1(1)*R(1,1)+P1(2)*R(1,2)+P1(3)*R(1,3)
       PTMP(2)=P1(1)*R(2,1)+P1(2)*R(2,2)+P1(3)*R(2,3)
       PTMP(3)=P1(1)*R(3,1)+P1(2)*R(3,2)+P1(3)*R(3,3)
       PTMP(4)=P1(4)
C
       P2(1)=PTMP(1)
       P2(2)=PTMP(2)
       P2(3)=PTMP(3)
       P2(4)=PTMP(4)
C
       RETURN
       END
