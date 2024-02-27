       SUBROUTINE PBOOST(P,Q,PB)
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION P(4),Q(4),PB(4)
C
       PQ=P(1)*Q(1)+P(2)*Q(2)+P(3)*Q(3)
       Q2=Q(1)**2+Q(2)**2+Q(3)**2
       Q1=SQRT(Q2)
C
       AM=SQRT( (Q(4)-Q1)*(Q(4)+Q1) )
        F=((Q(4)-AM)*PQ/Q2+P(4))/AM
       PB(1)= P(1)+Q(1)*F
       PB(2)= P(2)+Q(2)*F
       PB(3)= P(3)+Q(3)*F
       PB(4)=(P(4)*Q(4)+PQ)/AM
       RETURN
       END
