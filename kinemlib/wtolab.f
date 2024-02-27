       SUBROUTINE WTOLAB(P1,P2,P3, P4,P5)
* Debug 04/Jan/95 Y.Kurihara: To treat P3(1)**2+P3(2)**2=0
* Debug 04/Jan/95 Y.Kurihara: >1.d-15 --> >1.d-10
* Debug 25/Jan/95 Y.Kurihara: change to use rotmtx 
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION P1(4),P2(4),P3(4),P4(4),P5(4)
       DIMENSION ROT1(3,3),pt(4)
C DETERMIN ROTATION MATRIX 1
       pt(1)=p3(1)
       pt(2)=p3(2)
       pt(3)=p3(3)
       pt(4)=p3(4)
       call rotmtx(pt, rot1)
       CALL MINVR2(ROT1,ROT1)
C
       CALL MVMULT(ROT1,P1,P4)
       CALL MVMULT(ROT1,P2,P5)
       CALL PBOOST(P4,P3,P4)
       CALL PBOOST(P5,P3,P5)
C
       RETURN
       END
