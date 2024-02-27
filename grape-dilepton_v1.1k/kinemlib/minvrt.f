       SUBROUTINE MINVRT(R1,R2)
       REAL*8 R1(9),R2(9)
       REAL*8 RTMP1(9)
C
       DO 1 I=1,9
1       RTMP1(I)=R1(I)
C
C      CALL DINV(3,RTMP1,3,RTMP2,ICON)
       CALL MTXINV(RTMP1,3)
       DO 2 I=1,9
2       R2(I)=RTMP1(I)
C
       RETURN
       END
