       SUBROUTINE MMMULT(R1,R2,R3)
       REAL*8 R1(3,3),R2(3,3),R3(3,3),RDUMMY(3,3)
       CALL MEMCLR(RDUMMY,9,0)
       DO 1 I=1,3
       DO 2 J=1,3
       DO 3 K=1,3
3       RDUMMY(I,J)=RDUMMY(I,J)+R1(I,K)*R2(K,J)
2      CONTINUE
1      CONTINUE
       DO 4 I=1,3
       DO 5 J=1,3
5       R3(I,J)=RDUMMY(I,J)
4      CONTINUE
C
       RETURN
       END
