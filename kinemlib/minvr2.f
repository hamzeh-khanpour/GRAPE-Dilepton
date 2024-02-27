       SUBROUTINE MINVR2(M1,M2)
       REAL*8 M1(3,3),M2(3,3),M3(3,3)
       DO 1 I=1,3
       DO 1 J=1,3
        M3(I,J)=M1(J,I)
1      CONTINUE
       DO 2 I=1,3
       DO 2 J=1,3
        M2(I,J)=M3(I,J)
2      CONTINUE
       RETURN
       END
