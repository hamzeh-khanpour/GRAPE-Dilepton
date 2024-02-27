       SUBROUTINE MTXINV(A,N)
       REAL*8 A(100),B(100),DET,DETERM,MINOR
       IF(N.GT.10) THEN
        WRITE(6,*)' # ORDER IS GREATER THAN TEN ',N
        RETURN
       END IF
       CALL MTXCPY(A,B,N)
       DET=DETERM(A,N)
       DO 2 I3=1,N
       DO 2 I2=1,N
2        A(I2+N*(I3-1))=(1-MOD(I3+I2,2)*2)*MINOR(B,N,I3,I2)/DET
       RETURN
       END
