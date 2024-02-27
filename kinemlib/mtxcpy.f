       SUBROUTINE MTXCPY(A,B,N)
       REAL*8 A(100),B(100)
       IF(N.GT.10) THEN
        WRITE(6,*)' # ORDER IS GREATER THAN TEN ',N
        RETURN
       END IF
       DO 1 I12=1,N
       DO 1 I11=1,N
1       B(I11+N*(I12-1))=A(I11+N*(I12-1))
       RETURN
       END
