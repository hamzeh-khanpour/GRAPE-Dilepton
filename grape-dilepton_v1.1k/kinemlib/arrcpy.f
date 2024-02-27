       SUBROUTINE ARRCPY(A,B,N)
       REAL*8 A(100),B(100)
       IF(N.GT.100) THEN
        WRITE(6,*)' # ORDER IS GREATER THAN 100 ',N
        RETURN
       END IF
       DO 1 I11=1,N
1        A(I11)=B(I11)
       RETURN
       END
