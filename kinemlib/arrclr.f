       SUBROUTINE ARRCLR(A,N)
       REAL*8 A(100)
       IF(N.GT.100) THEN
        WRITE(6,*)' # ORDER IS GREATER THAN 100 ',N
        RETURN
       END IF
       DO 1 I12=1,N
1        A(I12)=0
       RETURN
       END
