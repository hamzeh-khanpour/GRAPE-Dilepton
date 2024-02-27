       FUNCTION MINOR(A,N,I,J)
       REAL*8 MINOR,A(100),B(100),DETERM
       IF(N.GT.10) THEN
        WRITE(6,*)' # ORDER IS GREATER THAN TEN ',N
        RETURN
       END IF
       DO 1 I1=1,N-1
        II1=I1
        IF(I1.GE.I)II1=I1+1
        DO 2 I2=1,N-1
         II2=I2
         IF(I2.GE.J)II2=I2+1
2        B(I1+(N-1)*(I2-1))=A(II1+N*(II2-1))
1       CONTINUE
       MINOR=DETERM(B,N-1)
       RETURN
       END
