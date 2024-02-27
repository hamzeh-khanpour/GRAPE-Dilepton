       FUNCTION DETERM(A,N)
       REAL*8 DETERM,A(100),B(10,10),SAVE
       IF(N.GT.10) THEN
        WRITE(6,*)' # ORDER IS GREATER THAN TEN ',N
        RETURN
       END IF
       DO 1 I12=1,N
       DO 1 I11=1,N
1        B(I11,I12)=A(I11+N*(I12-1))
       DETERM=1
       DO 2 I2=1,N
        IF(B(I2,I2).EQ.0)THEN
         DO 5 I5=I2,N
          IF(B(I2,I5).EQ.0)THEN
           DETERM=0
           RETURN
          ELSE
           DO 6 I6=I2,N
            SAVE=B(I6,I5)
            B(I6,I5)=B(I5,I2)
            B(I5,I2)=SAVE
6           CONTINUE
           DETERM=-DETERM
          END IF
5         CONTINUE
        END IF
        DETERM=DETERM*B(I2,I2)
        IF(I2.EQ.N) GOTO 2
        DO 3 I3=I2+1,N
        DO 3 I4=I2+1,N
3         B(I3,I4)=B(I3,I4)-B(I3,I2)*B(I2,I4)/B(I2,I2)
2       CONTINUE
       RETURN
       END
