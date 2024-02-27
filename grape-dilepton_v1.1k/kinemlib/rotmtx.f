       SUBROUTINE ROTMTX(PAX, ROT)
* Debug 27/DEC/94 Y.Kurihara: To treat PAX(1)**2+PAX(2)**2=0
* Debug 04/Jan/95 Y.Kurihara: den>1.d-15 -->den>1.d-10 
*       19/Jan/95 Y.Kurihara: To treat phi rotation    
*
       IMPLICIT REAL* 8(A-H,O-Z)
       REAL*8 PAX(4)
       REAL*8 ROT(3,3),ROTTMP(3,3)
       den=SQRT(PAX(1)*PAX(1)+PAX(2)*PAX(2))
       if(den.le.1.d-10) then
         CALL MEMCLR(ROT,9,0)
         rot(1,1)=sign(1.d0,pax(3))
         rot(2,2)=sign(1.d0,pax(3))
         rot(3,3)=sign(1.d0,pax(3))
         return
       end if
       COSA=PAX(1)/den
       SINA=PAX(2)/den
       CALL MEMCLR(ROTTMP,9,0)
       ROTTMP(1,1)= COSA
       ROTTMP(1,2)= SINA
       ROTTMP(2,1)=-SINA
       ROTTMP(2,2)= COSA
       ROTTMP(3,3)= 1
       CALL MVMULT(ROTTMP,PAX,PAX)
       CALL MEMCLR(ROT,9,0)
       COSB=PAX(3)/SQRT(PAX(1)*PAX(1)+PAX(3)*PAX(3))
       SINB=PAX(1)/SQRT(PAX(1)*PAX(1)+PAX(3)*PAX(3))
       ROT(1,1)= COSB
       ROT(1,3)=-SINB
       ROT(2,2)= sign(1.d0,pax(3))
       ROT(3,1)= SINB
       ROT(3,3)= COSB
       CALL MMMULT(ROT,ROTTMP,ROT)
       RETURN
       END
