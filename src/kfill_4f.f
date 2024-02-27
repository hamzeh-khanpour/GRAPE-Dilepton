* FILE kfill
************************************************************************
      SUBROUTINE KFILL_4f(NEXTRN, ndim, X, PE, PP, ANS)
      IMPLICIT REAL* 8(A-H,O-Z)
************************************************************************
      INTEGER NEXTRN
      PARAMETER ( MXDIM = 50 )
      REAL*8  X(MXDIM)
      REAL*8  PE(4,NEXTRN), PP(NEXTRN,NEXTRN)
      INCLUDE 'inclk.h'
      include 'inclp_4f.h'
      common/kmcntl_4f/iresns(4),icos3,icosq3,icos5,isr,iswap,ident
     &                ,iphi6,ieeee,i34,itag,isym
      real*8   mass45_2,mass45, mass56_2,mass56
*=======================================================================
*          Fill Histograms and Scatter plots
*=======================================================================
        anss=ans
        if (anss.lt.1.d-10)  anss=0
        DO 40 I = 1, NDIM
           CALL XHFILL( I, X(I), ANSs )
   40   CONTINUE
           CALL XHFILL(ndim+1,pe(4,3) , ANSs )
           CALL XHFILL(ndim+2,pe(4,4) , ANSs )
           CALL XHFILL(ndim+3,pe(4,5) , ANSs )
           CALL XHFILL(ndim+4,pe(4,6) , ANSs )
	   cos3=pe(3,3)/sqrt(pe(1,3)**2+pe(2,3)**2+pe(3,3)**2)
           CALL XHFILL(ndim+5,cos3, ANSs )
	   cos4=pe(3,4)/sqrt(pe(1,4)**2+pe(2,4)**2+pe(3,4)**2)
           CALL XHFILL(ndim+6,cos4, ANSs )
	   cos5=pe(3,5)/sqrt(pe(1,5)**2+pe(2,5)**2+pe(3,5)**2)
           CALL XHFILL(ndim+7,cos5, ANsS )
	   cos6=pe(3,6)/sqrt(pe(1,6)**2+pe(2,6)**2+pe(3,6)**2)
           CALL XHFILL(ndim+8,cos6, ANSs )
           mass45_2 = amass2(4)+amass2(5)+pp(4,5)*2d0
           if (mass45_2 .LT. 0) then
             write(6,*) '!!!Warning in kfill_4f!!!'
             write(6,*) '  ---> Mass of Particle4-5 squared =',mass45_2
     &                 ,' GeV < 0 (negative!)'
             write(6,*) '  ---> X =', (X(i),i=1,NDIM)
           endif
           mass45 = sqrt(  max( mass45_2, 0.D0 )  )
           CALL XHFILL(ndim+ 9, mass45, ANSs )
           mass56_2 = amass2(5)+amass2(6)+pp(5,6)*2d0
           if (mass56_2 .LT. 0) then
             write(6,*) '!!!Warning kfill_4f!!!'
             write(6,*) '  ---> Mass of particle5-6 squared =',mass56_2
     &                 ,' GeV < 0 (negative!)'
             write(6,*) '  ---> X = ', (X(i),i=1,NDIM)
           endif
           mass56 = sqrt(  max( mass56_2, 0.D0 )  )
           CALL XHFILL(ndim+10, mass56, ANSs)
      RETURN
      END
