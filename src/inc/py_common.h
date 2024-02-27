C ********************************************************************
      integer*4   MSTP(200),MSTI(200)    !@ Event Information
      real*8      PARP(200),PARI(200)    !@   on the latest one
       COMMON/PYPARS/MSTP,PARP,MSTI,PARI

      integer*4   N,NPAD, K(4000,5)      !@ Event Record
      real*8      P(4000,5), V(4000,5)   !@   on the current one
       common /PYJETS/N,NPAD,K,P,V

      integer*4   MSEL,MSELPD,MSUB(500)  !@ Process control
     &           ,KFIN(2,-40:40)         !@   variables
      real*8      CKIN(200)              !@     by user
       COMMON/PYSUBS/MSEL,MSELPD,MSUB,KFIN,CKIN

      integer*4   MINT(400)              !@ Variables
      real*8      VINT(400)              !@   used internally
       COMMON/PYINT1/MINT,VINT

      integer*4   MSTU(200),MSTJ(200)
      real*8      PARU(200),PARJ(200)
       COMMON/PYDAT1/MSTU,PARU,MSTJ,PARJ

      integer     KCHG(500,4)
      real*8      PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT2/KCHG,PMAS,PARF,VCKM
C ********************************************************************
