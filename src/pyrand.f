C*********************************************************************
C...PYRAND
C...Generates quantities characterizing the high-pT scattering at the
C...parton level according to the matrix elements. Chooses incoming,
C...reacting partons, their momentum fractions and one of the possible
C...subprocesses.
      SUBROUTINE PYRAND
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER PYK,PYCHGE,PYCOMP
C...Parameter statement to help give large particle numbers.
      PARAMETER (KSUSY1=1000000,KSUSY2=2000000,KEXCIT=4000000)
      logical           ROT_kinem_flag, ROT_pyrand_flag
       common /GEP_ROT/ ROT_kinem_flag, ROT_pyrand_flag
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(4000,2),BRAT(4000),KFDP(4000,5)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(500),KFPR(500,2),COEF(500,20),ICOL(40,4,2)
      COMMON/PYINT3/XSFX(2,-40:40),ISIG(1000,3),SIGH(1000)
      COMMON/PYINT4/MWID(500),WIDS(500,5)
      COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
      COMMON/PYINT7/SIGT(0:6,0:6,0:5)
      COMMON/PYUPPR/NUP,KUP(20,7),NFUP,IFUP(10,2),PUP(20,5),Q2UP(0:10)
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)
      SAVE /PYDAT1/,/PYDAT2/,/PYDAT3/,/PYSUBS/,/PYPARS/,/PYINT1/,
     &/PYINT2/,/PYINT3/,/PYINT4/,/PYINT5/,/PYINT7/,/PYUPPR/,/PYMSSM/
C...Local arrays.
      DIMENSION XPQ(-25:25),PMM(2),PDIF(4),BHAD(4),PMMN(2)
C...Parameters and data used in elastic/diffractive treatment.
      DATA EPS/0.0808D0/, ALP/0.25D0/, CRES/2D0/, PMRC/1.062D0/,
     &SMP/0.880D0/, BHAD/2.3D0,1.4D0,1.4D0,0.23D0/
C...Initial values, specifically for (first) semihard interaction.
      MINT(10)=0
      MINT(17)=0
      MINT(18)=0
      VINT(143)=1D0
      VINT(144)=1D0
      MFAIL=0
      IF(MSTP(171).EQ.1.AND.MSTP(172).EQ.2) MFAIL=1
      ISUB=0
      LOOP=0
  100 LOOP=LOOP+1
      MINT(51)=0
C...Choice of process type - first event of pileup.
      IF(MINT(82).EQ.1.AND.(ISUB.LE.90.OR.ISUB.GT.96)) THEN
C...For gamma-p or gamma-gamma first pick between alternatives.
        IF(MINT(121).GT.1) CALL PYSAVE(4,IGA)
        MINT(122)=IGA
C...For gamma + gamma with different nature, flip at random.
        IF(MINT(11).EQ.22.AND.MINT(12).EQ.22.AND.MINT(123).GE.4.AND.
     &  PYR(0).GT.0.5D0) THEN
          MINTSV=MINT(41)
          MINT(41)=MINT(42)
          MINT(42)=MINTSV
          MINTSV=MINT(45)
          MINT(45)=MINT(46)
          MINT(46)=MINTSV
          MINTSV=MINT(107)
          MINT(107)=MINT(108)
          MINT(108)=MINTSV
          IF(MINT(47).EQ.2.OR.MINT(47).EQ.3) MINT(47)=5-MINT(47)
        ENDIF
C...Pick process type.
        RSUB=XSEC(0,1)*PYR(0)
        DO 110 I=1,500
          IF(MSUB(I).NE.1) GOTO 110
          ISUB=I
          RSUB=RSUB-XSEC(I,1)
          IF(RSUB.LE.0D0) GOTO 120
  110   CONTINUE
  120   IF(ISUB.EQ.95) ISUB=96
        IF(ISUB.EQ.96) CALL PYMULT(2)
C...Choice of inclusive process type - pileup events.
      ELSEIF(MINT(82).GE.2.AND.ISUB.EQ.0) THEN
        RSUB=VINT(131)*PYR(0)
        ISUB=96
        IF(RSUB.GT.SIGT(0,0,5)) ISUB=94
        IF(RSUB.GT.SIGT(0,0,5)+SIGT(0,0,4)) ISUB=93
        IF(RSUB.GT.SIGT(0,0,5)+SIGT(0,0,4)+SIGT(0,0,3)) ISUB=92
        IF(RSUB.GT.SIGT(0,0,5)+SIGT(0,0,4)+SIGT(0,0,3)+SIGT(0,0,2))
     &  ISUB=91
        IF(ISUB.EQ.96) CALL PYMULT(2)
      ENDIF
      IF(MINT(82).EQ.1) NGEN(0,1)=NGEN(0,1)+1
      IF(MINT(82).EQ.1) NGEN(ISUB,1)=NGEN(ISUB,1)+1
      IF(ISUB.EQ.96.AND.LOOP.EQ.1.AND.MINT(82).EQ.1)
     &NGEN(97,1)=NGEN(97,1)+1
      MINT(1)=ISUB
      ISTSB=ISET(ISUB)
C...Random choice of flavour for some SUSY processes.
      IF(ISUB.GE.201.AND.ISUB.LE.280) THEN
C...~e_L ~nu_e or ~mu_L ~nu_mu.
        IF(ISUB.EQ.210) THEN
          KFPR(ISUB,1)=KSUSY1+11+2*INT(0.5D0+PYR(0))
          KFPR(ISUB,2)=KFPR(ISUB,1)+1
C...~nu_e ~nu_e(bar) or ~nu_mu ~nu_mu(bar).
        ELSEIF(ISUB.EQ.213) THEN
          KFPR(ISUB,1)=KSUSY1+12+2*INT(0.5D0+PYR(0))
          KFPR(ISUB,2)=KFPR(ISUB,1)
C...~q ~chi/~g; ~q = ~d, ~u, ~s, ~c or ~b.
        ELSEIF(ISUB.GE.246.AND.ISUB.LE.259) THEN
          IF(MOD(ISUB,2).EQ.0) THEN
            KFPR(ISUB,1)=KSUSY1+1+INT(5D0*PYR(0))
          ELSE
            KFPR(ISUB,1)=KSUSY2+1+INT(5D0*PYR(0))
          ENDIF
C...~q1 ~q2; ~q = ~d, ~u, ~s, ~c or ~b.
        ELSEIF(ISUB.GE.271.AND.ISUB.LE.276) THEN
          IF(ISUB.EQ.271.OR.ISUB.EQ.274) THEN
            KSU1=KSUSY1
            KSU2=KSUSY1
          ELSEIF(ISUB.EQ.272.OR.ISUB.EQ.275) THEN
            KSU1=KSUSY2
            KSU2=KSUSY2
          ELSEIF(PYR(0).LT.0.5D0) THEN
            KSU1=KSUSY1
            KSU2=KSUSY2
          ELSE
            KSU1=KSUSY2
            KSU2=KSUSY1
          ENDIF
          KFPR(ISUB,1)=KSU1+1+INT(5D0*PYR(0))
          KFPR(ISUB,2)=KSU2+1+INT(5D0*PYR(0))
C...~q ~q(bar);  ~q = ~d, ~u, ~s, ~c or ~b.
        ELSEIF(ISUB.EQ.277.OR.ISUB.EQ.279) THEN
          KFPR(ISUB,1)=KSUSY1+1+INT(5D0*PYR(0))
          KFPR(ISUB,2)=KFPR(ISUB,1)
        ELSEIF(ISUB.EQ.278.OR.ISUB.EQ.280) THEN
          KFPR(ISUB,1)=KSUSY2+1+INT(5D0*PYR(0))
          KFPR(ISUB,2)=KFPR(ISUB,1)
        ENDIF
      ENDIF
C...Find resonances (explicit or implicit in cross-section).
      MINT(72)=0
      KFR1=0
      IF(ISTSB.EQ.1.OR.ISTSB.EQ.3.OR.ISTSB.EQ.5) THEN
        KFR1=KFPR(ISUB,1)
      ELSEIF(ISUB.EQ.24.OR.ISUB.EQ.25.OR.ISUB.EQ.110.OR.ISUB.EQ.165.OR.
     &  ISUB.EQ.171.OR.ISUB.EQ.176) THEN
        KFR1=23
      ELSEIF(ISUB.EQ.23.OR.ISUB.EQ.26.OR.ISUB.EQ.166.OR.ISUB.EQ.172.OR.
     &  ISUB.EQ.177) THEN
        KFR1=24
      ELSEIF(ISUB.GE.71.AND.ISUB.LE.77) THEN
        KFR1=25
        IF(MSTP(46).EQ.5) THEN
          KFR1=30
          PMAS(30,1)=PARP(45)
          PMAS(30,2)=PARP(45)**3/(96D0*PARU(1)*PARP(47)**2)
        ENDIF
      ELSEIF(ISUB.EQ.194) THEN
        KFR1=54
      ENDIF
      CKMX=CKIN(2)
      IF(CKMX.LE.0D0) CKMX=VINT(1)
      KCR1=PYCOMP(KFR1)
      IF(KFR1.NE.0) THEN
        IF(CKIN(1).GT.PMAS(KCR1,1)+20D0*PMAS(KCR1,2).OR.
     &  CKMX.LT.PMAS(KCR1,1)-20D0*PMAS(KCR1,2)) KFR1=0
      ENDIF
      IF(KFR1.NE.0) THEN
        TAUR1=PMAS(KCR1,1)**2/VINT(2)
        GAMR1=PMAS(KCR1,1)*PMAS(KCR1,2)/VINT(2)
        MINT(72)=1
        MINT(73)=KFR1
        VINT(73)=TAUR1
        VINT(74)=GAMR1
      ENDIF
      IF(ISUB.EQ.141.OR.ISUB.EQ.194) THEN
        KFR2=23
        IF(ISUB.EQ.194) KFR2=56
        KCR2=PYCOMP(KFR2)
        TAUR2=PMAS(KCR2,1)**2/VINT(2)
        GAMR2=PMAS(KCR2,1)*PMAS(KCR2,2)/VINT(2)
        IF(CKIN(1).GT.PMAS(KCR2,1)+20D0*PMAS(KCR2,2).OR.
     &  CKMX.LT.PMAS(KCR2,1)-20D0*PMAS(KCR2,2)) KFR2=0
        IF(KFR2.NE.0.AND.KFR1.NE.0) THEN
          MINT(72)=2
          MINT(74)=KFR2
          VINT(75)=TAUR2
          VINT(76)=GAMR2
        ELSEIF(KFR2.NE.0) THEN
          KFR1=KFR2
          TAUR1=TAUR2
          GAMR1=GAMR2
          MINT(72)=1
          MINT(73)=KFR1
          VINT(73)=TAUR1
          VINT(74)=GAMR1
        ENDIF
      ENDIF
C...Find product masses and minimum pT of process,
C...optionally with broadening according to a truncated Breit-Wigner.
      VINT(63)=0D0
      VINT(64)=0D0
      MINT(71)=0
      VINT(71)=CKIN(3)
      IF(MINT(82).GE.2) VINT(71)=0D0
      VINT(80)=1D0
      IF(ISTSB.EQ.2.OR.ISTSB.EQ.4) THEN
        NBW=0
        DO 140 I=1,2
          PMMN(I)=0D0
          IF(KFPR(ISUB,I).EQ.0) THEN
          ELSEIF(MSTP(42).LE.0.OR.PMAS(PYCOMP(KFPR(ISUB,I)),2).LT.
     &      PARP(41)) THEN
            VINT(62+I)=PMAS(PYCOMP(KFPR(ISUB,I)),1)**2
          ELSE
            NBW=NBW+1
C...This prevents SUSY/t particles from becoming too light.
            KFLW=KFPR(ISUB,I)
            IF(KFLW/KSUSY1.EQ.1.OR.KFLW/KSUSY1.EQ.2) THEN
              KCW=PYCOMP(KFLW)
              PMMN(I)=PMAS(KCW,1)
              DO 130 IDC=MDCY(KCW,2),MDCY(KCW,2)+MDCY(KCW,3)-1
                IF(MDME(IDC,1).GT.0.AND.BRAT(IDC).GT.1E-4) THEN
                  PMSUM=PMAS(PYCOMP(KFDP(IDC,1)),1)+
     &            PMAS(PYCOMP(KFDP(IDC,2)),1)
                  IF(KFDP(IDC,3).NE.0) PMSUM=PMSUM+
     &            PMAS(PYCOMP(KFDP(IDC,3)),1)
                  PMMN(I)=MIN(PMMN(I),PMSUM)
                ENDIF
  130         CONTINUE
            ELSEIF(KFLW.EQ.6) THEN
              PMMN(I)=PMAS(24,1)+PMAS(5,1)
            ENDIF
          ENDIF
  140   CONTINUE
        IF(NBW.GE.1) THEN
          CKIN41=CKIN(41)
          CKIN43=CKIN(43)
          CKIN(41)=MAX(PMMN(1),CKIN(41))
          CKIN(43)=MAX(PMMN(2),CKIN(43))
          CALL PYOFSH(4,0,KFPR(ISUB,1),KFPR(ISUB,2),0D0,PQM3,PQM4)
          CKIN(41)=CKIN41
          CKIN(43)=CKIN43
          IF(MINT(51).EQ.1) THEN
            IF(MINT(121).GT.1) CALL PYSAVE(2,IGA)
            IF(MFAIL.EQ.1) THEN
              MSTI(61)=1
              RETURN
            ENDIF
            GOTO 100
          ENDIF
          VINT(63)=PQM3**2
          VINT(64)=PQM4**2
        ENDIF
        IF(MIN(VINT(63),VINT(64)).LT.CKIN(6)**2) MINT(71)=1
        IF(MINT(71).EQ.1) VINT(71)=MAX(CKIN(3),CKIN(5))
      ENDIF
C...Prepare for additional variable choices in 2 -> 3.
      IF(ISTSB.EQ.5) THEN
        VINT(201)=0D0
        IF(KFPR(ISUB,2).GT.0) VINT(201)=PMAS(PYCOMP(KFPR(ISUB,2)),1)
        VINT(206)=VINT(201)
        VINT(204)=PMAS(23,1)
        IF(ISUB.EQ.124) VINT(204)=PMAS(24,1)
        IF(ISUB.EQ.121.OR.ISUB.EQ.122.OR.ISUB.EQ.181.OR.ISUB.EQ.182.OR.
     &  ISUB.EQ.186.OR.ISUB.EQ.187) VINT(204)=VINT(201)
        VINT(209)=VINT(204)
      ENDIF
C...Select incoming VDM particle (rho/omega/phi/J/psi).
      IF(ISTSB.NE.0.AND.(MINT(101).GE.2.OR.MINT(102).GE.2).AND.
     &(MINT(123).EQ.2.OR.MINT(123).EQ.5.OR.MINT(123).EQ.7)) THEN
        VRN=PYR(0)*SIGT(0,0,5)
        IF(MINT(101).LE.1) THEN
          I1MN=0
          I1MX=0
        ELSE
          I1MN=1
          I1MX=MINT(101)
        ENDIF
        IF(MINT(102).LE.1) THEN
          I2MN=0
          I2MX=0
        ELSE
          I2MN=1
          I2MX=MINT(102)
        ENDIF
        DO 160 I1=I1MN,I1MX
          KFV1=110*I1+3
          DO 150 I2=I2MN,I2MX
            KFV2=110*I2+3
            VRN=VRN-SIGT(I1,I2,5)
            IF(VRN.LE.0D0) GOTO 170
  150     CONTINUE
  160   CONTINUE
  170   IF(MINT(101).GE.2) MINT(103)=KFV1
        IF(MINT(102).GE.2) MINT(104)=KFV2
      ENDIF
      IF(ISTSB.EQ.0) THEN
C...Elastic scattering or single or double diffractive scattering.
C...Select incoming particle (rho/omega/phi/J/psi for VDM) and mass.
        MINT(103)=MINT(11)
        MINT(104)=MINT(12)
        PMM(1)=VINT(3)
        PMM(2)=VINT(4)
        IF(MINT(101).GE.2.OR.MINT(102).GE.2) THEN
          JJ=ISUB-90
          VRN=PYR(0)*SIGT(0,0,JJ)
          IF(MINT(101).LE.1) THEN
            I1MN=0
            I1MX=0
          ELSE
            I1MN=1
            I1MX=MINT(101)
          ENDIF
          IF(MINT(102).LE.1) THEN
            I2MN=0
            I2MX=0
          ELSE
            I2MN=1
            I2MX=MINT(102)
          ENDIF
          DO 190 I1=I1MN,I1MX
            KFV1=110*I1+3
            DO 180 I2=I2MN,I2MX
              KFV2=110*I2+3
              VRN=VRN-SIGT(I1,I2,JJ)
              IF(VRN.LE.0D0) GOTO 200
  180       CONTINUE
  190     CONTINUE
  200     IF(MINT(101).GE.2) THEN
            MINT(103)=KFV1
            PMM(1)=PYMASS(KFV1)
          ENDIF
          IF(MINT(102).GE.2) THEN
            MINT(104)=KFV2
            PMM(2)=PYMASS(KFV2)
          ENDIF
        ENDIF
C...Side/sides of diffractive system.
        MINT(17)=0
        MINT(18)=0
        IF(ISUB.EQ.92.OR.ISUB.EQ.94) MINT(17)=1
        IF(ISUB.EQ.93.OR.ISUB.EQ.94) MINT(18)=1
C...Find masses of particles and minimal masses of diffractive states.
        DO 210 JT=1,2
          PDIF(JT)=PMM(JT)
          VINT(66+JT)=PDIF(JT)
          IF(MINT(16+JT).EQ.1) PDIF(JT)=PDIF(JT)+PARP(102)
  210   CONTINUE
        SH=VINT(2)
        SQM1=PMM(1)**2
        SQM2=PMM(2)**2
        SQM3=PDIF(1)**2
        SQM4=PDIF(2)**2
        SMRES1=(PMM(1)+PMRC)**2
        SMRES2=(PMM(2)+PMRC)**2
C...Find elastic slope and lower limit diffractive slope.
        IHA=MAX(2,IABS(MINT(103))/110)
        IF(IHA.GE.5) IHA=1
        IHB=MAX(2,IABS(MINT(104))/110)
        IF(IHB.GE.5) IHB=1
        IF(ISUB.EQ.91) THEN
          BMN=2D0*BHAD(IHA)+2D0*BHAD(IHB)+4D0*SH**EPS-4.2D0
        ELSEIF(ISUB.EQ.92) THEN
          BMN=MAX(2D0,2D0*BHAD(IHB))
        ELSEIF(ISUB.EQ.93) THEN
          BMN=MAX(2D0,2D0*BHAD(IHA))
        ELSEIF(ISUB.EQ.94) THEN
          BMN=2D0*ALP*4D0
        ENDIF
C...Determine maximum possible t range and coefficient of generation.
        SQLA12=(SH-SQM1-SQM2)**2-4D0*SQM1*SQM2
        SQLA34=(SH-SQM3-SQM4)**2-4D0*SQM3*SQM4
        THA=SH-(SQM1+SQM2+SQM3+SQM4)+(SQM1-SQM2)*(SQM3-SQM4)/SH
        THB=SQRT(MAX(0D0,SQLA12))*SQRT(MAX(0D0,SQLA34))/SH
        THC=(SQM3-SQM1)*(SQM4-SQM2)+(SQM1+SQM4-SQM2-SQM3)*
     &  (SQM1*SQM4-SQM2*SQM3)/SH
        THL=-0.5D0*(THA+THB)
        THU=THC/THL
        THRND=EXP(MAX(-50D0,BMN*(THL-THU)))-1D0
C...Select diffractive mass/masses according to dm^2/m^2.
  220   DO 230 JT=1,2
          IF(MINT(16+JT).EQ.0) THEN
            PDIF(2+JT)=PDIF(JT)
          ELSE
            PMMIN=PDIF(JT)
            PMMAX=MAX(VINT(2+JT),VINT(1)-PDIF(3-JT))
            PDIF(2+JT)=PMMIN*(PMMAX/PMMIN)**PYR(0)
          ENDIF
  230   CONTINUE
        SQM3=PDIF(3)**2
        SQM4=PDIF(4)**2
C..Additional mass factors, including resonance enhancement.
        IF(PDIF(3)+PDIF(4).GE.VINT(1)) GOTO 220
        IF(ISUB.EQ.92) THEN
          FSD=(1D0-SQM3/SH)*(1D0+CRES*SMRES1/(SMRES1+SQM3))
          IF(FSD.LT.PYR(0)*(1D0+CRES)) GOTO 220
        ELSEIF(ISUB.EQ.93) THEN
          FSD=(1D0-SQM4/SH)*(1D0+CRES*SMRES2/(SMRES2+SQM4))
          IF(FSD.LT.PYR(0)*(1D0+CRES)) GOTO 220
        ELSEIF(ISUB.EQ.94) THEN
          FDD=(1D0-(PDIF(3)+PDIF(4))**2/SH)*(SH*SMP/
     &    (SH*SMP+SQM3*SQM4))*(1D0+CRES*SMRES1/(SMRES1+SQM3))*
     &    (1D0+CRES*SMRES2/(SMRES2+SQM4))
          IF(FDD.LT.PYR(0)*(1D0+CRES)**2) GOTO 220
        ENDIF
C...Select t according to exp(Bmn*t) and correct to right slope.
        TH=THU+LOG(1D0+THRND*PYR(0))/BMN
        IF(ISUB.GE.92) THEN
          IF(ISUB.EQ.92) THEN
            BADD=2D0*ALP*LOG(SH/SQM3)
            IF(BHAD(IHB).LT.1D0) BADD=MAX(0D0,BADD+2D0*BHAD(IHB)-2D0)
          ELSEIF(ISUB.EQ.93) THEN
            BADD=2D0*ALP*LOG(SH/SQM4)
            IF(BHAD(IHA).LT.1D0) BADD=MAX(0D0,BADD+2D0*BHAD(IHA)-2D0)
          ELSEIF(ISUB.EQ.94) THEN
            BADD=2D0*ALP*(LOG(EXP(4D0)+SH/(ALP*SQM3*SQM4))-4D0)
          ENDIF
          IF(EXP(MAX(-50D0,BADD*(TH-THU))).LT.PYR(0)) GOTO 220
        ENDIF
C...Check whether m^2 and t choices are consistent.
        SQLA34=(SH-SQM3-SQM4)**2-4D0*SQM3*SQM4
        THA=SH-(SQM1+SQM2+SQM3+SQM4)+(SQM1-SQM2)*(SQM3-SQM4)/SH
        THB=SQRT(MAX(0D0,SQLA12))*SQRT(MAX(0D0,SQLA34))/SH
        IF(THB.LE.1D-8) GOTO 220
        THC=(SQM3-SQM1)*(SQM4-SQM2)+(SQM1+SQM4-SQM2-SQM3)*
     &  (SQM1*SQM4-SQM2*SQM3)/SH
        THLM=-0.5D0*(THA+THB)
        THUM=THC/THLM
        IF(TH.LT.THLM.OR.TH.GT.THUM) GOTO 220
C...Information to output.
        VINT(21)=1D0
        VINT(22)=0D0
        VINT(23)=MIN(1D0,MAX(-1D0,(THA+2D0*TH)/THB))
        VINT(45)=TH
        VINT(59)=2D0*SQRT(MAX(0D0,-(THC+THA*TH+TH**2)))/THB
        VINT(63)=PDIF(3)**2
        VINT(64)=PDIF(4)**2
C...Note: in the following, by In is meant the integral over the
C...quantity multiplying coefficient cn.
C...Choose tau according to h1(tau)/tau, where
C...h1(tau) = c1 + I1/I2*c2*1/tau + I1/I3*c3*1/(tau+tau_R) +
C...I1/I4*c4*tau/((s*tau-m^2)^2+(m*Gamma)^2) +
C...I1/I5*c5*1/(tau+tau_R') +
C...I1/I6*c6*tau/((s*tau-m'^2)^2+(m'*Gamma')^2) +
C...I1/I7*c7*tau/(1.-tau), and
C...c1 + c2 + c3 + c4 + c5 + c6 + c7 = 1.
      ELSEIF(ISTSB.GE.1.AND.ISTSB.LE.5) THEN
        CALL PYKLIM(1)
        IF(MINT(51).NE.0) THEN
          IF(MINT(121).GT.1) CALL PYSAVE(2,IGA)
          IF(MFAIL.EQ.1) THEN
            MSTI(61)=1
            RETURN
          ENDIF
          GOTO 100
        ENDIF
        RTAU=PYR(0)
        MTAU=1
        IF(RTAU.GT.COEF(ISUB,1)) MTAU=2
        IF(RTAU.GT.COEF(ISUB,1)+COEF(ISUB,2)) MTAU=3
        IF(RTAU.GT.COEF(ISUB,1)+COEF(ISUB,2)+COEF(ISUB,3)) MTAU=4
        IF(RTAU.GT.COEF(ISUB,1)+COEF(ISUB,2)+COEF(ISUB,3)+COEF(ISUB,4))
     &  MTAU=5
        IF(RTAU.GT.COEF(ISUB,1)+COEF(ISUB,2)+COEF(ISUB,3)+COEF(ISUB,4)+
     &  COEF(ISUB,5)) MTAU=6
        IF(RTAU.GT.COEF(ISUB,1)+COEF(ISUB,2)+COEF(ISUB,3)+COEF(ISUB,4)+
     &  COEF(ISUB,5)+COEF(ISUB,6)) MTAU=7
        CALL PYKMAP(1,MTAU,PYR(0))
C...2 -> 3, 4 processes:
C...Choose tau' according to h4(tau,tau')/tau', where
C...h4(tau,tau') = c1 + I1/I2*c2*(1 - tau/tau')^3/tau' +
C...I1/I3*c3*1/(1 - tau'), and c1 + c2 + c3 = 1.
        IF(ISTSB.GE.3.AND.ISTSB.LE.5) THEN
          CALL PYKLIM(4)
          IF(MINT(51).NE.0) THEN
            IF(MINT(121).GT.1) CALL PYSAVE(2,IGA)
            IF(MFAIL.EQ.1) THEN
              MSTI(61)=1
              RETURN
            ENDIF
            GOTO 100
          ENDIF
          RTAUP=PYR(0)
          MTAUP=1
          IF(RTAUP.GT.COEF(ISUB,18)) MTAUP=2
          IF(RTAUP.GT.COEF(ISUB,18)+COEF(ISUB,19)) MTAUP=3
          CALL PYKMAP(4,MTAUP,PYR(0))
        ENDIF
C...Choose y* according to h2(y*), where
C...h2(y*) = I0/I1*c1*(y*-y*min) + I0/I2*c2*(y*max-y*) +
C...I0/I3*c3*1/cosh(y*) + I0/I4*c4*1/(1-exp(y*-y*max)) +
C...I0/I5*c5*1/(1-exp(-y*-y*min)), I0 = y*max-y*min,
C...and c1 + c2 + c3 + c4 + c5 = 1.
        CALL PYKLIM(2)
        IF(MINT(51).NE.0) THEN
          IF(MINT(121).GT.1) CALL PYSAVE(2,IGA)
          IF(MFAIL.EQ.1) THEN
            MSTI(61)=1
            RETURN
          ENDIF
          GOTO 100
        ENDIF
        RYST=PYR(0)
        MYST=1
        IF(RYST.GT.COEF(ISUB,8)) MYST=2
        IF(RYST.GT.COEF(ISUB,8)+COEF(ISUB,9)) MYST=3
        IF(RYST.GT.COEF(ISUB,8)+COEF(ISUB,9)+COEF(ISUB,10)) MYST=4
        IF(RYST.GT.COEF(ISUB,8)+COEF(ISUB,9)+COEF(ISUB,10)+
     &  COEF(ISUB,11)) MYST=5
        CALL PYKMAP(2,MYST,PYR(0))
C...2 -> 2 processes:
C...Choose cos(theta-hat) (cth) according to h3(cth), where
C...h3(cth) = c0 + I0/I1*c1*1/(A - cth) + I0/I2*c2*1/(A + cth) +
C...I0/I3*c3*1/(A - cth)^2 + I0/I4*c4*1/(A + cth)^2,
C...A = 1 + 2*(m3*m4/sh)^2 (= 1 for massless products),
C...and c0 + c1 + c2 + c3 + c4 = 1.
        CALL PYKLIM(3)
        IF(MINT(51).NE.0) THEN
          IF(MINT(121).GT.1) CALL PYSAVE(2,IGA)
          IF(MFAIL.EQ.1) THEN
            MSTI(61)=1
            RETURN
          ENDIF
          GOTO 100
        ENDIF
        IF(ISTSB.EQ.2.OR.ISTSB.EQ.4) THEN
          RCTH=PYR(0)
          MCTH=1
          IF(RCTH.GT.COEF(ISUB,13)) MCTH=2
          IF(RCTH.GT.COEF(ISUB,13)+COEF(ISUB,14)) MCTH=3
          IF(RCTH.GT.COEF(ISUB,13)+COEF(ISUB,14)+COEF(ISUB,15)) MCTH=4
          IF(RCTH.GT.COEF(ISUB,13)+COEF(ISUB,14)+COEF(ISUB,15)+
     &    COEF(ISUB,16)) MCTH=5
          CALL PYKMAP(3,MCTH,PYR(0))
        ENDIF
C...2 -> 3 : select pT1, phi1, pT2, phi2, y3 for 3 outgoing.
        IF(ISTSB.EQ.5) THEN
          CALL PYKMAP(5,0,0D0)
          IF(MINT(51).NE.0) THEN
            IF(MINT(121).GT.1) CALL PYSAVE(2,IGA)
            IF(MFAIL.EQ.1) THEN
              MSTI(61)=1
              RETURN
            ENDIF
            GOTO 100
          ENDIF
        ENDIF
C...Low-pT or multiple interactions (first semihard interaction).
      ELSEIF(ISTSB.EQ.9) THEN
        CALL PYMULT(3)
        ISUB=MINT(1)
C...Generate user-defined process: kinematics plus weight.
      ELSEIF(ISTSB.EQ.11) THEN
        MSTI(51)=0
        CALL PYUPEV(ISUB,SIGS)
        IF(NUP.LE.0) THEN
          MINT(51)=2
          MSTI(51)=1
          IF(MINT(82).EQ.1) THEN
            NGEN(0,1)=NGEN(0,1)-1
            NGEN(0,2)=NGEN(0,2)-1
            NGEN(ISUB,1)=NGEN(ISUB,1)-1
          ENDIF
          IF(MINT(121).GT.1) CALL PYSAVE(2,IGA)
          RETURN
        ENDIF
C...Construct 'trivial' kinematical variables needed.
        KFL1=KUP(1,2)
        KFL2=KUP(2,2)
        VINT(41)=2D0*PUP(1,4)/VINT(1)
        VINT(42)=2D0*PUP(2,4)/VINT(1)
        VINT(21)=VINT(41)*VINT(42)
        VINT(22)=0.5D0*LOG(VINT(41)/VINT(42))
        VINT(44)=VINT(21)*VINT(2)
        VINT(43)=SQRT(MAX(0D0,VINT(44)))
        VINT(56)=Q2UP(0)
        VINT(55)=SQRT(MAX(0D0,VINT(56)))
C...Construct other kinematical variables needed (approximately).
        VINT(23)=0D0
        VINT(26)=VINT(21)
        VINT(45)=-0.5D0*VINT(44)
        VINT(46)=-0.5D0*VINT(44)
        VINT(49)=VINT(43)
        VINT(50)=VINT(44)
        VINT(51)=VINT(55)
        VINT(52)=VINT(56)
        VINT(53)=VINT(55)
        VINT(54)=VINT(56)
        VINT(25)=0D0
        VINT(48)=0D0
        DO 240 IUP=3,NUP
          IF(KUP(IUP,1).EQ.1) VINT(25)=VINT(25)+2D0*(PUP(IUP,5)**2+
     &    PUP(IUP,1)**2+PUP(IUP,2)**2)/VINT(1)
          IF(KUP(IUP,1).EQ.1) VINT(48)=VINT(48)+0.5D0*(PUP(IUP,1)**2+
     &    PUP(IUP,2)**2)
  240   CONTINUE
        VINT(47)=SQRT(VINT(48))
C...Calculate parton distribution weights.
        IF(MINT(47).GE.2) THEN
          DO 260 I=3-MIN(2,MINT(45)),MIN(2,MINT(46))
            MINT(105)=MINT(102+I)
            MINT(109)=MINT(106+I)
            IF(MSTP(57).LE.1) THEN
              CALL PYPDFU(MINT(10+I),VINT(40+I),Q2UP(0),XPQ)
            ELSE
              CALL PYPDFL(MINT(10+I),VINT(40+I),Q2UP(0),XPQ)
            ENDIF
            DO 250 KFL=-25,25
              XSFX(I,KFL)=XPQ(KFL)
  250       CONTINUE
  260     CONTINUE
        ENDIF
      ENDIF
C...Choose azimuthal angle.
      if (ROT_pyrand_flag) then      
         VINT(24)= PARU(2)*PYR(0)
       else
         VINT(24)= 0.D0
      endif
C...Check against user cuts on kinematics at parton level.
      MINT(51)=0
      IF((ISUB.LE.90.OR.ISUB.GT.100).AND.ISTSB.LE.10) CALL PYKLIM(0)
      IF(MINT(51).NE.0) THEN
        IF(MINT(121).GT.1) CALL PYSAVE(2,IGA)
        IF(MFAIL.EQ.1) THEN
          MSTI(61)=1
          RETURN
        ENDIF
        GOTO 100
      ENDIF
      IF(MINT(82).EQ.1.AND.MSTP(141).GE.1.AND.ISTSB.LE.10) THEN
        MCUT=0
        IF(MSUB(91)+MSUB(92)+MSUB(93)+MSUB(94)+MSUB(95).EQ.0)
     &  CALL PYKCUT(MCUT)
        IF(MCUT.NE.0) THEN
          IF(MINT(121).GT.1) CALL PYSAVE(2,IGA)
          IF(MFAIL.EQ.1) THEN
            MSTI(61)=1
            RETURN
          ENDIF
          GOTO 100
        ENDIF
      ENDIF
C...Calculate differential cross-section for different subprocesses.
      IF(ISTSB.LE.10) CALL PYSIGH(NCHN,SIGS)
      SIGSOR=SIGS
      SIGLPT=SIGT(0,0,5)
C...Multiply cross-section by user-defined weights.
      IF(MSTP(173).EQ.1) THEN
        SIGS=PARP(173)*SIGS
        DO 270 ICHN=1,NCHN
          SIGH(ICHN)=PARP(173)*SIGH(ICHN)
  270   CONTINUE
        SIGLPT=PARP(173)*SIGLPT
      ENDIF
      WTXS=1D0
      SIGSWT=SIGS
      VINT(99)=1D0
      VINT(100)=1D0
      IF(MINT(82).EQ.1.AND.MSTP(142).GE.1) THEN
        IF(ISUB.NE.96.AND.MSUB(91)+MSUB(92)+MSUB(93)+MSUB(94)+
     &  MSUB(95).EQ.0) CALL PYEVWT(WTXS)
        SIGSWT=WTXS*SIGS
        VINT(99)=WTXS
        IF(MSTP(142).EQ.1) VINT(100)=1D0/WTXS
      ENDIF
C...Calculations for Monte Carlo estimate of all cross-sections.
      IF(MINT(82).EQ.1.AND.ISUB.LE.90.OR.ISUB.GE.96) THEN
        IF(MSTP(142).LE.1) THEN
          XSEC(ISUB,2)=XSEC(ISUB,2)+SIGS
        ELSE
          XSEC(ISUB,2)=XSEC(ISUB,2)+SIGSWT
        ENDIF
      ELSEIF(MINT(82).EQ.1) THEN
        XSEC(ISUB,2)=XSEC(ISUB,2)+SIGS
      ENDIF
      IF((ISUB.EQ.95.OR.ISUB.EQ.96).AND.LOOP.EQ.1.AND.MINT(82).EQ.1)
     &XSEC(97,2)=XSEC(97,2)+SIGLPT
C...Multiple interactions: store results of cross-section calculation.
      IF(MINT(50).EQ.1.AND.MSTP(82).GE.3) THEN
        VINT(153)=SIGSOR
        CALL PYMULT(4)
      ENDIF
C...Check that weight not negative.
      VIOL=SIGSWT/XSEC(ISUB,1)
      IF(ISUB.EQ.96.AND.MSTP(173).EQ.1) VIOL=VIOL/PARP(174)
      IF(MSTP(123).LE.0) THEN
        IF(VIOL.LT.-1D-3) THEN
          WRITE(MSTU(11),5000) VIOL,NGEN(0,3)+1
          IF(MSTP(122).GE.1) WRITE(MSTU(11),5100) ISUB,VINT(21),
     &    VINT(22),VINT(23),VINT(26)
          STOP
        ENDIF
      ELSE
        IF(VIOL.LT.MIN(-1D-3,VINT(109))) THEN
          VINT(109)=VIOL
          WRITE(MSTU(11),5200) VIOL,NGEN(0,3)+1
          IF(MSTP(122).GE.1) WRITE(MSTU(11),5100) ISUB,VINT(21),
     &    VINT(22),VINT(23),VINT(26)
        ENDIF
      ENDIF
C...Weighting using estimate of maximum of differential cross-section.
      IF(MFAIL.EQ.0) THEN
        IF(VIOL.LT.PYR(0)) THEN
          IF(MINT(121).GT.1) CALL PYSAVE(2,IGA)
          GOTO 100
        ENDIF
      ELSEIF(ISUB.NE.95.AND.ISUB.NE.96) THEN
        IF(VIOL.LT.PYR(0)) THEN
          MSTI(61)=1
          IF(MINT(121).GT.1) CALL PYSAVE(2,IGA)
          RETURN
        ENDIF
      ELSE
        RATND=SIGLPT/XSEC(95,1)
        IF(LOOP.EQ.1.AND.RATND.LT.PYR(0)) THEN
          MSTI(61)=1
          IF(MINT(121).GT.1) CALL PYSAVE(2,IGA)
          RETURN
        ENDIF
        VIOL=VIOL/RATND
        IF(VIOL.LT.PYR(0)) THEN
          IF(MINT(121).GT.1) CALL PYSAVE(2,IGA)
          GOTO 100
        ENDIF
      ENDIF
C...Check for possible violation of estimated maximum of differential
C...cross-section used in weighting.
      IF(MSTP(123).LE.0) THEN
        IF(VIOL.GT.1D0) THEN
          WRITE(MSTU(11),5300) VIOL,NGEN(0,3)+1
          IF(MSTP(122).GE.2) WRITE(MSTU(11),5100) ISUB,VINT(21),
     &    VINT(22),VINT(23),VINT(26)
          STOP
        ENDIF
      ELSEIF(MSTP(123).EQ.1) THEN
        IF(VIOL.GT.VINT(108)) THEN
          VINT(108)=VIOL
          IF(VIOL.GT.1D0) THEN
            MINT(10)=1
            WRITE(MSTU(11),5400) VIOL,NGEN(0,3)+1
            IF(MSTP(122).GE.2) WRITE(MSTU(11),5100) ISUB,VINT(21),
     &      VINT(22),VINT(23),VINT(26)
          ENDIF
        ENDIF
      ELSEIF(VIOL.GT.VINT(108)) THEN
        VINT(108)=VIOL
        IF(VIOL.GT.1D0) THEN
          MINT(10)=1
          XDIF=XSEC(ISUB,1)*(VIOL-1D0)
          XSEC(ISUB,1)=XSEC(ISUB,1)+XDIF
          IF(MSUB(ISUB).EQ.1.AND.(ISUB.LE.90.OR.ISUB.GT.96))
     &    XSEC(0,1)=XSEC(0,1)+XDIF
          WRITE(MSTU(11),5400) VIOL,NGEN(0,3)+1
          IF(MSTP(122).GE.2) WRITE(MSTU(11),5100) ISUB,VINT(21),
     &    VINT(22),VINT(23),VINT(26)
          IF(ISUB.LE.9) THEN
            WRITE(MSTU(11),5500) ISUB,XSEC(ISUB,1)
          ELSEIF(ISUB.LE.99) THEN
            WRITE(MSTU(11),5600) ISUB,XSEC(ISUB,1)
          ELSE
            WRITE(MSTU(11),5700) ISUB,XSEC(ISUB,1)
          ENDIF
          VINT(108)=1D0
        ENDIF
      ENDIF
C...Multiple interactions: choose impact parameter.
      VINT(148)=1D0
      IF(MINT(50).EQ.1.AND.(ISUB.LE.90.OR.ISUB.GE.96).AND.
     &MSTP(82).GE.3) THEN
        CALL PYMULT(5)
        IF(VINT(150).LT.PYR(0)) THEN
          IF(MINT(121).GT.1) CALL PYSAVE(2,IGA)
          IF(MFAIL.EQ.1) THEN
            MSTI(61)=1
            RETURN
          ENDIF
          GOTO 100
        ENDIF
      ENDIF
      IF(MINT(82).EQ.1) NGEN(0,2)=NGEN(0,2)+1
      IF(MINT(82).EQ.1.AND.MSUB(95).EQ.1) THEN
        IF(ISUB.LE.90.OR.ISUB.GE.95) NGEN(95,1)=NGEN(95,1)+1
        IF(ISUB.LE.90.OR.ISUB.GE.96) NGEN(96,2)=NGEN(96,2)+1
      ENDIF
      IF(ISUB.LE.90.OR.ISUB.GE.96) MINT(31)=MINT(31)+1
C...Choose flavour of reacting partons (and subprocess).
      IF(ISTSB.GE.11) GOTO 290
      RSIGS=SIGS*PYR(0)
      QT2=VINT(48)
      RQQBAR=PARP(87)*(1D0-(QT2/(QT2+(PARP(88)*PARP(82))**2))**2)
      IF(ISUB.NE.95.AND.(ISUB.NE.96.OR.MSTP(82).LE.1.OR.
     &PYR(0).GT.RQQBAR)) THEN
        DO 280 ICHN=1,NCHN
          KFL1=ISIG(ICHN,1)
          KFL2=ISIG(ICHN,2)
          MINT(2)=ISIG(ICHN,3)
          RSIGS=RSIGS-SIGH(ICHN)
          IF(RSIGS.LE.0D0) GOTO 290
  280   CONTINUE
C...Multiple interactions: choose qqbar preferentially at small pT.
      ELSEIF(ISUB.EQ.96) THEN
        MINT(105)=MINT(103)
        MINT(109)=MINT(107)
        CALL PYSPLI(MINT(11),21,KFL1,KFLDUM)
        MINT(105)=MINT(104)
        MINT(109)=MINT(108)
        CALL PYSPLI(MINT(12),21,KFL2,KFLDUM)
        MINT(1)=11
        MINT(2)=1
        IF(KFL1.EQ.KFL2.AND.PYR(0).LT.0.5D0) MINT(2)=2
C...Low-pT: choose string drawing configuration.
      ELSE
        KFL1=21
        KFL2=21
        RSIGS=6D0*PYR(0)
        MINT(2)=1
        IF(RSIGS.GT.1D0) MINT(2)=2
        IF(RSIGS.GT.2D0) MINT(2)=3
      ENDIF
C...Reassign QCD process. Partons before initial state radiation.
  290 IF(MINT(2).GT.10) THEN
        MINT(1)=MINT(2)/10
        MINT(2)=MOD(MINT(2),10)
      ENDIF
      IF(MINT(82).EQ.1.AND.MSTP(111).GE.0) NGEN(MINT(1),2)=
     &NGEN(MINT(1),2)+1
      MINT(15)=KFL1
      MINT(16)=KFL2
      MINT(13)=MINT(15)
      MINT(14)=MINT(16)
      VINT(141)=VINT(41)
      VINT(142)=VINT(42)
      VINT(151)=0D0
      VINT(152)=0D0
C...Calculate x value of photon for parton inside photon inside e.
      DO 320 JT=1,2
        MINT(18+JT)=0
        VINT(154+JT)=0D0
        MSPLI=0
        IF(JT.EQ.1.AND.MINT(43).LE.2) MSPLI=1
        IF(JT.EQ.2.AND.MOD(MINT(43),2).EQ.1) MSPLI=1
        IF(IABS(MINT(14+JT)).LE.8.OR.MINT(14+JT).EQ.21) MSPLI=MSPLI+1
        IF(MSPLI.EQ.2) THEN
          KFLH=MINT(14+JT)
          XHRD=VINT(140+JT)
          Q2HRD=VINT(54)
          MINT(105)=MINT(102+JT)
          MINT(109)=MINT(106+JT)
          IF(MSTP(57).LE.1) THEN
            CALL PYPDFU(22,XHRD,Q2HRD,XPQ)
          ELSE
            CALL PYPDFL(22,XHRD,Q2HRD,XPQ)
          ENDIF
          WTMX=4D0*XPQ(KFLH)
          IF(MSTP(13).EQ.2) THEN
            Q2PMS=Q2HRD/PMAS(11,1)**2
            WTMX=WTMX*LOG(MAX(2D0,Q2PMS*(1D0-XHRD)/XHRD**2))
          ENDIF
  300     XE=XHRD**PYR(0)
          XG=MIN(0.999999D0,XHRD/XE)
          IF(MSTP(57).LE.1) THEN
            CALL PYPDFU(22,XG,Q2HRD,XPQ)
          ELSE
            CALL PYPDFL(22,XG,Q2HRD,XPQ)
          ENDIF
          WT=(1D0+(1D0-XE)**2)*XPQ(KFLH)
          IF(MSTP(13).EQ.2) WT=WT*LOG(MAX(2D0,Q2PMS*(1D0-XE)/XE**2))
          IF(WT.LT.PYR(0)*WTMX) GOTO 300
          MINT(18+JT)=1
          VINT(154+JT)=XE
          DO 310 KFLS=-25,25
            XSFX(JT,KFLS)=XPQ(KFLS)
  310     CONTINUE
        ENDIF
  320 CONTINUE
C...Pick scale where photon is resolved.
      IF(MINT(107).EQ.3) VINT(283)=PARP(15)**2*
     &(VINT(54)/PARP(15)**2)**PYR(0)
      IF(MINT(108).EQ.3) VINT(284)=PARP(15)**2*
     &(VINT(54)/PARP(15)**2)**PYR(0)
      IF(MINT(121).GT.1) CALL PYSAVE(2,IGA)
C...Format statements for differential cross-section maximum violations.
 5000 FORMAT(/1X,'Error: negative cross-section fraction',1P,D11.3,1X,
     &'in event',1X,I7,'D0'/1X,'Execution stopped!')
 5100 FORMAT(1X,'ISUB = ',I3,'; Point of violation:'/1X,'tau =',1P,
     &D11.3,', y* =',D11.3,', cthe = ',0P,F11.7,', tau'' =',1P,D11.3)
 5200 FORMAT(/1X,'Warning: negative cross-section fraction',1P,D11.3,1X,
     &'in event',1X,I7)
 5300 FORMAT(/1X,'Error: maximum violated by',1P,D11.3,1X,
     &'in event',1X,I7,'D0'/1X,'Execution stopped!')
 5400 FORMAT(/1X,'Advisory warning: maximum violated by',1P,D11.3,1X,
     &'in event',1X,I7)
 5500 FORMAT(1X,'XSEC(',I1,',1) increased to',1P,D11.3)
 5600 FORMAT(1X,'XSEC(',I2,',1) increased to',1P,D11.3)
 5700 FORMAT(1X,'XSEC(',I3,',1) increased to',1P,D11.3)
      RETURN
      END
