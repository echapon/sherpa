C*********************************************************************
 
C...SPPDFU
C...Gives electron, muon, tau, photon, pi+, neutron, proton and hyperon
C...parton distributions according to a few different parametrizations.
C...Note that what is coded is x times the probability distribution,
C...i.e. xq(x,Q2) etc.
 
      SUBROUTINE SPPDFU(KF,X,Q2,XPQ)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER SPK,SPCHGE,SPCOMP
C...Commonblocks.
      COMMON/SPDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/SPDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/SPPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/SPINT1/MINT(400),VINT(400)
      COMMON/SPINT8/XPVMD(-6:6),XPANL(-6:6),XPANH(-6:6),XPBEH(-6:6),
     &XPDIR(-6:6)
      SAVE /SPDAT1/,/SPDAT2/,/SPPARS/,/SPINT1/,/SPINT8/
C...Local arrays.
      DIMENSION XPQ(-25:25),XPEL(-25:25),XPGA(-6:6),VXPGA(-6:6),
     &XPPI(-6:6),XPPR(-6:6)
 
C...Interface to PDFLIB.
      COMMON/W50513/XMIN,XMAX,Q2MIN,Q2MAX
      SAVE /W50513/
      DOUBLE PRECISION XX,QQ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU,
     &VALUE(20),XMIN,XMAX,Q2MIN,Q2MAX
      CHARACTER*20 PARM(20)
      DATA VALUE/20*0D0/,PARM/20*' '/
 
C...Data related to Schuler-Sjostrand photon distributions.
      DATA ALAMGA/0.2D0/, PMCGA/1.3D0/, PMBGA/4.6D0/
 
C...Reset parton distributions.
      MINT(92)=0
      DO 100 KFL=-25,25
        XPQ(KFL)=0D0
  100 CONTINUE
 
C...Check x and particle species.
      IF(X.LE.0D0.OR.X.GE.1D0) THEN
        WRITE(MSTU(11),5000) X
        RETURN
      ENDIF
      KFA=IABS(KF)
      IF(KFA.NE.11.AND.KFA.NE.13.AND.KFA.NE.15.AND.KFA.NE.22.AND.
     &KFA.NE.211.AND.KFA.NE.2112.AND.KFA.NE.2212.AND.KFA.NE.3122.AND.
     &KFA.NE.3112.AND.KFA.NE.3212.AND.KFA.NE.3222.AND.KFA.NE.3312.AND.
     &KFA.NE.3322.AND.KFA.NE.3334.AND.KFA.NE.111.AND.KFA.NE.321.AND.
     &KFA.NE.310.AND.KFA.NE.130) THEN
        WRITE(MSTU(11),5100) KF
        RETURN
      ENDIF
 
C...Electron (or muon or tau) parton distribution call.
      IF(KFA.EQ.11.OR.KFA.EQ.13.OR.KFA.EQ.15) THEN
        CALL SPPDEL(KFA,X,Q2,XPEL)
        DO 110 KFL=-25,25
          XPQ(KFL)=XPEL(KFL)
  110   CONTINUE
 
C...Photon parton distribution call (VDM+anomalous).
      ELSEIF(KFA.EQ.22.AND.MINT(109).LE.1) THEN
        IF(MSTP(56).EQ.1.AND.MSTP(55).EQ.1) THEN
          CALL SPPDGA(X,Q2,XPGA)
          DO 120 KFL=-6,6
            XPQ(KFL)=XPGA(KFL)
  120     CONTINUE
        ELSEIF(MSTP(56).EQ.1.AND.MSTP(55).GE.5.AND.MSTP(55).LE.8) THEN
          Q2MX=Q2
          P2MX=0.36D0
          IF(MSTP(55).GE.7) P2MX=4.0D0
          IF(MSTP(57).EQ.0) Q2MX=P2MX
          P2=0D0
          IF(VINT(120).LT.0D0) P2=VINT(120)**2
          CALL SPGGAM(MSTP(55)-4,X,Q2MX,P2,MSTP(60),F2GAM,XPGA)
          DO 130 KFL=-6,6
            XPQ(KFL)=XPGA(KFL)
  130     CONTINUE
          VINT(231)=P2MX
        ELSEIF(MSTP(56).EQ.1.AND.MSTP(55).GE.9.AND.MSTP(55).LE.12) THEN
          Q2MX=Q2
          P2MX=0.36D0
          IF(MSTP(55).GE.11) P2MX=4.0D0
          IF(MSTP(57).EQ.0) Q2MX=P2MX
          P2=0D0
          IF(VINT(120).LT.0D0) P2=VINT(120)**2
          CALL SPGGAM(MSTP(55)-8,X,Q2MX,P2,MSTP(60),F2GAM,XPGA)
          DO 140 KFL=-6,6
            XPQ(KFL)=XPVMD(KFL)+XPANL(KFL)+XPBEH(KFL)+XPDIR(KFL)
  140     CONTINUE
          VINT(231)=P2MX
        ELSEIF(MSTP(56).EQ.2) THEN
C...Call PDFLIB parton distributions.
          PARM(1)='NPTYPE'
          VALUE(1)=3
          PARM(2)='NGROUP'
          VALUE(2)=MSTP(55)/1000
          PARM(3)='NSET'
          VALUE(3)=MOD(MSTP(55),1000)
          IF(MINT(93).NE.3000000+MSTP(55)) THEN
            CALL SDFSET(PARM,VALUE)
            MINT(93)=3000000+MSTP(55)
          ENDIF
          XX=X
          QQ2=MAX(0D0,Q2MIN,Q2)
          IF(MSTP(57).EQ.0) QQ2=Q2MIN
          P2=0D0
          IF(VINT(120).LT.0D0) P2=VINT(120)**2
          IP2=MSTP(60)
          IF(MSTP(55).EQ.5004) THEN
            IF(5D0*P2.LT.QQ2.AND.
     &      QQ2.GT.0.6D0.AND.QQ2.LT.5D4.AND.
     &      P2.GE.0D0.AND.P2.LT.10D0.AND.
     &      XX.GT.1D-4.AND.XX.LT.1D0) THEN
              CALL PTRUCTP(XX,QQ2,P2,IP2,UPV,DNV,USEA,DSEA,STR,CHM,
     &        BOT,TOP,GLU)
            ELSE
              UPV=0D0
              DNV=0D0
              USEA=0D0
              DSEA=0D0
              STR=0D0
              CHM=0D0
              BOT=0D0
              TOP=0D0
              GLU=0D0
            ENDIF
          ELSE
            IF(P2.LT.QQ2) THEN
              CALL PTRUCTP(XX,QQ2,P2,IP2,UPV,DNV,USEA,DSEA,STR,CHM,
     &        BOT,TOP,GLU)
            ELSE
              UPV=0D0
              DNV=0D0
              USEA=0D0
              DSEA=0D0
              STR=0D0
              CHM=0D0
              BOT=0D0
              TOP=0D0
              GLU=0D0
            ENDIF
          ENDIF
          VINT(231)=Q2MIN
          XPQ(0)=GLU
          XPQ(1)=DNV
          XPQ(-1)=DNV
          XPQ(2)=UPV
          XPQ(-2)=UPV
          XPQ(3)=STR
          XPQ(-3)=STR
          XPQ(4)=CHM
          XPQ(-4)=CHM
          XPQ(5)=BOT
          XPQ(-5)=BOT
          XPQ(6)=TOP
          XPQ(-6)=TOP
        ELSE
          WRITE(MSTU(11),5200) KF,MSTP(56),MSTP(55)
        ENDIF
 
C...Pion/gammaVDM parton distribution call.
      ELSEIF(KFA.EQ.211.OR.KFA.EQ.111.OR.KFA.EQ.321.OR.KFA.EQ.130.OR.
     &KFA.EQ.310.OR.(KFA.EQ.22.AND.MINT(109).EQ.2)) THEN
        IF(KFA.EQ.22.AND.MSTP(56).EQ.1.AND.MSTP(55).GE.5.AND.
     &  MSTP(55).LE.12) THEN
          ISET=1+MOD(MSTP(55)-1,4)
          Q2MX=Q2
          P2MX=0.36D0
          IF(ISET.GE.3) P2MX=4.0D0
          IF(MSTP(57).EQ.0) Q2MX=P2MX
          P2=0D0
          IF(VINT(120).LT.0D0) P2=VINT(120)**2
          CALL SPGGAM(ISET,X,Q2MX,P2,MSTP(60),F2GAM,XPGA)
          DO 150 KFL=-6,6
            XPQ(KFL)=XPVMD(KFL)
  150     CONTINUE
          VINT(231)=P2MX
        ELSEIF(MSTP(54).EQ.1.AND.MSTP(53).GE.1.AND.MSTP(53).LE.3) THEN
          CALL SPPDPI(X,Q2,XPPI)
          DO 160 KFL=-6,6
            XPQ(KFL)=XPPI(KFL)
  160     CONTINUE
        ELSEIF(MSTP(54).EQ.2) THEN
C...Call PDFLIB parton distributions.
          PARM(1)='NPTYPE'
          VALUE(1)=2
          PARM(2)='NGROUP'
          VALUE(2)=MSTP(53)/1000
          PARM(3)='NSET'
          VALUE(3)=MOD(MSTP(53),1000)
          IF(MINT(93).NE.2000000+MSTP(53)) THEN
            CALL SDFSET(PARM,VALUE)
            MINT(93)=2000000+MSTP(53)
          ENDIF
          XX=X
          QQ=SQRT(MAX(0D0,Q2MIN,Q2))
          IF(MSTP(57).EQ.0) QQ=SQRT(Q2MIN)
          CALL PTRUCTM(XX,QQ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU)
          VINT(231)=Q2MIN
          XPQ(0)=GLU
          XPQ(1)=DSEA
          XPQ(-1)=UPV+DSEA
          XPQ(2)=UPV+USEA
          XPQ(-2)=USEA
          XPQ(3)=STR
          XPQ(-3)=STR
          XPQ(4)=CHM
          XPQ(-4)=CHM
          XPQ(5)=BOT
          XPQ(-5)=BOT
          XPQ(6)=TOP
          XPQ(-6)=TOP
        ELSE
          WRITE(MSTU(11),5200) KF,MSTP(54),MSTP(53)
        ENDIF
 
C...Anomalous photon parton distribution call.
      ELSEIF(KFA.EQ.22.AND.MINT(109).EQ.3) THEN
        Q2MX=Q2
        P2MX=PARP(15)**2
        IF(MSTP(56).EQ.1.AND.MSTP(55).LE.8) THEN
          IF(MSTP(55).EQ.5.OR.MSTP(55).EQ.6) P2MX=0.36D0
          IF(MSTP(55).EQ.7.OR.MSTP(55).EQ.8) P2MX=4.0D0
          IF(MSTP(57).EQ.0) Q2MX=P2MX
          P2=0D0
          IF(VINT(120).LT.0D0) P2=VINT(120)**2
          CALL SPGGAM(MSTP(55)-4,X,Q2MX,P2,MSTP(60),F2GM,XPGA)
          DO 170 KFL=-6,6
            XPQ(KFL)=XPANL(KFL)+XPANH(KFL)
  170     CONTINUE
          VINT(231)=P2MX
        ELSEIF(MSTP(56).EQ.1) THEN
          IF(MSTP(55).EQ.9.OR.MSTP(55).EQ.10) P2MX=0.36D0
          IF(MSTP(55).EQ.11.OR.MSTP(55).EQ.12) P2MX=4.0D0
          IF(MSTP(57).EQ.0) Q2MX=P2MX
          P2=0D0
          IF(VINT(120).LT.0D0) P2=VINT(120)**2
          CALL SPGGAM(MSTP(55)-8,X,Q2MX,P2,MSTP(60),F2GM,XPGA)
          DO 180 KFL=-6,6
            XPQ(KFL)=MAX(0D0,XPANL(KFL)+XPBEH(KFL)+XPDIR(KFL))
  180     CONTINUE
          VINT(231)=P2MX
        ELSEIF(MSTP(56).EQ.2) THEN
          IF(MSTP(57).EQ.0) Q2MX=P2MX
          CALL SPGANO(0,X,Q2MX,P2MX,ALAMGA,XPGA,VXPGA)
          DO 190 KFL=-6,6
            XPQ(KFL)=XPGA(KFL)
  190     CONTINUE
          VINT(231)=P2MX
        ELSEIF(MSTP(55).GE.1.AND.MSTP(55).LE.5) THEN
          IF(MSTP(57).EQ.0) Q2MX=P2MX
          CALL SPGVMD(0,MSTP(55),X,Q2MX,P2MX,PARP(1),XPGA,VXPGA)
          DO 200 KFL=-6,6
            XPQ(KFL)=XPGA(KFL)
  200     CONTINUE
          VINT(231)=P2MX
        ELSE
  210     RKF=11D0*SPR(0)
          KFR=1
          IF(RKF.GT.1D0) KFR=2
          IF(RKF.GT.5D0) KFR=3
          IF(RKF.GT.6D0) KFR=4
          IF(RKF.GT.10D0) KFR=5
          IF(KFR.EQ.4.AND.Q2.LT.PMCGA**2) GOTO 210
          IF(KFR.EQ.5.AND.Q2.LT.PMBGA**2) GOTO 210
          IF(MSTP(57).EQ.0) Q2MX=P2MX
          CALL SPGVMD(0,KFR,X,Q2MX,P2MX,PARP(1),XPGA,VXPGA)
          DO 220 KFL=-6,6
            XPQ(KFL)=XPGA(KFL)
  220     CONTINUE
          VINT(231)=P2MX
        ENDIF
 
C...Proton parton distribution call.
      ELSE
        IF(MSTP(52).EQ.1.AND.MSTP(51).GE.1.AND.MSTP(51).LE.20) THEN
C*sh* you may use Sherpa interfaces to PDF here
C          CALL SPPDPR(X,Q2,XPPR)
          CALL SHPDPR(X,Q2,XPPR)
C*end*
          DO 230 KFL=-6,6
            XPQ(KFL)=XPPR(KFL)
  230     CONTINUE
        ELSEIF(MSTP(52).EQ.2) THEN
C...Call PDFLIB parton distributions.
          PARM(1)='NPTYPE'
          VALUE(1)=1
          PARM(2)='NGROUP'
          VALUE(2)=MSTP(51)/1000
          PARM(3)='NSET'
          VALUE(3)=MOD(MSTP(51),1000)
          IF(MINT(93).NE.1000000+MSTP(51)) THEN
            CALL SDFSET(PARM,VALUE)
            MINT(93)=1000000+MSTP(51)
          ENDIF
          XX=X
          QQ=SQRT(MAX(0D0,Q2MIN,Q2))
          IF(MSTP(57).EQ.0) QQ=SQRT(Q2MIN)
          CALL PTRUCTM(XX,QQ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU)
          VINT(231)=Q2MIN
          XPQ(0)=GLU
          XPQ(1)=DNV+DSEA
          XPQ(-1)=DSEA
          XPQ(2)=UPV+USEA
          XPQ(-2)=USEA
          XPQ(3)=STR
          XPQ(-3)=STR
          XPQ(4)=CHM
          XPQ(-4)=CHM
          XPQ(5)=BOT
          XPQ(-5)=BOT
          XPQ(6)=TOP
          XPQ(-6)=TOP
        ELSE
          WRITE(MSTU(11),5200) KF,MSTP(52),MSTP(51)
        ENDIF
      ENDIF
 
C...Isospin average for pi0/gammaVDM.
      IF(KFA.EQ.111.OR.(KFA.EQ.22.AND.MINT(109).EQ.2)) THEN
        IF(KFA.EQ.22.AND.MSTP(55).GE.5.AND.MSTP(55).LE.12) THEN
          XPV=XPQ(2)-XPQ(1)
          XPQ(2)=XPQ(1)
          XPQ(-2)=XPQ(-1)
        ELSE
          XPS=0.5D0*(XPQ(1)+XPQ(-2))
          XPV=0.5D0*(XPQ(2)+XPQ(-1))-XPS
          XPQ(2)=XPS
          XPQ(-1)=XPS
        ENDIF
        IF(KFA.EQ.22.AND.MINT(105).LE.223) THEN
          XPQ(1)=XPQ(1)+0.2D0*XPV
          XPQ(-1)=XPQ(-1)+0.2D0*XPV
          XPQ(2)=XPQ(2)+0.8D0*XPV
          XPQ(-2)=XPQ(-2)+0.8D0*XPV
        ELSEIF(KFA.EQ.22.AND.MINT(105).EQ.333) THEN
          XPQ(3)=XPQ(3)+XPV
          XPQ(-3)=XPQ(-3)+XPV
        ELSEIF(KFA.EQ.22.AND.MINT(105).EQ.443) THEN
          XPQ(4)=XPQ(4)+XPV
          XPQ(-4)=XPQ(-4)+XPV
          IF(MSTP(55).GE.9) THEN
            DO 240 KFL=-6,6
              XPQ(KFL)=0D0
  240       CONTINUE
          ENDIF
        ELSE
          XPQ(1)=XPQ(1)+0.5D0*XPV
          XPQ(-1)=XPQ(-1)+0.5D0*XPV
          XPQ(2)=XPQ(2)+0.5D0*XPV
          XPQ(-2)=XPQ(-2)+0.5D0*XPV
        ENDIF
 
C...Rescale for gammaVDM by effective gamma -> rho coupling.
C+++Do not rescale?
        IF(KFA.EQ.22.AND.MINT(109).EQ.2.AND..NOT.(MSTP(56).EQ.1
     &  .AND.MSTP(55).GE.5.AND.MSTP(55).LE.12)) THEN
          DO 250 KFL=-6,6
            XPQ(KFL)=VINT(281)*XPQ(KFL)
  250     CONTINUE
          VINT(232)=VINT(281)*XPV
        ENDIF
 
C...Simple recipes for kaons.
      ELSEIF(KFA.EQ.321) THEN
        XPQ(-3)=XPQ(-3)+XPQ(-1)-XPQ(1)
        XPQ(-1)=XPQ(1)
      ELSEIF(KFA.EQ.130.OR.KFA.EQ.310) THEN
        XPS=0.5D0*(XPQ(1)+XPQ(-2))
        XPV=0.5D0*(XPQ(2)+XPQ(-1))-XPS
        XPQ(2)=XPS
        XPQ(-1)=XPS
        XPQ(1)=XPQ(1)+0.5D0*XPV
        XPQ(-1)=XPQ(-1)+0.5D0*XPV
        XPQ(3)=XPQ(3)+0.5D0*XPV
        XPQ(-3)=XPQ(-3)+0.5D0*XPV
 
C...Isospin conjugation for neutron.
      ELSEIF(KFA.EQ.2112) THEN
        XPS=XPQ(1)
        XPQ(1)=XPQ(2)
        XPQ(2)=XPS
        XPS=XPQ(-1)
        XPQ(-1)=XPQ(-2)
        XPQ(-2)=XPS
 
C...Simple recipes for hyperon (average valence parton distribution).
      ELSEIF(KFA.EQ.3122.OR.KFA.EQ.3112.OR.KFA.EQ.3212.OR.KFA.EQ.3222
     &  .OR.KFA.EQ.3312.OR.KFA.EQ.3322.OR.KFA.EQ.3334) THEN
        XPVAL=(XPQ(1)+XPQ(2)-XPQ(-1)-XPQ(-2))/3D0
        XPSEA=0.5D0*(XPQ(-1)+XPQ(-2))
        XPQ(1)=XPSEA
        XPQ(2)=XPSEA
        XPQ(-1)=XPSEA
        XPQ(-2)=XPSEA
        XPQ(KFA/1000)=XPQ(KFA/1000)+XPVAL
        XPQ(MOD(KFA/100,10))=XPQ(MOD(KFA/100,10))+XPVAL
        XPQ(MOD(KFA/10,10))=XPQ(MOD(KFA/10,10))+XPVAL
      ENDIF
 
C...Charge conjugation for antiparticle.
      IF(KF.LT.0) THEN
        DO 260 KFL=1,25
          IF(KFL.EQ.21.OR.KFL.EQ.22.OR.KFL.EQ.23.OR.KFL.EQ.25) GOTO 260
          XPS=XPQ(KFL)
          XPQ(KFL)=XPQ(-KFL)
          XPQ(-KFL)=XPS
  260   CONTINUE
      ENDIF
 
C...Allow gluon also in position 21.
      XPQ(21)=XPQ(0)
 
C...Check positivity and reset above maximum allowed flavour.
      DO 270 KFL=-25,25
        XPQ(KFL)=MAX(0D0,XPQ(KFL))
        IF(IABS(KFL).GT.MSTP(58).AND.IABS(KFL).LE.8) XPQ(KFL)=0D0
  270 CONTINUE
 
C...Formats for error printouts.
 5000 FORMAT(' Error: x value outside physical range; x =',1P,D12.3)
 5100 FORMAT(' Error: illegal particle code for parton distribution;',
     &' KF =',I5)
 5200 FORMAT(' Error: unknown parton distribution; KF, library, set =',
     &3I5)
 
      RETURN
      END
 
