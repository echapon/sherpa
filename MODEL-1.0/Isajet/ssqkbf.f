CDECK  ID>, SSQKBF. 
        SUBROUTINE SSQKBF
C-----------------------------------------------------------------------
C
C        This program gives squark branching fractions to gauginos
C        according to Baer,Barger,Karatas,Tata (Phys.Rev.D36,96(1987)
C        Updated for b_1,b_2 and non-degenerate sq masses 8/13/96
C        Baer's SQUBF
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
      COMMON/SSLUN/LOUT
      INTEGER LOUT
      SAVE /SSLUN/
C          MXSS         =  maximum number of modes
C          NSSMOD       = number of modes
C          ISSMOD       = initial particle
C          JSSMOD       = final particles
C          GSSMOD       = width
C          BSSMOD       = branching ratio
C          MSSMOD       = decay matrix element pointer
C          LSSMOD       = logical flag used internally by SSME3
      INTEGER MXSS
      PARAMETER (MXSS=1000)
      COMMON/SSMODE/NSSMOD,ISSMOD(MXSS),JSSMOD(5,MXSS),GSSMOD(MXSS)
     $,BSSMOD(MXSS),MSSMOD(MXSS),LSSMOD
      INTEGER NSSMOD,ISSMOD,JSSMOD,MSSMOD
      REAL GSSMOD,BSSMOD
      LOGICAL LSSMOD
      SAVE /SSMODE/
C          Standard model parameters
C          AMUP,...,AMTP        = quark masses
C          AME,AMMU,AMTAU       = lepton masses
C          AMW,AMZ              = W,Z masses
C          GAMW,GAMZ            = W,Z widths
C          ALFAEM,SN2THW,ALFA3  = SM couplings
C          ALQCD4               = 4 flavor lambda
      COMMON/SSSM/AMUP,AMDN,AMST,AMCH,AMBT,AMTP,AME,AMMU,AMTAU
     $,AMW,AMZ,GAMW,GAMZ,ALFAEM,SN2THW,ALFA2,ALFA3,ALQCD4
      REAL AMUP,AMDN,AMST,AMCH,AMBT,AMTP,AME,AMMU,AMTAU
     $,AMW,AMZ,GAMW,GAMZ,ALFAEM,SN2THW,ALFA2,ALFA3,ALQCD4
      SAVE /SSSM/
C          SUSY parameters
C          AMGLSS               = gluino mass
C          AMULSS               = up-left squark mass
C          AMELSS               = left-selectron mass
C          AMERSS               = right-slepton mass
C          AMNiSS               = sneutrino mass for generation i
C          TWOM1                = Higgsino mass = - mu
C          RV2V1                = ratio v2/v1 of vev's
C          AMTLSS,AMTRSS        = left,right stop masses
C          AMT1SS,AMT2SS        = light,heavy stop masses
C          AMBLSS,AMBRSS        = left,right sbottom masses
C          AMB1SS,AMB2SS        = light,heavy sbottom masses
C          AMLLSS,AMLRSS        = left,right stau masses
C          AML1SS,AML2SS        = light,heavy stau masses
C          AMZiSS               = signed mass of Zi
C          ZMIXSS               = Zi mixing matrix
C          AMWiSS               = signed Wi mass
C          GAMMAL,GAMMAR        = Wi left, right mixing angles
C          AMHL,AMHH,AMHA       = neutral Higgs h0, H0, A0 masses
C          AMHC                 = charged Higgs H+ mass
C          ALFAH                = Higgs mixing angle
C          AAT                  = stop trilinear term
C          THETAT               = stop mixing angle
C          AAB                  = sbottom trilinear term
C          THETAB               = sbottom mixing angle
C          AAL                  = stau trilinear term
C          THETAL               = stau mixing angle
C          AMGVSS               = gravitino mass
C          MTQ                  = top mass at MSUSY
C          MBQ                  = bottom mass at MSUSY
C          MLQ                  = tau mass at MSUSY
C          FBMA                 = b-Yukawa at mA scale
C          VUQ                  = Hu vev at MSUSY
C          VDQ                  = Hd vev at MSUSY
      COMMON/SSPAR/AMGLSS,AMULSS,AMURSS,AMDLSS,AMDRSS,AMSLSS
     $,AMSRSS,AMCLSS,AMCRSS,AMBLSS,AMBRSS,AMB1SS,AMB2SS
     $,AMTLSS,AMTRSS,AMT1SS,AMT2SS,AMELSS,AMERSS,AMMLSS,AMMRSS
     $,AMLLSS,AMLRSS,AML1SS,AML2SS,AMN1SS,AMN2SS,AMN3SS
     $,TWOM1,RV2V1,AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,ZMIXSS(4,4)
     $,AMW1SS,AMW2SS
     $,GAMMAL,GAMMAR,AMHL,AMHH,AMHA,AMHC,ALFAH,AAT,THETAT
     $,AAB,THETAB,AAL,THETAL,AMGVSS,MTQ,MBQ,MLQ,FBMA,
     $VUQ,VDQ
      REAL AMGLSS,AMULSS,AMURSS,AMDLSS,AMDRSS,AMSLSS
     $,AMSRSS,AMCLSS,AMCRSS,AMBLSS,AMBRSS,AMB1SS,AMB2SS
     $,AMTLSS,AMTRSS,AMT1SS,AMT2SS,AMELSS,AMERSS,AMMLSS,AMMRSS
     $,AMLLSS,AMLRSS,AML1SS,AML2SS,AMN1SS,AMN2SS,AMN3SS
     $,TWOM1,RV2V1,AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,ZMIXSS
     $,AMW1SS,AMW2SS
     $,GAMMAL,GAMMAR,AMHL,AMHH,AMHA,AMHC,ALFAH,AAT,THETAT
     $,AAB,THETAB,AAL,THETAL,AMGVSS,MTQ,MBQ,MLQ,FBMA,VUQ,VDQ
      REAL AMZISS(4)
      EQUIVALENCE (AMZISS(1),AMZ1SS)
      SAVE /SSPAR/
C          SM ident code definitions. These are standard ISAJET but
C          can be changed.
      INTEGER IDUP,IDDN,IDST,IDCH,IDBT,IDTP
      INTEGER IDNE,IDE,IDNM,IDMU,IDNT,IDTAU
      INTEGER IDGL,IDGM,IDW,IDZ,IDH
      PARAMETER (IDUP=1,IDDN=2,IDST=3,IDCH=4,IDBT=5,IDTP=6)
      PARAMETER (IDNE=11,IDE=12,IDNM=13,IDMU=14,IDNT=15,IDTAU=16)
      PARAMETER (IDGL=9,IDGM=10,IDW=80,IDZ=90,IDH=81)
C          SUSY ident code definitions. They are chosen to be similar
C          to those in versions < 6.50 but may be changed.
      INTEGER ISUPL,ISDNL,ISSTL,ISCHL,ISBT1,ISTP1
      INTEGER ISNEL,ISEL,ISNML,ISMUL,ISNTL,ISTAU1
      INTEGER ISUPR,ISDNR,ISSTR,ISCHR,ISBT2,ISTP2
      INTEGER ISNER,ISER,ISNMR,ISMUR,ISNTR,ISTAU2
      INTEGER ISZ1,ISZ2,ISZ3,ISZ4,ISW1,ISW2,ISGL
      INTEGER ISHL,ISHH,ISHA,ISHC
      INTEGER ISGRAV
      INTEGER IDTAUL,IDTAUR
      PARAMETER (ISUPL=21,ISDNL=22,ISSTL=23,ISCHL=24,ISBT1=25,ISTP1=26)
      PARAMETER (ISNEL=31,ISEL=32,ISNML=33,ISMUL=34,ISNTL=35,ISTAU1=36)
      PARAMETER (ISUPR=41,ISDNR=42,ISSTR=43,ISCHR=44,ISBT2=45,ISTP2=46)
      PARAMETER (ISNER=51,ISER=52,ISNMR=53,ISMUR=54,ISNTR=55,ISTAU2=56)
      PARAMETER (ISGL=29)
      PARAMETER (ISZ1=30,ISZ2=40,ISZ3=50,ISZ4=60,ISW1=39,ISW2=49)
      PARAMETER (ISHL=82,ISHH=83,ISHA=84,ISHC=86)
      PARAMETER (ISGRAV=91)
      PARAMETER (IDTAUL=10016,IDTAUR=20016)
C
      COMPLEX ZI,ZONE,ZA,ZB,ZAUIZ,ZADIZ,ZBUIZ,ZBDIZ
      DOUBLE PRECISION SSALFS
      REAL SSXLAM,WID,AUIZS,ADIZS,BUIZS,BDIZS
      REAL PI,SR2,G,GP,COSA,SINA,SNZI,THIZ
     $,TANB,COTB,XM,YM,THX,THY,FT,FB
      REAL MZIZ,CS2THW,TN2THW,BETA,BH,A,AS
      INTEGER IZ
      REAL MW1,MW2,SNW1,SNW2,COST,SINT,COSB,SINB
      REAL AWD(2),AWU(2),BW(2),BWP(2)
      INTEGER ISZIZ(4)
      DATA ZI/(0.,1.)/,ZONE/(1.,0.)/
C
C          Partly duplicated from SSMASS
C
      PI=4.*ATAN(1.)
      SR2=SQRT(2.)
      G=SQRT(4*PI*ALFAEM/SN2THW)
      GP=G*SQRT(SN2THW/(1.-SN2THW))
      CS2THW=1.-SN2THW
      TN2THW=SN2THW/CS2THW
      TANB=1./RV2V1
      COTB=RV2V1
      BETA=ATAN(TANB)
C          Reconstruct masses from SSMASS
      MW1=ABS(AMW1SS)
      MW2=ABS(AMW2SS)
      COST=COS(THETAT)
      SINT=SIN(THETAT)
      COSB=COS(THETAB)
      SINB=SIN(THETAB)
      COSA=COS(ALFAH)
      SINA=SIN(ALFAH)
      SNW1=SIGN(1.,AMW1SS)
      SNW2=SIGN(1.,AMW2SS)
      XM=1./TAN(GAMMAL)
      YM=1./TAN(GAMMAR)
      THX=SIGN(1.,XM)
      THY=SIGN(1.,YM)
      FB=G*MBQ/SR2/AMW/COS(BETA)
      FT=G*MTQ/SR2/AMW/SIN(BETA)
      AWD(1)=-G*SNW1*SIN(GAMMAR)
      AWD(2)=-G*SNW2*THY*COS(GAMMAR)
      AWU(1)=-G*SIN(GAMMAL)
      AWU(2)=-G*THX*COS(GAMMAL)
      BW(1)=-FT*SNW1*COS(GAMMAR)
      BW(2)=FT*SNW2*THY*SIN(GAMMAR)
      BWP(1)=-FB*COS(GAMMAL)
      BWP(2)=FB*THX*SIN(GAMMAL)
C
C          Compute squark branching fractions to zi
C
      ISZIZ(1)=ISZ1
      ISZIZ(2)=ISZ2
      ISZIZ(3)=ISZ3
      ISZIZ(4)=ISZ4
      DO 100 IZ=1,4
        MZIZ=ABS(AMZISS(IZ))
        SNZI=SIGN(1.,AMZISS(IZ))
        IF (SNZI.EQ.1.) THEN
           THIZ=0.
        ELSE
           THIZ=1.
        END IF
        ZAUIZ=ZI**(THIZ-1.)*(-1)*SNZI
     $  *(G/SR2*ZMIXSS(3,IZ)+GP/3./SR2*ZMIXSS(4,IZ))
        ZBUIZ=ZI**(THIZ-1.)*4*GP*ZMIXSS(4,IZ)/3./SR2
        ZADIZ=ZI**(THIZ-1.)*(-1)*SNZI
     $  *(-G/SR2*ZMIXSS(3,IZ)+GP/3./SR2*ZMIXSS(4,IZ))
        ZBDIZ=ZI**(THIZ-1.)*(-2)*GP*ZMIXSS(4,IZ)/3./SR2
        AUIZS=ZAUIZ*CONJG(ZAUIZ)
        ADIZS=ZADIZ*CONJG(ZADIZ)
        BUIZS=ZBUIZ*CONJG(ZBUIZ)
        BDIZS=ZBDIZ*CONJG(ZBDIZ)
C          squark --> q + qb + zi, q = u, d, s
        IF (AMULSS.GT.MZIZ) THEN
          WID=AUIZS*AMULSS*(1.-MZIZ**2/AMULSS**2)**2/16./PI
          CALL SSSAVE(ISUPL,WID,ISZIZ(IZ),IDUP,0,0,0)
        END IF
        IF (AMDLSS.GT.MZIZ) THEN
          WID=ADIZS*AMDLSS*(1.-MZIZ**2/AMDLSS**2)**2/16./PI
          CALL SSSAVE(ISDNL,WID,ISZIZ(IZ),IDDN,0,0,0)
        END IF
        IF (AMSLSS.GT.MZIZ) THEN
          WID=ADIZS*AMSLSS*(1.-MZIZ**2/AMSLSS**2)**2/16./PI
          CALL SSSAVE(ISSTL,WID,ISZIZ(IZ),IDST,0,0,0)
        END IF
        IF (AMURSS.GT.MZIZ) THEN
          WID=BUIZS*AMURSS*(1.-MZIZ**2/AMURSS**2)**2/16./PI
          CALL SSSAVE(ISUPR,WID,ISZIZ(IZ),IDUP,0,0,0)
        END IF
        IF (AMDRSS.GT.MZIZ) THEN
          WID=BDIZS*AMDRSS*(1.-MZIZ**2/AMDRSS**2)**2/16./PI
          CALL SSSAVE(ISDNR,WID,ISZIZ(IZ),IDDN,0,0,0)
        END IF
        IF (AMSRSS.GT.MZIZ) THEN
          WID=BDIZS*AMSRSS*(1.-MZIZ**2/AMSRSS**2)**2/16./PI
          CALL SSSAVE(ISSTR,WID,ISZIZ(IZ),IDST,0,0,0)
        END IF
C          squark --> q + zi, q = c
        IF (AMCLSS.GT.(MZIZ+AMCH)) THEN
          WID=AUIZS*AMCLSS*(1.-MZIZ**2/AMCLSS**2-AMCH**2/AMCLSS**2)
     $    *SQRT(SSXLAM(1.,MZIZ**2/AMCLSS**2,AMCH**2/AMCLSS**2))/16./PI
          CALL SSSAVE(ISCHL,WID,ISZIZ(IZ),IDCH,0,0,0)
        END IF
        IF (AMCRSS.GT.(MZIZ+AMCH)) THEN
          WID=BUIZS*AMCRSS*(1.-MZIZ**2/AMCRSS**2-AMCH**2/AMCRSS**2)
     $    *SQRT(SSXLAM(1.,MZIZ**2/AMCRSS**2,AMCH**2/AMCRSS**2))/16./PI
          CALL SSSAVE(ISCHR,WID,ISZIZ(IZ),IDCH,0,0,0)
        END IF
C          sbottom_1 --> b + zi
        IF (AMB1SS.GT.(MZIZ+AMBT)) THEN
          ZA=(ZI*ZADIZ-FB*ZMIXSS(2,IZ)*ZI**THIZ)*COSB/2.-
     $       (ZI*ZBDIZ-FB*ZMIXSS(2,IZ)*(-ZI)**THIZ)*SINB/2.
          ZB=(-ZI*ZADIZ-FB*ZMIXSS(2,IZ)*ZI**THIZ)*COSB/2.-
     $       (ZI*ZBDIZ+FB*ZMIXSS(2,IZ)*(-ZI)**THIZ)*SINB/2.
          WID=(ZA*CONJG(ZA)*(AMB1SS**2-(AMBT+MZIZ)**2)+
     $     ZB*CONJG(ZB)*(AMB1SS**2-(MZIZ-AMBT)**2))/8./PI/AMB1SS
     $    *SQRT(SSXLAM(1.,MZIZ**2/AMB1SS**2,AMBT**2/AMB1SS**2))
          CALL SSSAVE(ISBT1,WID,ISZIZ(IZ),IDBT,0,0,0)
        END IF
C          sbottom_2 --> b + zi
        IF (AMB2SS.GT.(MZIZ+AMBT)) THEN
          ZA=(ZI*ZADIZ-FB*ZMIXSS(2,IZ)*ZI**THIZ)*SINB/2.+
     $       (ZI*ZBDIZ-FB*ZMIXSS(2,IZ)*(-ZI)**THIZ)*COSB/2.
          ZB=(-ZI*ZADIZ-FB*ZMIXSS(2,IZ)*ZI**THIZ)*SINB/2.+
     $       (ZI*ZBDIZ+FB*ZMIXSS(2,IZ)*(-ZI)**THIZ)*COSB/2.
          WID=(ZA*CONJG(ZA)*(AMB2SS**2-(AMBT+MZIZ)**2)+
     $     ZB*CONJG(ZB)*(AMB2SS**2-(MZIZ-AMBT)**2))/8./PI/AMB2SS
     $    *SQRT(SSXLAM(1.,MZIZ**2/AMB2SS**2,AMBT**2/AMB2SS**2))
          CALL SSSAVE(ISBT2,WID,ISZIZ(IZ),IDBT,0,0,0)
        END IF
100   CONTINUE
C
C          Compute squark branching fractions to gluinos
C
      IF (AMULSS.GT.AMGLSS) THEN
        WID=2*SSALFS(DBLE(AMULSS**2))*AMULSS*
     $  (1.-AMGLSS**2/AMULSS**2)**2/3.
        CALL SSSAVE(ISUPL,WID,ISGL,IDUP,0,0,0)
      END IF
      IF (AMDLSS.GT.AMGLSS) THEN
        WID=2*SSALFS(DBLE(AMDLSS**2))*AMDLSS*
     $  (1.-AMGLSS**2/AMDLSS**2)**2/3.
        CALL SSSAVE(ISDNL,WID,ISGL,IDDN,0,0,0)
      END IF
      IF (AMSLSS.GT.AMGLSS) THEN
        WID=2*SSALFS(DBLE(AMSLSS**2))*AMSLSS*
     $  (1.-AMGLSS**2/AMSLSS**2)**2/3.
        CALL SSSAVE(ISSTL,WID,ISGL,IDST,0,0,0)
      END IF
      IF (AMURSS.GT.AMGLSS) THEN
        WID=2*SSALFS(DBLE(AMURSS**2))*AMURSS*
     $  (1.-AMGLSS**2/AMURSS**2)**2/3.
        CALL SSSAVE(ISUPR,WID,ISGL,IDUP,0,0,0)
      END IF
      IF (AMDRSS.GT.AMGLSS) THEN
        WID=2*SSALFS(DBLE(AMDRSS**2))*AMDRSS*
     $  (1.-AMGLSS**2/AMDRSS**2)**2/3.
        CALL SSSAVE(ISDNR,WID,ISGL,IDDN,0,0,0)
      END IF
      IF (AMSRSS.GT.AMGLSS) THEN
        WID=2*SSALFS(DBLE(AMSRSS**2))*AMSRSS*
     $  (1.-AMGLSS**2/AMSRSS**2)**2/3.
        CALL SSSAVE(ISSTR,WID,ISGL,IDST,0,0,0)
      END IF
C
      IF (AMCLSS.GT.(AMGLSS+AMCH)) THEN
        WID=2*SSALFS(DBLE(AMCLSS**2))*AMCLSS*(1.-AMGLSS**2/AMCLSS**2-
     $  AMCH**2/AMCLSS**2)*SQRT(SSXLAM(1.,AMGLSS**2/AMCLSS**2,
     $  AMCH**2/AMCLSS**2))/3.
        CALL SSSAVE(ISCHL,WID,ISGL,IDCH,0,0,0)
      END IF
      IF (AMCRSS.GT.(AMGLSS+AMCH)) THEN
        WID=2*SSALFS(DBLE(AMCRSS**2))*AMCRSS*(1.-AMGLSS**2/AMCRSS**2-
     $  AMCH**2/AMCRSS**2)*SQRT(SSXLAM(1.,AMGLSS**2/AMCRSS**2,
     $  AMCH**2/AMCRSS**2))/3.
        CALL SSSAVE(ISCHR,WID,ISGL,IDCH,0,0,0)
      END IF
C
      IF (AMB1SS.GT.(AMGLSS+AMBT)) THEN
        WID=2*SSALFS(DBLE(AMB1SS**2))*AMB1SS*(1.-AMGLSS**2/AMB1SS**2-
     $  AMBT**2/AMB1SS**2)*SQRT(SSXLAM(1.,AMGLSS**2/AMB1SS**2,
     $  AMBT**2/AMB1SS**2))/3.
        CALL SSSAVE(ISBT1,WID,ISGL,IDBT,0,0,0)
      END IF
C
      IF (AMB2SS.GT.(AMGLSS+AMBT)) THEN
        WID=2*SSALFS(DBLE(AMB2SS**2))*AMB2SS*(1.-AMGLSS**2/AMB2SS**2-
     $  AMBT**2/AMB2SS**2)*SQRT(SSXLAM(1.,AMGLSS**2/AMB2SS**2,
     $  AMBT**2/AMB2SS**2))/3.
        CALL SSSAVE(ISBT2,WID,ISGL,IDBT,0,0,0)
      END IF
C
C           Compute branching fractions to wi --- theta-C = 0
C
      IF (AMULSS.GT.MW1) THEN
        WID=G**2*SIN(GAMMAR)**2*AMULSS*(1.-MW1**2/AMULSS**2)**2/16./PI
        CALL SSSAVE(ISUPL,WID,ISW1,IDDN,0,0,0)
      END IF
      IF (AMCLSS.GT.MW1) THEN
        WID=G**2*SIN(GAMMAR)**2*AMCLSS*(1.-MW1**2/AMCLSS**2)**2/16./PI
        CALL SSSAVE(ISCHL,WID,ISW1,IDST,0,0,0)
      END IF
      IF (AMDLSS.GT.MW1) THEN
        WID=G**2*SIN(GAMMAL)**2*AMDLSS*(1.-MW1**2/AMDLSS**2)**2/16./PI
        CALL SSSAVE(ISDNL,WID,-ISW1,IDUP,0,0,0)
      END IF
C
      IF (AMSLSS.GT.(MW1+AMCH)) THEN
        WID=G**2*SIN(GAMMAL)**2*AMSLSS*(1.-MW1**2/AMSLSS**2
     $  -AMCH**2/AMSLSS**2)
     $  *SQRT(SSXLAM(1.,MW1**2/AMSLSS**2,AMCH**2/AMSLSS**2))/16./PI
        CALL SSSAVE(ISSTL,WID,-ISW1,IDCH,0,0,0)
      ENDIF
C
       IF (AMB1SS.GT.(MW1+AMTP)) THEN
         A=AWU(1)*COSB-BWP(1)*SINB
         AS=A*A
         WID=AMB1SS*((AS+BW(1)**2*COSB**2)*(1.-MW1**2/AMB1SS**2
     $  -AMTP**2/AMB1SS**2)-4*AMTP*MW1*BW(1)*A*COSB/AMB1SS**2)
     $   *SQRT(SSXLAM(1.,MW1**2/AMB1SS**2,AMTP**2/AMB1SS**2))/16./PI
        CALL SSSAVE(ISBT1,WID,-ISW1,IDTP,0,0,0)
      ENDIF
C
       IF (AMB2SS.GT.(MW1+AMTP)) THEN
         A=AWU(1)*SINB+BWP(1)*COSB
         AS=A*A
         WID=AMB2SS*((AS+BW(1)**2*SINB**2)*(1.-MW1**2/AMB2SS**2
     $  -AMTP**2/AMB2SS**2)-4*AMTP*MW1*BW(1)*A*SINB/AMB2SS**2)
     $   *SQRT(SSXLAM(1.,MW1**2/AMB2SS**2,AMTP**2/AMB2SS**2))/16./PI
        CALL SSSAVE(ISBT2,WID,-ISW1,IDTP,0,0,0)
      ENDIF
C
      IF (AMULSS.GT.MW2) THEN
        WID=G**2*COS(GAMMAR)**2*AMULSS*(1.-MW2**2/AMULSS**2)**2/16./PI
        CALL SSSAVE(ISUPL,WID,ISW2,IDDN,0,0,0)
      END IF
      IF (AMCLSS.GT.MW2) THEN
        WID=G**2*COS(GAMMAR)**2*AMCLSS*(1.-MW2**2/AMCLSS**2)**2/16./PI
        CALL SSSAVE(ISCHL,WID,ISW2,IDST,0,0,0)
      END IF
      IF (AMDLSS.GT.MW2) THEN
        WID=G**2*COS(GAMMAL)**2*AMDLSS*(1.-MW2**2/AMDLSS**2)**2/16./PI
        CALL SSSAVE(ISDNL,WID,-ISW2,IDUP,0,0,0)
      END IF
C
      IF (AMSLSS.GT.(MW2+AMCH)) THEN
        WID=G**2*COS(GAMMAL)**2*AMSLSS*(1.-MW2**2/AMSLSS**2
     $  -AMCH**2/AMSLSS**2)
     $  *SQRT(SSXLAM(1.,MW2**2/AMSLSS**2,AMCH**2/AMSLSS**2))/16./PI
        CALL SSSAVE(ISSTL,WID,-ISW2,IDCH,0,0,0)
      ENDIF
C
      IF (AMB1SS.GT.(MW2+AMTP)) THEN
         A=AWU(2)*COSB-BWP(2)*SINB
         AS=A*A
         WID=AMB1SS*((AS+BW(2)**2*COSB**2)*(1.-MW2**2/AMB1SS**2
     $  -AMTP**2/AMB1SS**2)-4*AMTP*MW2*BW(2)*A*COSB/AMB1SS**2)
     $   *SQRT(SSXLAM(1.,MW2**2/AMB1SS**2,AMTP**2/AMB1SS**2))/16./PI
        CALL SSSAVE(ISBT1,WID,-ISW2,IDTP,0,0,0)
      ENDIF
C
      IF (AMB2SS.GT.(MW2+AMTP)) THEN
         A=AWU(2)*SINB+BWP(2)*COSB
         AS=A*A
         WID=AMB2SS*((AS+BW(2)**2*SINB**2)*(1.-MW2**2/AMB2SS**2
     $  -AMTP**2/AMB2SS**2)-4*AMTP*MW2*BW(2)*A*SINB/AMB2SS**2)
     $   *SQRT(SSXLAM(1.,MW2**2/AMB2SS**2,AMTP**2/AMB2SS**2))/16./PI
        CALL SSSAVE(ISBT2,WID,-ISW2,IDTP,0,0,0)
      ENDIF
C
      IF (AMB1SS.GT.(AMW+AMT1SS)) THEN
        WID=G**2*COST**2*COSB**2*(SSXLAM(AMB1SS**2,AMW**2,
     $   AMT1SS**2))**1.5/32./PI/AMB1SS**3/AMW**2
        CALL SSSAVE(ISBT1,WID,-IDW,ISTP1,0,0,0)
      ENDIF
C
      IF (AMB1SS.GT.(AMW+AMT2SS)) THEN
        WID=G**2*SINT**2*COSB**2*(SSXLAM(AMB1SS**2,AMW**2,
     $   AMT2SS**2))**1.5/32./PI/AMB1SS**3/AMW**2
        CALL SSSAVE(ISBT1,WID,-IDW,ISTP2,0,0,0)
      ENDIF
C
      IF (AMB2SS.GT.(AMW+AMT1SS)) THEN
        WID=G**2*COST**2*SINB**2*(SSXLAM(AMB2SS**2,AMW**2,
     $   AMT1SS**2))**1.5/32./PI/AMB2SS**3/AMW**2
        CALL SSSAVE(ISBT2,WID,-IDW,ISTP1,0,0,0)
      ENDIF
C
      IF (AMB2SS.GT.(AMW+AMT2SS)) THEN
        WID=G**2*SINT**2*SINB**2*(SSXLAM(AMB2SS**2,AMW**2,
     $   AMT2SS**2))**1.5/32./PI/AMB2SS**3/AMW**2
        CALL SSSAVE(ISBT2,WID,-IDW,ISTP2,0,0,0)
      ENDIF
C
      IF (AMB2SS.GT.(AMZ+AMB1SS)) THEN
        WID=G**2*COSB**2*SINB**2*(SSXLAM(AMB2SS**2,AMZ**2,
     $   AMB1SS**2))**1.5/64./PI/AMB2SS**3/AMZ**2/CS2THW
        CALL SSSAVE(ISBT2,WID,IDZ,ISBT1,0,0,0)
      ENDIF
C
      IF (AMB2SS.GT.(AMHL+AMB1SS)) THEN
        BH=G*AMW*SIN(BETA-ALFAH)*(-1.+TN2THW/3.)*SINB*COSB/2.+G*
     $  AMBT*(TWOM1*COSA+AAB*SINA)*COS(2*THETAB)/2./AMW/COS(BETA)
        WID=BH**2*SQRT(SSXLAM(AMB2SS**2,AMHL**2,AMB1SS**2))/
     $      16./PI/AMB2SS**3
        CALL SSSAVE(ISBT2,WID,ISHL,ISBT1,0,0,0)
      ENDIF
C
      IF (AMB2SS.GT.(AMHA+AMB1SS)) THEN
        BH=G*AMBT*(TWOM1-AAB*TANB)/2./AMW
        WID=BH**2*SQRT(SSXLAM(AMB2SS**2,AMHA**2,AMB1SS**2))/
     $      16./PI/AMB2SS**3
        CALL SSSAVE(ISBT2,WID,ISHA,ISBT1,0,0,0)
      ENDIF
C
      IF (AMB2SS.GT.(AMHH+AMB1SS)) THEN
        BH=-G*AMW*COS(BETA-ALFAH)*(-1.+TN2THW/3.)*SINB*COSB/2.+G*
     $  AMBT*(-TWOM1*SINA+AAB*COSA)*COS(2*THETAB)/2./AMW/COS(BETA)
        WID=BH**2*SQRT(SSXLAM(AMB2SS**2,AMHH**2,AMB1SS**2))/
     $      16./PI/AMB2SS**3
        CALL SSSAVE(ISBT2,WID,ISHH,ISBT1,0,0,0)
      ENDIF
C
C     b_i -> H^- t_i
C
      IF (AMB1SS.GT.(AMT1SS+AMHC)) THEN
        A=G/SR2/AMW*(AMTP*AMBT*(COTB+TANB)*SINT*SINB+
     $(AMBT**2*TANB+AMTP**2*COTB-AMW**2*SIN(2*BETA))*
     $COST*COSB-AMTP*(TWOM1-AAT*COTB)*SINT*COSB-AMBT*
     $(TWOM1-AAB*TANB)*SINB*COST)
        AS=A*A
        WID=AS*SQRT(SSXLAM(AMB1SS**2,AMT1SS**2,AMHC**2))/
     $      16./PI/AMB1SS**3
        CALL SSSAVE(ISBT1,WID,-ISHC,ISTP1,0,0,0)
      END IF
C
      IF (AMB1SS.GT.(AMT2SS+AMHC)) THEN
        A=G/SR2/AMW*(-AMTP*AMBT*(COTB+TANB)*COST*SINT+
     $(AMBT**2*TANB+AMTP**2*COTB-AMW**2*SIN(2*BETA))*
     $SINT*COSB+AMTP*(TWOM1-AAT*COTB)*COST*COSB-AMBT*
     $(TWOM1-AAB*TANB)*SINT*SINB)
        AS=A*A
        WID=AS*SQRT(SSXLAM(AMB1SS**2,AMT2SS**2,AMHC**2))/
     $      16./PI/AMB1SS**3
        CALL SSSAVE(ISBT1,WID,-ISHC,ISTP2,0,0,0)
      END IF
C
      IF (AMB2SS.GT.(AMT1SS+AMHC)) THEN
        A=G/SR2/AMW*(-AMTP*AMBT*(COTB+TANB)*SINT*COSB+
     $(AMBT**2*TANB+AMTP**2*COTB-AMW**2*SIN(2*BETA))*
     $COST*SINB-AMTP*(TWOM1-AAT*COTB)*SINT*SINB+AMBT*
     $(TWOM1-AAB*TANB)*COST*COSB)
        AS=A*A
        WID=AS*SQRT(SSXLAM(AMB2SS**2,AMT1SS**2,AMHC**2))/
     $      16./PI/AMB2SS**3
        CALL SSSAVE(ISBT2,WID,-ISHC,ISTP1,0,0,0)
      END IF
C
      IF (AMB2SS.GT.(AMT2SS+AMHC)) THEN
        A=G/SR2/AMW*(AMTP*AMBT*(COTB+TANB)*COST*COSB+
     $(AMBT**2*TANB+AMTP**2*COTB-AMW**2*SIN(2*BETA))*
     $SINT*SINB+AMTP*(TWOM1-AAT*COTB)*SINB*COST+AMBT*
     $(TWOM1-AAB*TANB)*COSB*SINT)
        AS=A*A
        WID=AS*SQRT(SSXLAM(AMB2SS**2,AMT2SS**2,AMHC**2))/
     $      16./PI/AMB2SS**3
        CALL SSSAVE(ISBT2,WID,-ISHC,ISTP2,0,0,0)
      END IF
C
C          Normalize branching ratios
C
      CALL SSNORM(ISUPL)
      CALL SSNORM(ISDNL)
      CALL SSNORM(ISSTL)
      CALL SSNORM(ISCHL)
      CALL SSNORM(ISBT1)
      CALL SSNORM(ISUPR)
      CALL SSNORM(ISDNR)
      CALL SSNORM(ISSTR)
      CALL SSNORM(ISCHR)
      CALL SSNORM(ISBT2)
C
       RETURN
       END
