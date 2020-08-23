      SUBROUTINE ZPZE(NBH,nc,NHE,NKHE,NPF,NPNE,nr,NVHE,NW,nx,
     '  CURVCORRECT,SE,ZA,ZE,ZP,ERROR,*)

C#### Subroutine: ZPZE
C###  Description:
C###    ZPZE transfers global node parameters ZP(nk,nv,nh,np,nc) and
C###    auxiliary element parameters ZA(na,nh,nc,ne) to element node
C###    parameters ZE(ns,nhx) for nc.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBH(NHM,NCM),nc,NHE,NKHE(NKM,NNM,NHM),NPF(9),
     '  NPNE(NNM,NBFM),nr,NVHE(NNM,NBFM,NHM),NW,nx
      REAL*8 CURVCORRECT(2,2,NNM),SE(NSM,NBFM),ZA(NAM,NHM,NCM),
     '  ZE(NSM,NHM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER MK,na,nb,nh,nhx,NITB,nk,nkk,NK1,NKTB,nn,NNK,NNS,NNTB,np,
     '  ns,ns1,ns2,ns3,ns4,ns5,ns6,ns7,ns8,ns9,ns10,ns11,ns12,nv
      REAL*8 SEN,SUM,THETA1,THETA2,TOLERANCE,Y,YC,YS
      LOGICAL INCR1,INCR12,INCR3,INCR5,INCR56,INCR7

      DATA TOLERANCE/1.0D-10/

      CALL ENTERS('ZPZE',*9999)

C***  Subroutine as of 25/11/94 ajp,cpb.

      IF(ITYP1(nr,nx).EQ.4.AND.(NW.EQ.7.OR.NW.EQ.8)) THEN !Lin. Elas.
        DO nhx=1,NHE
          nh=NH_LOC(nhx,nx)
          nb=NBH(nh,nc)
          ns=0
          DO nn=1,NNT(nb)
            np=NPNE(nn,nb)
            nv=NVHE(nn,nb,nh)
            IF(nhx.LE.2) THEN
              NNK=nhx+1
C !!! ajp 25/11/94 This is not the correct way to calculate ns!
              nns=(nn-1)*NKT(0,NBH(3,nc))+NNK
              SEN=SE(nns,NBH(3,nc))
              WRITE(OP_STRING,'('' Dubious code here (ZPZE)'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ELSE
              SEN=1.0d0
            ENDIF
            DO nk=1,NKT(0,nb)
              ns=ns+1
              ZE(ns,nhx)=ZP(nk,nv,nh,np,nc)*SE(ns,nb)*SEN
            ENDDO !nk
          ENDDO !nn
        ENDDO !nhx

      ELSE !not linear elasticity (non FE40)
        DO nhx=1,NHE
          nh=NH_LOC(nhx,nx)
          nb=NBH(nh,nc)
          ns=0
          DO nn=1,NNT(nb)
            np=NPNE(nn,nb)
            nv=NVHE(nn,nb,nh)
            DO mk=1,NKT(nn,nb)
              nk=NKHE(mk,nn,nh)
              ns=ns+1

C LKC 11-MAR-98 Problem if solving with a bi-cubic field
C               (crosses set to 0) and then fitting to a bi-linear field
C               and then using EVELEC on fitted field. KTYP93 is still
C               set and the statement below tries to set crosses
C               to 0 when there really are no longer any. The problem
C               is KTYP93 does not know which basis you are now using
C
C cpb 22/3/97   Fixing below
C GMH 13/3/97   We wish to set cross implicitly to zero
C              IF((mk.EQ.NKT(0,nb)).AND.(KTYP93(nc,nr).EQ.1)) THEN
C
              IF((mk.EQ.NKT(0,nb)).AND.
     '          ((KTYP93(nc,nr).EQ.1).AND.(NKT(0,nb).GT.3))) THEN
                ZE(ns,nhx)=0.0d0
              ELSE IF(BEMCURVATURECORRECTION.AND.nc.EQ.2.AND.
     '            (nk.EQ.2.OR.nk.EQ.3)) THEN
                SUM=0.0d0
                DO nkk=1,NKT(nn,nb)-1
                  SUM=SUM+CURVCORRECT(nkk,nk-1,nn)*ZP(nkk+1,nv,nhx,np,1)
                ENDDO !nkk
                ZE(ns,nhx)=(ZP(nk,nv,nh,np,nc)+SUM)*SE(ns,nb)
              ELSE
                ZE(ns,nhx)=ZP(nk,nv,nh,np,nc)*SE(ns,nb)
              ENDIF
            ENDDO !mk
          ENDDO !nn
          DO na=1,NAT(nb)
            ZE(ns+na,nhx)=ZA(na,nh,nc)
          ENDDO !na
        ENDDO !nhx
      ENDIF !FE40
      IF(ITYP11(nr).GT.1.AND.NHE.GT.1.AND.ITYP1(nr,nx).NE.3
     '  .AND.ITYP1(nr,nx).NE.9) THEN
        IF(JTYP10.GE.2) THEN
          nb=NBH(NH_LOC(1,nx),nc)
          NK1=NKT(0,nb)
          IF(NK1.GT.1) THEN
C****       18-Feb-89: Further modification needed if nodal parameters
C****       include second derivs, in which case we need IDO(nk,nn,0,nb)
            DO nn=1,NNT(nb)
              ns1=1+(nn-1)*NK1
              DO nk=2,NK1
                ns2=nk+(nn-1)*NK1
                IF(ITYP11(nr).LT.4) THEN
                  ZE(ns2,1)=ITYP11(nr)*ZE(ns1,1)**(ITYP11(nr)-1)*
     '              ZE(ns2,1)
                ELSE IF(ITYP11(nr).EQ.4) THEN
                  YC=DCOSH(ZE(ns1,1))
                  YS=DSINH(ZE(ns1,1))
                  IF(JTYP10.EQ.2) THEN
                    ZE(ns2,1)=2.0d0*FOCUS*FOCUS*YC*YS*ZE(ns2,1)
                  ELSE IF(JTYP10.EQ.3) THEN
                    ZE(ns2,1)=FOCUS**3*YS*(3.0d0*YC*YC-1.0d0)*
     '                ZE(ns2,1)
                  ENDIF
                ENDIF
              ENDDO !nk
            ENDDO !nn
          ENDIF
          DO nn=1,NNT(nb)
            ns=1+(nn-1)*NK1
            IF(ITYP11(nr).LT.4) THEN
              ZE(ns,1)=ZE(ns,1)**ITYP10(1)
            ELSE IF(ITYP11(nr).EQ.4) THEN
              Y=DCOSH(ZE(ns,1))
              IF(JTYP10.EQ.2) THEN
                ZE(ns,1)=FOCUS*FOCUS*(Y*Y-1.0d0)
              ELSE IF(JTYP10.EQ.3) THEN
                ZE(ns,1)=FOCUS**3*Y*(Y*Y-1.0d0)
              ENDIF
            ENDIF
          ENDDO !nn
        ENDIF !jtyp10.ge.2

        IF(ITYP11(nr).LE.3) THEN      !cyl. or sph. polar
          nhx=2   !theta
        ELSE IF(ITYP11(nr).EQ.4) THEN !prolate
          nhx=3   !theta
        ENDIF
        nh=NH_LOC(nhx,nx)
        nb=NBH(nh,nc)
        NNTB=NNT(nb)
        NITB=NIT(nb)
        NKTB=NKT(0,nb)
        ns1=1
        ns2=1+  NKTB
        ns3=1+2*NKTB
        ns4=1+3*NKTB
        IF(NITB.EQ.1.OR.NITB.EQ.3.OR.(NITB.EQ.2.AND.NPF(1).EQ.1)) THEN
          IF(nntb.NE.9) THEN
C ***       Apply corrections for angular coordinates
            IF(ITYP11(nr).EQ.2.OR.ITYP11(nr).EQ.3) THEN !cyl or sph
C             Check for phi=90 (ie on axis)
              IF(ITYP11(nr).EQ.3) THEN !spherical
                IF(NPNE(1,nb).EQ.NPNE(2,nb) !1st & 2nd in common
     '            .AND.DABS(ZE(ns1,3)-PI/2.0d0).LT.TOLERANCE !phi=90
     '            .AND.NPNE(3,nb).NE.NPNE(4,nb)) THEN
                  ZE(ns1,nhx)=ZE(ns3,nhx)
                  ZE(ns2,nhx)=ZE(ns4,nhx)
                ELSE IF(NPNE(3,nb).EQ.NPNE(4,nb) !3rd & 4th in common
     '              .AND.DABS(ZE(ns3,3)-PI/2.0d0).LT.TOLERANCE !phi=90
     '              .AND.NPNE(1,nb).NE.NPNE(2,nb)) THEN
                  ZE(ns3,nhx)=ZE(ns1,nhx)
                  ZE(ns4,nhx)=ZE(ns2,nhx)
                ENDIF
              ENDIF
C             Make sure theta is increasing with Xi1
              IF(ZE(ns1,nhx).GE.ZE(ns2,nhx)-TOLERANCE) THEN
                ZE(ns2,nhx)=ZE(ns2,nhx)+2.0d0*PI
              ENDIF
              IF(ZE(ns3,nhx).GE.ZE(ns4,nhx)-TOLERANCE) THEN
                ZE(ns4,nhx)=ZE(ns4,nhx)+2.0d0*PI
              ENDIF
            ELSE IF(ITYP11(nr).EQ.4.OR.ITYP11(nr).EQ.5)THEN !prol or obl
C             Make sure theta is decreasing with Xi1
              IF(ZE(ns2,nhx).GE.ZE(ns1,nhx)-TOLERANCE) THEN
                INCR1=.TRUE. !1st vertex increased by 2*PI
                ZE(ns1,nhx)=ZE(ns1,nhx)+2.0d0*PI
              ELSE
                INCR1=.FALSE.
              ENDIF
              IF(ZE(ns4,nhx).GE.ZE(ns3,nhx)-TOLERANCE) THEN
                INCR3=.TRUE. !3rd vertex increased by 2*PI
                ZE(ns3,nhx)=ZE(ns3,nhx)+2.0d0*PI
              ELSE
                INCR3=.FALSE.
              ENDIF
              INCR12=.FALSE.
              IF(INCR1.AND..NOT.INCR3.AND.ZE(ns3,nhx).LT.PI) THEN
                ZE(ns3,nhx)=ZE(ns3,nhx)+2.0d0*PI
                ZE(ns4,nhx)=ZE(ns4,nhx)+2.0d0*PI
                INCR12=.TRUE.
              ELSE IF(INCR3.AND..NOT.INCR1.AND.ZE(ns1,nhx).LT.PI) THEN
                ZE(ns1,nhx)=ZE(ns1,nhx)+2.0d0*PI
                ZE(ns2,nhx)=ZE(ns2,nhx)+2.0d0*PI
                INCR12=.TRUE.
              ENDIF
C             Check whether top of element too different from bottom
              THETA1=0.5d0*(ZE(ns1,nhx)+ZE(ns2,nhx))
              THETA2=0.5d0*(ZE(ns3,nhx)+ZE(ns4,nhx))
              IF(THETA1.GT.THETA2+PI) THEN
                ZE(ns3,nhx)=ZE(ns3,nhx)+2.0d0*PI
                ZE(ns4,nhx)=ZE(ns4,nhx)+2.0d0*PI
                INCR12=.TRUE.
              ELSE IF(THETA2.GT.THETA1+PI) THEN
                ZE(ns1,nhx)=ZE(ns1,nhx)+2.0d0*PI
                ZE(ns2,nhx)=ZE(ns2,nhx)+2.0d0*PI
                INCR12=.TRUE.
              ENDIF
            ENDIF !coord system ITYP11(nr)

            IF(NITB.EQ.3) THEN !3D elements
              ns5=1+4*NKTB
              ns6=1+5*NKTB
              ns7=1+6*NKTB
              ns8=1+7*NKTB
              ns9=1+8*NKTB
              ns10=1+9*NKTB
              ns11=1+10*NKTB
              ns12=1+11*NKTB
C ***         Apply corrections for angular coordinates
              IF(ITYP11(nr).EQ.2.OR.ITYP11(nr).EQ.3) THEN !cyl or sph
                IF(ITYP11(nr).EQ.3) THEN !spherical
C                 Check for phi=90 (ie on axis)
                  IF(NPNE(5,nb).EQ.NPNE(6,nb) !5th & 6th in common
     '              .AND.DABS(ZE(ns5,3)-PI/2.0d0).LT.TOLERANCE !phi=90
     '              .AND.NPNE(7,nb).NE.NPNE(8,nb)) THEN
                    ZE(ns5,nhx)=ZE(ns7,nhx)
                    ZE(ns6,nhx)=ZE(ns8,nhx)
                  ELSE IF(NPNE(7,nb).EQ.NPNE(8,nb) !7th & 8th in common
     '                .AND.DABS(ZE(ns7,3)-PI/2.0d0).LT.TOLERANCE !phi=90
     '                .AND.NPNE(5,nb).NE.NPNE(6,nb)) THEN
                    ZE(ns7,nhx)=ZE(ns5,nhx)
                    ZE(ns8,nhx)=ZE(ns6,nhx)
                  ENDIF
                ENDIF
C               Make sure theta is increasing with Xi1
                IF(ZE(ns5,nhx).GE.ZE(ns6,nhx)-TOLERANCE) THEN
                  ZE(ns6,nhx)=ZE(ns6,nhx)+2.0d0*PI
                ENDIF
                IF(ZE(ns7,nhx).GE.ZE(ns8,nhx)-TOLERANCE) THEN
                  ZE(ns8,nhx)=ZE(ns8,nhx)+2.0d0*PI
                ENDIF
                IF(NNTB.EQ.12) THEN  !bilin/biCubicHermite-quadr. elem.
                  IF(ZE(ns9,nhx).GE.ZE(ns10,nhx)-TOLERANCE)
     '              ZE(ns10,nhx)=ZE(ns10,nhx)+2.0d0*PI
                  IF(ZE(ns11,nhx).GE.ZE(ns12,nhx)-TOLERANCE)
     '              ZE(ns12,nhx)=ZE(ns12,nhx)+2.0d0*PI
                ENDIF
              ELSE IF(ITYP11(nr).EQ.4.OR.ITYP11(nr).EQ.5) THEN !prol/obl
C               Make sure theta is decreasing with Xi1
                IF(ZE(ns6,nhx).GE.ZE(ns5,nhx)-TOLERANCE) THEN
                  INCR5=.TRUE. !5th vertex increased by 2*PI
                  ZE(ns5,nhx)=ZE(ns5,nhx)+2.0d0*PI
                ELSE
                  INCR5=.FALSE.
                ENDIF
                IF(ZE(ns8,nhx).GE.ZE(ns7,nhx)-TOLERANCE) THEN
                  INCR7=.TRUE. !7th vertex increased by 2*PI
                  ZE(ns7,nhx)=ZE(ns7,nhx)+2.0d0*PI
                ELSE
                  INCR7=.FALSE.
                ENDIF
                INCR56=.FALSE.
                IF(INCR5.AND..NOT.INCR7.AND.ZE(ns7,nhx).LT.PI) THEN
                  ZE(ns7,nhx)=ZE(ns7,nhx)+2.0d0*PI
                  ZE(ns8,nhx)=ZE(ns8,nhx)+2.0d0*PI
                  INCR56=.TRUE.
                ELSE IF(INCR7.AND..NOT.INCR5.AND.ZE(ns5,nhx).LT.PI)
     '            THEN
                  ZE(ns5,nhx)=ZE(ns5,nhx)+2.0d0*PI
                  ZE(ns6,nhx)=ZE(ns6,nhx)+2.0d0*PI
                  INCR56=.TRUE.
                ENDIF
C               Check whether top of element too different from bottom
                THETA1=0.5d0*(ZE(ns5,nhx)+ZE(ns6,nhx))
                THETA2=0.5d0*(ZE(ns7,nhx)+ZE(ns8,nhx))
                IF(THETA1.GT.THETA2+PI) THEN
                  ZE(ns7,nhx)=ZE(ns7,nhx)+2.0d0*PI
                  ZE(ns8,nhx)=ZE(ns8,nhx)+2.0d0*PI
                  INCR56=.TRUE.
                ELSE IF(THETA2.GT.THETA1+PI) THEN
                  ZE(ns5,nhx)=ZE(ns5,nhx)+2.0d0*PI
                  ZE(ns6,nhx)=ZE(ns6,nhx)+2.0d0*PI
                  INCR56=.TRUE.
                ENDIF
                IF(INCR12.AND..NOT.INCR56) THEN
                  ZE(ns5,nhx)=ZE(ns5,nhx)+2.0d0*PI
                  ZE(ns6,nhx)=ZE(ns6,nhx)+2.0d0*PI
                  ZE(ns7,nhx)=ZE(ns7,nhx)+2.0d0*PI
                  ZE(ns8,nhx)=ZE(ns8,nhx)+2.0d0*PI
                ELSE IF(INCR56.AND..NOT.INCR12) THEN
                  ZE(ns1,nhx)=ZE(ns1,nhx)+2.0d0*PI
                  ZE(ns2,nhx)=ZE(ns2,nhx)+2.0d0*PI
                  ZE(ns3,nhx)=ZE(ns3,nhx)+2.0d0*PI
                  ZE(ns4,nhx)=ZE(ns4,nhx)+2.0d0*PI
                ENDIF
C               Check whether outer vertices too different from inner
                THETA1=0.25d0*(ZE(ns1,nhx)+ZE(ns2,nhx)+ZE(ns3,nhx)
     '            +ZE(ns4,nhx))
                THETA2=0.25d0*(ZE(ns5,nhx)+ZE(ns6,nhx)+ZE(ns7,nhx)
     '            +ZE(ns8,nhx))
                IF(THETA1.GT.THETA2+PI) THEN
                  ZE(ns5,nhx)=ZE(ns5,nhx)+2.0d0*PI
                  ZE(ns6,nhx)=ZE(ns6,nhx)+2.0d0*PI
                  ZE(ns7,nhx)=ZE(ns7,nhx)+2.0d0*PI
                  ZE(ns8,nhx)=ZE(ns8,nhx)+2.0d0*PI
                ELSE IF(THETA2.GT.THETA1+PI) THEN
                  ZE(ns1,nhx)=ZE(ns1,nhx)+2.0d0*PI
                  ZE(ns2,nhx)=ZE(ns2,nhx)+2.0d0*PI
                  ZE(ns3,nhx)=ZE(ns3,nhx)+2.0d0*PI
                  ZE(ns4,nhx)=ZE(ns4,nhx)+2.0d0*PI
                ENDIF
              ENDIF !coord system ITYP11(nr)
            ENDIF !nitb=3

          ELSE IF(NNTB.EQ.9) THEN
            IF(ZE(1,nhx).GT.ZE(2,nhx)) ZE(1,nhx)=ZE(1,nhx)-2.0d0*PI
            IF(ZE(4,nhx).GT.ZE(5,nhx)) ZE(4,nhx)=ZE(4,nhx)-2.0d0*PI
            IF(ZE(7,nhx).GT.ZE(8,nhx)) ZE(7,nhx)=ZE(7,nhx)-2.0d0*PI
            IF(ZE(3,nhx).LT.ZE(2,nhx)) ZE(3,nhx)=ZE(3,nhx)+2.0d0*PI
            IF(ZE(6,nhx).LT.ZE(5,nhx)) ZE(6,nhx)=ZE(6,nhx)+2.0d0*PI
            IF(ZE(9,nhx).LT.ZE(8,nhx)) ZE(9,nhx)=ZE(9,nhx)+2.0d0*PI
            IF(NPNE(1,nb).EQ.NPNE(2,nb).AND.DABS(ZE(1,1)).LT.TOLERANCE
     '        .AND.NPNE(7,nb).NE.NPNE(9,nb)) THEN
              ZE(1,nhx)=ZE(7,nhx)
              ZE(2,nhx)=ZE(8,nhx)
              ZE(3,nhx)=ZE(9,nhx)
            ENDIF
          ENDIF
        ELSE IF(NITB.EQ.2.AND.NPF(1).EQ.3) THEN
          IF(ZE(ns1,nhx).GE.ZE(ns3,nhx)-TOLERANCE) THEN
            ZE(ns3,nhx)=ZE(ns3,nhx)+2.0d0*PI
          ENDIF
          IF(ZE(ns2,nhx).GE.ZE(ns4,nhx)-TOLERANCE) THEN
            ZE(ns4,nhx)=ZE(ns4,nhx)+2.0d0*PI
          ENDIF
        ENDIF !nitb

      ENDIF !curvilinear coords etc

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        DO nhx=1,NHE
          nh=NH_LOC(nhx,nx)
          WRITE(OP_STRING,
     '      '('' ZE(ns,'',I2,''): '',8D11.3,/(12X,8D11.3))')
     '      nhx,(ZE(ns,nhx),ns=1,NST(NBH(nh,nc))+NAT(NBH(nh,nc)))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !nhx
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('ZPZE')
      RETURN
 9999 CALL ERRORS('ZPZE',ERROR)
      CALL EXITS('ZPZE')
      RETURN 1
      END


