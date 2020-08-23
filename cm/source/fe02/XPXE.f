      SUBROUTINE XPXE(NBJ,NKJE,NPF,NPNE,nr,NVJE,SE,XA,XE,XP,ERROR,*)

C#### Subroutine: XPXE
C###  Description:
C###    XPXE transfers global parameters XP and XA to element parameters
C###    XE (jtyp5=1) or monomial coefficients (jtyp5=2).

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBJ(NJM),NKJE(NKM,NNM,NJM),NPF(9),NPNE(NNM,NBFM),
     '  nr,NVJE(NNM,NBFM,NJM)
      REAL*8 SE(NSM,NBFM),XA(NAM,NJM),XE(NSM,NJM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER mk,ms,na,nb,NITB,nj,njj1,njj2,NJJSTEP,nk,NK1,
     '  NKTB,nn,NNTB,np,
     '  ns,ns1,ns2,ns3,ns4,ns5,ns6,ns7,ns8,ns9,ns10,ns11,ns12,NSTB,nv
      REAL*8 A(8,8),SUM,THETA1,THETA2,TOLERANCE,XM(8),Y,YC,YS
      LOGICAL INCR1,INCR12,INCR3,INCR5,INCR56,INCR7

      DATA A/1.0D0,-1.0D0,-1.0D0,1.0D0,60*0.0D0/
      DATA TOLERANCE/1.0D-10/

      CALL ENTERS('XPXE',*9999)

      IF(KTYP8.EQ.1.OR.KTYP8.EQ.5) THEN
        NJJSTEP=2               !Loop over geometry and field
      ELSE IF(KTYP8.EQ.2.OR.KTYP8.EQ.3) THEN
        NJJSTEP=1               !Loop over geometry, fibres and field
      ELSE IF(KTYP8.EQ.6) THEN
        NJJSTEP=1               !Loop over geometry
      ELSE
        NJJSTEP=1               !Loop over geometry, fibres and field
      ENDIF
      DO njj1=1,3,NJJSTEP
        DO njj2=1,NJ_LOC(njj1,0,nr)
          nj=NJ_LOC(njj1,njj2,nr)
          nb=NBJ(nj)
          IF(nb.GT.0) THEN
            ns=0
            IF(NNT(nb).GT.0) THEN
              DO nn=1,NNT(nb)
                np=NPNE(nn,nb)
                nv=NVJE(nn,nb,nj)
                DO mk=1,NKT(nn,nb)
                  ns=ns+1
C KAT 2Aug98: IF not necesary
C                  IF(JTYP2.EQ.0.OR.JTYP2.EQ.2.OR.JTYP2.EQ.3) THEN
                  nk=NKJE(mk,nn,nj)
C                  ENDIF
                  IF(NBI(nb).EQ.4) THEN !Scales calc.d from angle change
                    XE(ns,nj)=XP(nk,nv,nj,np)*SE(ns,nb) !needs addition
                  ELSE
                    XE(ns,nj)=XP(nk,nv,nj,np)*SE(ns,nb)
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
c            DO na=1,NAT(nb) !for each auxilliary parameter
c              XE(ns+na,nj)=XA(na,nj)
c            ENDDO !na
          ENDIF !nb.GT.0
        ENDDO !njj2
      ENDDO !njj1

C     IF(ITYP10(nr).GT.1.AND.NIT(NBJ(1)).GT.1) THEN
      IF(ITYP10(nr).GT.1) THEN
        IF(JTYP10.GE.2) THEN
          nb=NBJ(1)
          NK1=NKT(0,nb)
          IF(NK1.GT.1) THEN
C****       18-Feb-89: Further modification needed if nodal parameters
C****       include second derivs, in which case we need IDO(nk,nn,0,nb)
            DO nn=1,NNT(nb)
              ns1=1+(nn-1)*NK1
              DO nk=2,NK1
                ns2=nk+(nn-1)*NK1
                IF(ITYP10(nr).LT.4) THEN
                  XE(ns2,1)=ITYP10(nr)*XE(ns1,1)**(ITYP10(nr)-1)
     '                                *XE(ns2,1)
                ELSE IF(ITYP10(nr).EQ.4) THEN
                  YC=DCOSH(XE(ns1,1))
                  YS=DSINH(XE(ns1,1))
                  IF(JTYP10.EQ.2) THEN
                    XE(ns2,1)=2.0D0*FOCUS*FOCUS*YC*YS*XE(ns2,1)
                  ELSE IF(JTYP10.EQ.3) THEN
                    XE(ns2,1)=FOCUS**3*YS*(3.0D0*YC*YC-1.0D0)*XE(ns2,1)
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDIF
          DO nn=1,NNT(nb)
            ns=1+(nn-1)*NK1
            IF(ITYP10(nr).LT.4) THEN
              XE(ns,1)=XE(ns,1)**ITYP10(nr)
            ELSE IF(ITYP10(nr).EQ.4) THEN
              Y=DCOSH(XE(ns,1))
              IF(JTYP10.EQ.2) THEN
                XE(ns,1)=FOCUS*FOCUS*(Y*Y-1.0D0)
              ELSE IF(JTYP10.EQ.3) THEN
                XE(ns,1)=FOCUS**3*Y*(Y*Y-1.0D0)
              ENDIF
            ELSE IF(ITYP10(nr).EQ.5) THEN
            ENDIF
          ENDDO
        ENDIF

        IF(ITYP10(nr).LE.3) THEN      !cyl. or sph. polar
          nj=2   !theta
        ELSE IF(ITYP10(nr).EQ.4) THEN !prolate
          nj=3   !theta
        ENDIF
        nb=NBJ(nj)
        NNTB=NNT(nb)
        NITB=NIT(nb)
        NKTB=NKT(0,nb)
C!!! KAT 2001-06-26: No check for NNTB >= 4 here
        ns1=1
        ns2=1+  NKTB
        ns3=1+2*NKTB
        ns4=1+3*NKTB
C !!! KAT 2001-06-26: If NITB is 1 then we don't know what xi direction
C       this really corresponds to as this could be any edge of a 3d
C       element.  There seems to be a lot of using nj's and ns's without
C       checking that they are defined.
        IF(NITB.EQ.1.OR.NITB.EQ.3.OR.(NITB.EQ.2.AND.NPF(1).EQ.1)) THEN
          IF(NNTB.NE.9) THEN
C ***       Apply corrections for angular coordinates
            IF(ITYP10(nr).EQ.2.OR.ITYP10(nr).EQ.3) THEN !cyl or sph
              IF(ITYP10(nr).EQ.3) THEN !spherical
C               Check for phi=90 (ie on axis)
                IF(NPNE(1,nb).EQ.NPNE(2,nb) !1st & 2nd nodes common
     '            .AND.DABS(XE(ns1,3)-PI/2.0D0).LT.TOLERANCE !phi=90
     '            .AND.NPNE(3,nb).NE.NPNE(4,nb)) THEN
                  XE(ns1,nj)=XE(ns3,nj)
                  XE(ns2,nj)=XE(ns4,nj)
                ELSE IF(NPNE(3,nb).EQ.NPNE(4,nb) !3rd & 4th nodes common
     '              .AND.DABS(XE(ns3,3)-PI/2.0D0).LT.TOLERANCE !phi=90
     '              .AND.NPNE(1,nb).NE.NPNE(2,nb)) THEN
                  XE(ns3,nj)=XE(ns1,nj)
                  XE(ns4,nj)=XE(ns2,nj)
                ENDIF
              ENDIF
C             Make sure theta is increasing with Xi1
              IF(XE(ns1,nj).GE.XE(ns2,nj)-TOLERANCE)
     '          XE(ns2,nj)=XE(ns2,nj)+2.0D0*PI
              IF(NITB.GE.2) THEN
                IF(XE(ns3,nj).GE.XE(ns4,nj)-TOLERANCE)
     '            XE(ns4,nj)=XE(ns4,nj)+2.0D0*PI
              ENDIF
            ELSE IF(ITYP10(nr).EQ.4.OR.ITYP10(nr).EQ.5)THEN !prol or obl
C             Make sure theta is decreasing with Xi1
              IF(XE(ns2,nj).GE.XE(ns1,nj)-TOLERANCE) THEN
                INCR1=.TRUE. !1st vertex increased by 2*PI
                XE(ns1,nj)=XE(ns1,nj)+2.0D0*PI
              ELSE
                INCR1=.FALSE.
              ENDIF
              IF(XE(ns4,nj).GE.XE(ns3,nj)-TOLERANCE) THEN
                INCR3=.TRUE. !3rd vertex increased by 2*PI
                XE(ns3,nj)=XE(ns3,nj)+2.0D0*PI
              ELSE
                INCR3=.FALSE.
              ENDIF
              INCR12=.FALSE.
              IF(INCR1.AND..NOT.INCR3.AND.XE(ns3,nj).LT.PI) THEN
                XE(ns3,nj)=XE(ns3,nj)+2.0D0*PI
                XE(ns4,nj)=XE(ns4,nj)+2.0D0*PI
                INCR12=.TRUE.
              ELSE IF(INCR3.AND..NOT.INCR1.AND.XE(ns1,nj).LT.PI) THEN
                XE(ns1,nj)=XE(ns1,nj)+2.0D0*PI
                XE(ns2,nj)=XE(ns2,nj)+2.0D0*PI
                INCR12=.TRUE.
              ENDIF
C             Check whether top of element too different from bottom
              THETA1=0.5D0*(XE(ns1,nj)+XE(ns2,nj))
              THETA2=0.5D0*(XE(ns3,nj)+XE(ns4,nj))
              IF(THETA1.GT.THETA2+PI) THEN
                XE(ns3,nj)=XE(ns3,nj)+2.0D0*PI
                XE(ns4,nj)=XE(ns4,nj)+2.0D0*PI
                INCR12=.TRUE.
              ELSE IF(THETA2.GT.THETA1+PI) THEN
                XE(ns1,nj)=XE(ns1,nj)+2.0D0*PI
                XE(ns2,nj)=XE(ns2,nj)+2.0D0*PI
                INCR12=.TRUE.
              ENDIF
            ENDIF !coord system ITYP10(nr)

            IF(NITB.EQ.3) THEN !3D elements
              ns5=1+4*NKTB
              ns6=1+5*NKTB
              ns7=1+6*NKTB
              ns8=1+7*NKTB
              ns9=1+8*NKTB
              ns10=1+9*NKTB
              ns11=1+10*NKTB
              ns12=1+11*NKTB
!MPN 17-Nov-94:correct version (nv) of the njth coord of global node
!              is picked up via NVJP,NVJE
!              which are set up in IPNODE,IPELEM resp.
!old
c!             Check for zero radius (cyl or sph) or mu (prolate)
c!             to set theta
c              IF(NPNE(5,nb).EQ.NPNE(6,nb)      !5th & 6th nodes common
c     '          .AND.DABS(XE(ns5,njj)).LT.TOLERANCE  !radius=0
c     '          .AND.NPNE(7,nb).NE.NPNE(8,nb)) THEN
c                XE(ns5,nj)=XE(ns7,nj)
c                XE(ns6,nj)=XE(ns8,nj)
c              ELSE IF(NPNE(7,nb).EQ.NPNE(8,nb) !7th & 8th nodes common
c     '          .AND.DABS(XE(ns7,njj)).LT.TOLERANCE !radius=0
c     '          .AND.NPNE(5,nb).NE.NPNE(6,nb)) THEN
c                XE(ns7,nj)=XE(ns5,nj)
c                XE(ns8,nj)=XE(ns6,nj)
c              ENDIF
!endold

C ***         Apply corrections for angular coordinates
              IF(ITYP10(nr).EQ.2.OR.ITYP10(nr).EQ.3) THEN !cyl or sph
                IF(ITYP10(nr).EQ.3) THEN !spherical
C                 Check for phi=90 (ie on axis)
                  IF(NPNE(5,nb).EQ.NPNE(6,nb)     !5th/6th nodes common
     '              .AND.DABS(XE(ns5,3)-PI/2.0D0).LT.TOLERANCE !phi=90
     '              .AND.NPNE(7,nb).NE.NPNE(8,nb)) THEN
                    XE(ns5,nj)=XE(ns7,nj)
                    XE(ns6,nj)=XE(ns8,nj)
                  ELSE IF(NPNE(7,nb).EQ.NPNE(8,nb) !7th/8th nodes common
     '                .AND.DABS(XE(ns7,3)-PI/2.0D0).LT.TOLERANCE !phi=90
     '                .AND.NPNE(5,nb).NE.NPNE(6,nb)) THEN
                    XE(ns7,nj)=XE(ns5,nj)
                    XE(ns8,nj)=XE(ns6,nj)
                  ENDIF
                ENDIF
C               Make sure theta is increasing with Xi1
                IF(XE(ns5,nj).GE.XE(ns6,nj)-TOLERANCE)
     '            XE(ns6,nj)=XE(ns6,nj)+2.0D0*PI
                IF(XE(ns7,nj).GE.XE(ns8,nj)-TOLERANCE)
     '            XE(ns8,nj)=XE(ns8,nj)+2.0D0*PI
                IF(NNTB.EQ.12) THEN  !bilin/biCubicHermite-quadr. elem.
                  IF(XE(ns9,nj).GE.XE(ns10,nj)-TOLERANCE)
     '              XE(ns10,nj)=XE(ns10,nj)+2.0D0*PI
                  IF(XE(ns11,nj).GE.XE(ns12,nj)-TOLERANCE)
     '              XE(ns12,nj)=XE(ns12,nj)+2.0D0*PI
                ENDIF
              ELSE IF(ITYP10(nr).EQ.4.OR.ITYP10(nr).EQ.5) THEN !prol/obl
C               Make sure theta is decreasing with Xi1
                IF(XE(ns6,nj).GE.XE(ns5,nj)-TOLERANCE) THEN
                  INCR5=.TRUE. !5th vertex increased by 2*PI
                  XE(ns5,nj)=XE(ns5,nj)+2.0D0*PI
                ELSE
                  INCR5=.FALSE.
                ENDIF
                IF(XE(ns8,nj).GE.XE(ns7,nj)-TOLERANCE) THEN
                  INCR7=.TRUE. !7th vertex increased by 2*PI
                  XE(ns7,nj)=XE(ns7,nj)+2.0D0*PI
                ELSE
                  INCR7=.FALSE.
                ENDIF
                INCR56=.FALSE.
                IF(INCR5.AND..NOT.INCR7.AND.XE(ns7,nj).LT.PI) THEN
                  XE(ns7,nj)=XE(ns7,nj)+2.0D0*PI
                  XE(ns8,nj)=XE(ns8,nj)+2.0D0*PI
                  INCR56=.TRUE.
                ELSE IF(INCR7.AND..NOT.INCR5.AND.XE(ns5,nj).LT.PI) THEN
                  XE(ns5,nj)=XE(ns5,nj)+2.0D0*PI
                  XE(ns6,nj)=XE(ns6,nj)+2.0D0*PI
                  INCR56=.TRUE.
                ENDIF
C               Check whether top of element too different from bottom
                THETA1=0.5D0*(XE(ns5,nj)+XE(ns6,nj))
                THETA2=0.5D0*(XE(ns7,nj)+XE(ns8,nj))
                IF(THETA1.GT.THETA2+PI) THEN
                  XE(ns7,nj)=XE(ns7,nj)+2.0D0*PI
                  XE(ns8,nj)=XE(ns8,nj)+2.0D0*PI
                  INCR56=.TRUE.
                ELSE IF(THETA2.GT.THETA1+PI) THEN
                  XE(ns5,nj)=XE(ns5,nj)+2.0D0*PI
                  XE(ns6,nj)=XE(ns6,nj)+2.0D0*PI
                  INCR56=.TRUE.
                ENDIF
                IF(INCR12.AND..NOT.INCR56) THEN
                  XE(ns5,nj)=XE(ns5,nj)+2.0D0*PI
                  XE(ns6,nj)=XE(ns6,nj)+2.0D0*PI
                  XE(ns7,nj)=XE(ns7,nj)+2.0D0*PI
                  XE(ns8,nj)=XE(ns8,nj)+2.0D0*PI
                ELSE IF(INCR56.AND..NOT.INCR12) THEN
                  XE(ns1,nj)=XE(ns1,nj)+2.0D0*PI
                  XE(ns2,nj)=XE(ns2,nj)+2.0D0*PI
                  XE(ns3,nj)=XE(ns3,nj)+2.0D0*PI
                  XE(ns4,nj)=XE(ns4,nj)+2.0D0*PI
                ENDIF
C               Check whether outer vertices too different from inner
                THETA1=0.25D0*(XE(ns1,nj)+XE(ns2,nj)+XE(ns3,nj)
     '            +XE(ns4,nj))
                THETA2=0.25D0*(XE(ns5,nj)+XE(ns6,nj)+XE(ns7,nj)
     '            +XE(ns8,nj))
                IF(THETA1.GT.THETA2+PI) THEN
                  XE(ns5,nj)=XE(ns5,nj)+2.0D0*PI
                  XE(ns6,nj)=XE(ns6,nj)+2.0D0*PI
                  XE(ns7,nj)=XE(ns7,nj)+2.0D0*PI
                  XE(ns8,nj)=XE(ns8,nj)+2.0D0*PI
                ELSE IF(THETA2.GT.THETA1+PI) THEN
                  XE(ns1,nj)=XE(ns1,nj)+2.0D0*PI
                  XE(ns2,nj)=XE(ns2,nj)+2.0D0*PI
                  XE(ns3,nj)=XE(ns3,nj)+2.0D0*PI
                  XE(ns4,nj)=XE(ns4,nj)+2.0D0*PI
                ENDIF
              ENDIF !coord system ITYP10(nr)
            ENDIF !nitb=3

          ELSE IF(NNTB.EQ.9) THEN
            IF(XE(1,nj).GT.XE(2,nj)) XE(1,nj)=XE(1,nj)-2.0D0*PI
            IF(XE(4,nj).GT.XE(5,nj)) XE(4,nj)=XE(4,nj)-2.0D0*PI
            IF(XE(7,nj).GT.XE(8,nj)) XE(7,nj)=XE(7,nj)-2.0D0*PI
            IF(XE(3,nj).LT.XE(2,nj)) XE(3,nj)=XE(3,nj)+2.0D0*PI
            IF(XE(6,nj).LT.XE(5,nj)) XE(6,nj)=XE(6,nj)+2.0D0*PI
            IF(XE(9,nj).LT.XE(8,nj)) XE(9,nj)=XE(9,nj)+2.0D0*PI
            IF(NPNE(1,nb).EQ.NPNE(2,nb).AND.DABS(XE(1,1)).LT.TOLERANCE
     '        .AND.NPNE(7,nb).NE.NPNE(9,nb)) THEN
              XE(1,nj)=XE(7,nj)
              XE(2,nj)=XE(8,nj)
              XE(3,nj)=XE(9,nj)
            ENDIF
          ENDIF

        ELSE IF(NITB.EQ.2.AND.NPF(1).EQ.3) THEN
          IF(XE(ns1,nj).GE.XE(ns3,nj)-TOLERANCE)
     '      XE(ns3,nj)=XE(ns3,nj)+2.0D0*PI
          IF(XE(ns2,nj).GE.XE(ns4,nj)-TOLERANCE)
     '      XE(ns4,nj)=XE(ns4,nj)+2.0D0*PI
        ENDIF

! PJH 4Aug93
!       IF(ITYP10(nr).EQ.3) THEN !spherical polar
!         nb=NBJ(3)
!         nk2=NKT(0,nb)
!         DO nn=1,2
!           ns =1+(nn-1)*nk2
!           ns1=1+(nn+1)*nk2
!           IF(XE(ns,3).GE.XE(ns1,3)-TOLERANCE)
!     '       XE(ns1,3)=XE(ns1,3)+2.0D0*PI
!           IF(NNT(nb).GT.4) THEN
!             ns4=1+(nn+3)*nk2
!             ns5=1+(nn+5)*nk2
!             IF(XE(ns4,3).GE.XE(ns5,3)-TOLERANCE)
!     '         XE(ns5,3)=XE(ns5,3)+2.0D0*PI
!           ENDIF
!         ENDDO
!       ENDIF
      ENDIF

      IF(JTYP5.EQ.2) THEN
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
C           DO nj=1,NJTOT
            nb=NBJ(nj)
            NSTB=NST(nb)
            DO ms=1,NSTB
              SUM=0.0D0
              DO ns=1,NSTB
                SUM=SUM+A(ms,ns)*XE(ns,nj)
              ENDDO
              XM(ms)=SUM
            ENDDO
            DO ms=1,NSTB
              XE(ms,nj)=XM(ms)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

C     IF(JTYP9.EQ.1) THEN
C
C       nb=NBJ(nj)
C       DO nn=1,NNT(nb)
C         ns=1+(nn-1)*NKT(0,nb)
C         IF(XE(ns,nj).LT.0.0) THEN
C           XE(ns,nj)=XE(ns,nj)+PI
C         ELSE IF(XE(ns,nj).GT.PI) THEN
C           XE(ns,nj)=XE(ns,nj)-PI
C         ENDIF
C       ENDDO
C     ENDIF

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        DO njj1=1,3,NJJSTEP
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            IF(NBJ(nj).GT.0) THEN
              WRITE(OP_STRING,'('' XE(ns,'',I2,''): '',4(1X,D12.4),'
     '          //'/(12X,4(1X,D12.4)))')
     '          nj,(XE(ns,nj),ns=1,NST(NBJ(nj)))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO
        ENDDO
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('XPXE')
      RETURN
 9999 CALL ERRORS('XPXE',ERROR)
      CALL EXITS('XPXE')
      RETURN 1
      END


