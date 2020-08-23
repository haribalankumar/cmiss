      SUBROUTINE CHMESH_ALTER(IBT,IDO,INP,NBJ,NEELEM,NEL,NENP,
     '  NP_INTERFACE,NPL,NPNE,NPNODE,nr,NVJL,DL,SE,TOTAL,XP,ERROR,*)

C#### Subroutine: CHMESH_ALTER
C###  Description:
C###    CHMESH_ALTER is J.Crocombes alter subroutine.  Alters the mesh
C###    values.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'chmesh0.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NEL(0:NELM,NLM),
     '  NENP(NPM,0:NEPM,0:NRM),NP_INTERFACE(0:NPM,0:3),
     '  NPL(5,0:3,NLM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),nr,
     '  NVJL(4,NJM,NLM)
      CHARACTER ERROR*(*)
      REAL*8 DL(3,NLM),SE(NSM,NBFM,NEM),TOTAL,XP(NKM,NVM,NJM,NPM)
!     Local Variables
      REAL*8 ARC,CHMESH_DFDS,CHMESH_DGDS,CHMESH_FFUNC,
     '  CHMESH_GFUNC,DEPTH,ETA,LENGTH,MIN,MAX,TEMP,TEMP1,TEMP2,THETA,
     '  TRANS(3,4),WIDTH,XAVG,XX(3),YAVG
      INTEGER i,j,k,nb,nj,nk,nl,noline,nonode,np,nt,nu

      CALL ENTERS('CHMESH_ALTER',*9999)

      !Sorts length meas into order of incr. values if variable scal
      IF(SCLTYPE.EQ.2.AND.NJT.GT.2) THEN
        j=NTL
        DO i=1,NTL
          k=j-1
          DO j=1,k
            IF(TORSO_LENGTHS(1,j).GT.TORSO_LENGTHS(1,j+1))THEN
              TEMP1=TORSO_LENGTHS(1,j)
              TEMP2=TORSO_LENGTHS(2,j)
              TORSO_LENGTHS(1,j)=TORSO_LENGTHS(1,j+1)
              TORSO_LENGTHS(2,j)=TORSO_LENGTHS(2,j+1)
              TORSO_LENGTHS(1,j+1)=TEMP1
              TORSO_LENGTHS(2,j+1)=TEMP2
            ENDIF
          ENDDO !j
        ENDDO !i
C AJPs 2 July 1998
C Old code that did not work properly on the generic pig.
C Also did not understand the nu loop.
C        DO nt=1,NTL-1
C          DO nonode=1,NPNODE(0,nr)
C            np=NPNODE(nonode,nr)
C            IF((XP(1,1,3,np).GE.TORSO_LENGTHS(1,nt)).AND.(XP(1,1,3,np)
C     '        .LT.TORSO_LENGTHS(1,nt+1))) THEN
C              DO nu=1,3
C              XP(nu,1,3,np)=XP(nu,1,3,np)*((TORSO_LENGTHS(2,nt)/
C     '          TORSO_LENGTHS(1,nt))*(1-((XP(1,1,3,np)-
C     '          TORSO_LENGTHS(1,nt))/(TORSO_LENGTHS(1,nt+1)-
C     '          TORSO_LENGTHS(1,nt))))+((TORSO_LENGTHS(2,nt+1)/
C     '          TORSO_LENGTHS(1,nt+1))*((XP(1,1,3,np)-
C     '          TORSO_LENGTHS(1,nt))/(TORSO_LENGTHS(1,nt+1)-
C     '          TORSO_LENGTHS(1,nt)))))
C              ENDDO
C            ELSEIF(XP(1,1,3,np).GE.TORSO_LENGTHS(1,NTL).AND.nt.EQ.NTL-1)
C     '        THEN
C              DO nu=1,3
C                XP(nu,1,3,np)=XP(nu,1,3,np)*(TORSO_LENGTHS(2,NTL)/
C     '          TORSO_LENGTHS(1,NTL))
C              ENDDO
C            ELSEIF(XP(1,1,3,np).LT.TORSO_LENGTHS(1,1).AND.nt.EQ.1) THEN
C              DO nu=1,3
C                XP(nu,1,3,np)=XP(nu,1,3,np)*(TORSO_LENGTHS(2,nt)/
C     '          TORSO_LENGTHS(1,nt))
C              ENDDO
C            ENDIF
C          ENDDO !np
C        ENDDO !nt
C New code that allows for a vertical translation as well (and may be
C more efficient since it only loops over nodes once).
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          IF(nr.EQ.NP_INTERFACE(np,1)) THEN
            nt=1
            DO WHILE((nt.LT.NTL).AND.
     '        (XP(1,1,3,np).GE.TORSO_LENGTHS(1,nt)))
              nt=nt+1
            ENDDO !nt
            !Current z elevation is below TORSO_LENGTH(1,nt) OR above
             !top landmark.
            IF((nt.GT.1).AND.(nt.LT.NTL)) THEN
              !Point is between torso_length( ,nt-1) and
              !torso_length( ,nt).
              XP(1,1,3,np)=(XP(1,1,3,np)-TORSO_LENGTHS(1,nt-1))*
     '          (TORSO_LENGTHS(2,nt)-TORSO_LENGTHS(2,nt-1))/
     '          (TORSO_LENGTHS(1,nt)-TORSO_LENGTHS(1,nt-1))+
     '          TORSO_LENGTHS(2,nt-1)
            ELSEIF(nt.EQ.1) THEN !Point below bottom marker
              XP(1,1,3,np)=(XP(1,1,3,np)-TORSO_LENGTHS(1,1))*
     '          (TORSO_LENGTHS(2,2)-TORSO_LENGTHS(2,1))/
     '          (TORSO_LENGTHS(1,2)-TORSO_LENGTHS(1,1))+
     '          TORSO_LENGTHS(2,1)
            ELSEIF(nt.EQ.NTL) THEN !Point in final interval or above top
              XP(1,1,3,np)=(XP(1,1,3,np)-TORSO_LENGTHS(1,NTL-1))*
     '          (TORSO_LENGTHS(2,NTL)-TORSO_LENGTHS(2,NTL-1))/
     '          (TORSO_LENGTHS(1,NTL)-TORSO_LENGTHS(1,NTL-1))+
     '          TORSO_LENGTHS(2,NTL-1)
            ELSE
              ERROR='>>Should not get here - no interval found??'
              GOTO 9999
            ENDIF !nt
          ENDIF !np_interface
        ENDDO !np
      ENDIF

      !Calculates the total length of the mesh.
      IF(NJT.GT.2) THEN
        MAX = XP(1,1,3,1)
        MIN = XP(1,1,3,1)
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          IF(XP(1,1,3,np).GT.MAX) THEN
            MAX = XP(1,1,3,np)
          ENDIF
          IF(XP(1,1,3,np).LT.MIN) THEN
            MIN = XP(1,1,3,np)
          ENDIF
        ENDDO
        LENGTH = MAX - MIN

        IF(SCLTYPE.EQ.1) THEN
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            DO nu=1,3 !???WHat is this loop for AJP 2 July 1998
              XP(nu,1,3,np)=XP(nu,1,3,np)*(TORSO_LENGTHS(0,1)/LENGTH)
            ENDDO
          ENDDO
        ENDIF
      ELSEIF(NJT.EQ.2) THEN
        ETA=0.0d0
        LENGTH=0.0d0
      ENDIF

      IF(CUSTOMISATION_TYPE.NE.6) THEN
      !Alters the coordinate and derivative values using a polynomial
      !function in the z direction and cosine functions on the x,y plane
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
C GMH 8/1/97 Update cmgui link
          CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
          !Sets any very small values to be very small i.e effectively
              !zero but positive so the correct theta value is found.
          DO nj = 1,2
            IF(DABS(XP(1,1,nj,np)).LT.ZERO_TOL)THEN
              XP(1,1,nj,np) = 0.0d0
            ENDIF
          ENDDO

        !Finds values of eta, the normalised vertical distance, and
        !theta ,the angle the radius at the point makes with the
        !positve x axis.
          IF(NJT.GT.2) THEN
            ETA = XP(1,1,3,np)/LENGTH
            IF(DABS(XP(1,1,1,np)).GT.ZERO_TOL)THEN
              THETA = DATAN(XP(1,1,2,np)/XP(1,1,1,np))
              IF(XP(1,1,1,np).LT.-ZERO_TOL)THEN
                THETA = THETA + PI
              ELSE IF (XP(1,1,2,np).LT.-ZERO_TOL) THEN
                THETA = THETA + 2*PI
              ENDIF
            ELSE IF(XP(1,1,2,np).GT.ZERO_TOL) THEN
              THETA = PI*0.5d0
            ELSE
              THETA = PI*0.5d0+PI
            ENDIF
          ELSEIF(NJT.EQ.2)THEN !theta is arclength
            IF(np.EQ.NPNODE(1,nr)) THEN
              THETA=0.0d0
              ARC=0.0d0
            ELSE
              DO noline=1,NENP(np,0,nr)
                nl=NENP(np,noline,nr)
                IF(np.EQ.NPNE(2,NBJ(2,nl),nl)) THEN
                  ARC=ARC+DL(3,nl)
                  THETA=(ARC/TOTAL)*2*PI
                ENDIF
              ENDDO
            ENDIF
          ENDIF

          !Calculates the x & y derivitives with respect to s1 & s2.
          IF((DABS(XP(1,1,1,np)).GT.ZERO_TOL).OR.
     '      (DABS(XP(1,1,2,np)).GT.ZERO_TOL))THEN
            DO nk = 2,NJT
              TEMP = XP(nk,1,1,np)
              XP(nk,1,1,np) =
     '        XP(nk,1,1,np)*CHMESH_FFUNC(ETA)*
     '        CHMESH_GFUNC(THETA)
     '        +XP(1,1,1,np)*CHMESH_GFUNC(THETA)
     '        *CHMESH_DFDS(ETA,LENGTH,XP,np,nk)+XP(1,1,1,np)
     '        *CHMESH_FFUNC(ETA)*
     '         CHMESH_DGDS(THETA,XP,np,nk,TEMP)

              XP(nk,1,2,np) =
     '        XP(nk,1,2,np)*CHMESH_FFUNC(ETA)*
     '        CHMESH_GFUNC(THETA)
     '        +XP(1,1,2,np)*CHMESH_GFUNC(THETA)
     '        *CHMESH_DFDS(ETA,LENGTH,XP,np,nk)+XP(1,1,2,np)
     '        *CHMESH_FFUNC(ETA)*
     '        CHMESH_DGDS(THETA,XP,np,nk,TEMP)
            ENDDO
          ENDIF
          !Multiplies the old x and y coordinate values by f(zeta) and
          !g(theta) to give a new shape.
          DO nj = 1,2
            XP(1,1,nj,np) = CHMESH_FFUNC(ETA)*
     '        CHMESH_GFUNC(THETA)*XP(1,1,nj,np)
          ENDDO

          !Normalises the derivative vectors.

          IF(NJT.GT.2)  CALL CHMESH_NORMALISE(XP,3,np,ERROR,*9999)
          CALL CHMESH_NORMALISE(XP,2,np,ERROR,*9999)
C           CALL NORMALISE(3,XP(1,3,nj,np),ERROR,*9999)
C           CALL NORMALISE(3,XP(1,2,nj,np),ERROR,*9999)

        ENDDO !nonode
      ELSEIF(CUSTOMISATION_TYPE.EQ.6) THEN
        XAVG=0.0d0
        YAVG=0.0d0
        DO nt=1,NTDM
          CALL MESHXY(IBT,IDO,INP,NBJ,nr,NEELEM,NPNE,DEPTH,SE,
     '             XP,WIDTH,WDMEASURE(3,nt),ERROR,*9999)
          XAVG=(WDMEASURE(1,nt)/WIDTH)+XAVG
          YAVG=(WDMEASURE(2,nt)/DEPTH)+YAVG
        ENDDO
        XAVG=XAVG/NTDM
        YAVG=YAVG/NTDM
        CALL RESET(TRANS)
        CALL SCALE1(XAVG,TRANS)
        CALL SCALE2(YAVG,TRANS)

        XX(1)=0.d0
        XX(2)=0.d0
        XX(3)=0.d0
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
C GMH 8/1/97 Update cmgui link
          CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
C             transform nodal values
          DO nj=1,NJT
            XX(nj)=XP(1,1,nj,np)
          ENDDO
          CALL ZZ(XX,XX,TRANS)
          DO nj=1,NJT
            XP(1,1,nj,np)=XX(nj)
          ENDDO
        ENDDO
        DO nb=1,NBFT
          IF(NBI(nb).EQ.3.OR.(NBI(nb).GE.5.AND.NBI(nb).LE.7)) THEN
            DO nl=1,NLT
              CALL ARCSCA(IDO,0,0,0,NBJ,NEL(0,nl),nl,
     '          NPL(1,0,nl),NPNE,NVJL(1,1,nl),DL,1.0D-6,XP,ERROR,*9999)
            ENDDO
          ENDIF
        ENDDO

      ENDIF

      CALL EXITS('CHMESH_ALTER')
      RETURN
 9999 CALL ERRORS('CHMESH_ALTER',ERROR)
      CALL EXITS('CHMESH_ALTER')
      RETURN 1
      END


