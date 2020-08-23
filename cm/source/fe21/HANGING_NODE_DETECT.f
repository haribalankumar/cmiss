      SUBROUTINE HANGING_NODE_DETECT(IBT,IDO,IDRN,
     '  INP,NBJ,NEELEM,NELIST,
     '  NENP,NKJE,NPF,NPLIST,NPNE,NPNODE,
     '  nr,NVJE,NVJP,NWP,SE,SP,XA,XE,XP,ERROR,*)

C#### Subroutine: HANGING_NODE_DETECT
C###  Description:
C###    <HTML>
C###    <PRE>
C###    HANGING_NODE_DETECT detects hanging nodes and flags them in NWP
C###    NWP(np,1)=element node hangs in
C###    NWP(np,2)=local line/face of element that node hangs on
C###    </PRE>
C###    </HTML>
C**** Created by Carey Stevens 13 May 1998

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'mxch.inc'
      INCLUDE 'parameters.inc'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  IDRN,INP(NNM,NIM,NBFM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NENP(NPM,0:NEPM,0:NRM),NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),
     '  NPLIST(0:NPM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  nr,NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),NWP(NPM,2)
      REAL*8 SE(NSM,NBFM,NEM),SP(NKM,NBFM,NPM),XA(NAM,NJM,NEM),
     '  XE(NSM,NJM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER LD,nb,nb2,nb3,ne,ni,NITB,NITE,nj,nj1,nk,
     '  nk2,NKJNEED,nn,noelem,
     '  no_nelist,nonode,np,ns,nu,NUNK(8),pos_inlist
      REAL*8 AA,DSDXI(3),G1,G3,PXI,R,RC,RR,RRC,
     '  SCALEFACTOR(3),SLX,SMX,SUM,X(11,20),XD(3),XI(3),XS(8),XTEMP
      LOGICAL  CALCDSDXI,FOUND,FOUND_LINEAR_BASIS,INLIST,NOCROSS
C      LOGICAL SECTOR

      DATA NUNK/1,2,4,6,7,9,10,11/

      CALL ENTERS('HANGING_NODE_DETECT',*9999)

C     Check for any hanging nodes
C     If the global node of ny is in any element other
C     than the ones it is conected to it is 'hanging'
      DO nonode=1,NPNODE(0,nr)
        np=NPNODE(nonode,nr)

      IF((NVJP(1,np).EQ.1).AND.(NWP(np,1).NE.-1)) THEN

C ??? Don't know why I thought this was good idea. Restoring above.
C new CS 8/7/99 only check nodes in list
C      DO nonode=1,NPLIST(0)
C        np=NPLIST(nonode)

C putting this back in since a node may no longer hang after a refine
C        NWP(np,1)=0
        FOUND=.FALSE.
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '      NPNE(1,1,ne),nr,NVJE(1,1,1,ne),SE(1,1,ne),
     '      XA(1,1,ne),XE,XP,ERROR,*9999)
          NITB=NIT(NBJ(1,ne))
          XD(1)=XP(1,1,1,np)
          XD(2)=XP(1,1,2,np)
          XD(3)=XP(1,1,3,np)
          DO ni=1,NITB
            XI(ni)=0.5d0
          ENDDO
          CALL DEXI_POINT(IBT,IDO,INP,LD,NBJ,ne,NITB,nr,
     '      0.d0,XE,XI,XI,XD,.FALSE.,ERROR,*9999)

          nb=NBJ(1,ne)
          IF((LD.NE.0).AND.(NIT(nb).EQ.2).AND.(NJT.EQ.3)) THEN
            ! 2D element in 3D space, check 3rd component is correct
            XD(3)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
     '        XE(1,3))
            IF ((XD(3).GT.(XP(1,1,3,np)+SE(1,nb,ne)/4.0d0))
     '        .OR.(XD(3).LT.(XP(1,1,3,np)-SE(1,nb,ne)/4.0d0))) THEN
              LD=0
            ENDIF
          ENDIF

          DO no_nelist=1,NENP(np,0,nr)
            NELIST(no_nelist)=NENP(np,no_nelist,nr)
          ENDDO
          IF(((LD.NE.0).AND..NOT.
     '      (INLIST(ne,NELIST(1),NENP(np,0,nr),
     '      pos_inlist))).AND.(.NOT.FOUND)) THEN

            FOUND=.TRUE.
            IF(ne.NE.NWP(np,1)) THEN

            NWP(np,1)=ne !element np hangs in

C           Find local line/face np hangs on
C           SECTOR=.FALSE.
C           DO ni=1,NIT(nb)
C             IF(IBT(1,ni,nb).EQ.5.OR.IBT(1,ni,nb).EQ.6) SECTOR=.TRUE.
C           ENDDO !ni

            IF((XI(1).LT.0.05d0).OR.(XI(1).GT.0.95d0)) THEN
              ! on face 1 or 2
              IF(XI(1).LT.0.1d0)THEN
                NWP(np,2)=1
              ELSE
                NWP(np,2)=2
              ENDIF
            ENDIF

            IF((XI(2).LT.0.05d0).OR.(XI(2).GT.0.95d0)) THEN
              ! on face 3 or 4
              IF(XI(2).LT.0.1d0)THEN
                NWP(np,2)=3
              ELSE
                NWP(np,2)=4
              ENDIF
            ENDIF

            IF((XI(3).LT.0.05d0).OR.(XI(3).GT.0.95d0)) THEN
              ! on face 5 or 6
              IF(XI(3).LT.0.1d0)THEN
                NWP(np,2)=5
              ELSE
                NWP(np,2)=6
              ENDIF
            ENDIF

C            DISTANCE=0.5d0
C            DO ni=1,NIT(NBJ(1,ne)) !ni is direction normal to face/line
C              IF(XI(ni).GT.0.5d0) THEN
C                IF((1.0d0-XI(ni)).LT.DISTANCE) THEN
C                  DISTANCE=1.0d0-XI(ni)
C                  DIRECTION=ni
C                  CONST=1
C                ENDIF
C              ELSE
C                IF(XI(ni).LT.DISTANCE) THEN
C                  DISTANCE=XI(ni)
C                  DIRECTION=ni
C                  CONST=0
C                ENDIF
C              ENDIF
C            ENDDO ! ni
C**** Should do this a little better
C**** not right for sectors anyway
C            IF(DIRECTION.EQ.1) THEN
C              IF(CONST.EQ.0) NWP(np,2)=3
C              IF(CONST.EQ.1) NWP(np,2)=4
C            ELSE IF(DIRECTION.EQ.2) THEN
C              IF(SECTOR) THEN
C                NWP(np,2)=1
C              ELSE
C                IF(CONST.EQ.0) NWP(np,2)=1
C                IF(CONST.EQ.1) NWP(np,2)=2
C              ENDIF
C            ELSE
C              IF(CONST.EQ.0) NWP(np,2)=5
C              IF(CONST.EQ.1) NWP(np,2)=6
C            ENDIF

C This is taken from REFINE_SETNODE
C Sets scale factors by calculation
            NITE=NIT(nb)
            CALCDSDXI=.FALSE.

            nj=1    !just set scale factors for nj=1 for now
            nb=NBJ(nj,ne)
            IF(nb.GT.0) THEN
              IF(NBI(nb).GE.5.AND.NBI(nb).LE.7) CALCDSDXI=.TRUE.
            ENDIF !nb>0
            IF(CALCDSDXI) THEN
              IF(ITYP10(nr).GE.2) THEN
                nb=NBJ(1,ne)
                X(1,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '              1,XI,XE(1,1))
C               If an angle is 2*pi set it to zero
                IF(ITYP10(nr).EQ.3) THEN
                  nb=NBJ(3,ne)
                  X(1,3)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '              1,XI,XE(1,3))
                  IF(DABS(X(1,3)-2.0d0*PI).LT.LOOSE_TOL) X(1,3)=0.0d0
                ENDIF
                IF(ITYP10(nr).EQ.4) THEN
                  nb=NBJ(2,ne)
                  X(1,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '              1,XI,XE(1,2))
                  IF(DABS(X(1,2)-2.0d0*PI).LT.LOOSE_TOL) X(1,2)=0.0d0
                ENDIF
              ENDIF
              DO ni=1,NITE
                nu=1+ni*(1+ni)/2
                DO nj1=1,NJT
                  nb=NBJ(nj1,ne)
                  X(nu,nj1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '              INP(1,1,nb),nb,nu,XI,XE(1,nj1))
                ENDDO !nj1
                SUM=0.0d0
                IF(ITYP10(nr).EQ.1) THEN
                  SUM=X(nu,1)**2+X(nu,2)**2
                  IF(NJT.GT.2) SUM=SUM+X(nu,3)**2
                ELSE IF(ITYP10(nr).EQ.2) THEN
                  R=X(1,1)
                  RR=R*R
                  SUM=SUM+X(nu,1)**2+RR*X(nu,2)**2
                  IF(NJT.GT.2) SUM=SUM+X(nu,3)**2
                ELSE IF(ITYP10(nr).EQ.3) THEN
                  R=X(1,1)
                  RR=R*R
                  RC=R*DCOS(X(1,3))
                  RRC=RC*RC
                  SUM=SUM+X(nu,1)**2+RRC*X(nu,2)**2+RR*X(nu,3)**2
                ELSE IF(ITYP10(nr).EQ.4) THEN
                  AA=FOCUS*FOCUS
                  SLX=DSINH(X(1,1))
                  SMX=DSIN(X(1,2))
                  G1=AA*(SLX*SLX+SMX*SMX)
                  G3=AA* SLX*SLX*SMX*SMX
                  SUM=SUM+G1*(X(nu,1)**2+X(nu,2)**2)
                  IF(NJT.GT.2) SUM=SUM+G3*X(nu,3)**2
                ENDIF
                DSDXI(ni)=DSQRT(SUM)
                IF(DABS(DSDXI(ni)).LT.LOOSE_TOL) THEN
                  DSDXI(ni)=1.0d0
                  WRITE(OP_STRING,
     '              '('' >>Warning: Arc length derivative is zero'')')
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDIF
              ENDDO !ni
              IF(DOP) THEN
CC$              call mp_setlock()
                WRITE(OP_STRING,'('' DSDXI(ni):'',3(X,D12.5))')
     '            (DSDXI(ni),ni=1,NITE)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
              ENDIF
            ENDIF !CALCDSDXI

C           Number of derivatives needed at the new node
            NKJNEED=NKT(1,nb)  ! this will do for now
            DO nk=1,NKJNEED
              nu=NUNK(nk)
              X(nu,nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '          nb,nu,XI,XE(1,nj))
              IF((nj.EQ.2.AND.ITYP10(nr).GE.2).OR.
     '          (nj.EQ.3.AND.ITYP10(nr).GE.3)) THEN
C KAT 18Oct00:  XPXE often modifies nodal angles usually by adding 2pi.
C               Here we make a half hearted attempt to keep them
C               within [0,2pi).
                XTEMP=X(1,nj)-2d0*PI
                IF(XTEMP.GE.0d0) X(1,nj)=XTEMP
CC***  If an angle > 2*pi set it to zero
C                IF(X(1,nj)-2.0d0*PI).LT.LOOSE_TOL) X(1,nj)=0.0d0
              ENDIF
C              ENDIF
            ENDDO !nk

            IF(NBI(nb).GE.5.AND.NBI(nb).LE.7) THEN
C             arc-length or ave. arc-length scale factors
              DO ni=1,NITE
                SCALEFACTOR(ni)=DSDXI(ni)
              ENDDO !ni
            ELSE
              DO ni=1,3
                SCALEFACTOR(ni)=1.0d0
              ENDDO !ni
            ENDIF !nbi
            NOCROSS=.FALSE.     !this will do for now

C For the refined direction set the scale factor from interpolation
            FOUND_LINEAR_BASIS=.FALSE.
            nb2=0
            DO WHILE((.NOT.FOUND_LINEAR_BASIS)
     '        .AND.(nb2.LE.NBFT))
              nb2=nb2+1
              IF(NIT(nb).EQ.3) THEN
                IF((IBT(1,1,nb2).EQ.1).AND.
     '            (IBT(1,2,nb2).EQ.1).AND.(IBT(1,3,nb2).EQ.1)) THEN
                  FOUND_LINEAR_BASIS=.TRUE.
                ENDIF
              ELSE
                IF((IBT(1,1,nb2).EQ.1).AND.
     '            (IBT(1,2,nb2).EQ.1)) THEN
                  FOUND_LINEAR_BASIS=.TRUE.
                ENDIF
              ENDIF
            ENDDO
            CALL ASSERT(nb2.LE.NBFT,
     '        '>>Define linear basis',ERROR,*9999)

            DO nb3=1,NBFT
              DO nk=1,NKT(1,nb3)
                SP(nk,nb3,np)=1.0D0
              ENDDO !ns
C !!! not setting up faces or lower bases
              IF(NIT(nb3).EQ.3) THEN

              IF(NKT(0,nb3).GT.1.AND.NBI(nb3).NE.1) THEN !scale fac not unit
                DO nk2=2,NKT(1,nb3)
                  ns=0
                  DO nn=1,NNT(nb3)
                    ns=ns+1
                    DO nk=2,NKT(nn,nb3)
                      ns=ns+1
                      IF(nk.EQ.nk2) XS(nn)=SE(ns,nb3,ne)
                    ENDDO !nk
                  ENDDO !nn
                  SP(nk2,nb3,np)=
     '              PXI(IBT(1,1,nb2),IDO(1,1,0,nb2),
     '              INP(1,1,nb2),nb2,1,XI,XS)
                  !only set in refined direction
C!!! setting using basis 1 for now
                  IF(IDRN.EQ.1) THEN
                    IF((nk2.EQ.2).AND.(nb3.EQ.1)) THEN
                       SCALEFACTOR(1)=SP(nk2,nb3,np)
                    ENDIF
                  ELSE IF(IDRN.EQ.2) THEN
                    IF((nk2.EQ.3).AND.(nb3.EQ.1)) THEN
                       SCALEFACTOR(2)=SP(nk2,nb3,np)
                    ENDIF
                  ELSE
                    IF((nk2.EQ.5).AND.(nb3.EQ.1)) THEN
                       SCALEFACTOR(3)=SP(nk2,nb3,np)
                    ENDIF
                  ENDIF
                ENDDO !nk2
              ENDIF
              ENDIF
            ENDDO !nb3

C Set up the scale factors to use at the hanging node
            SP(1,nb,np)=1.0d0
            DO nk=2,NKJNEED
              nu=NUNK(nk)
              IF(nu.EQ.2) THEN !dXi1
                SP(nk,nb,np)=SCALEFACTOR(1)
              ELSE IF(nu.EQ.4) THEN !dXi2
                SP(nk,nb,np)=SCALEFACTOR(2)
              ELSE IF(nu.EQ.7) THEN !dXi3
                SP(nk,nb,np)=SCALEFACTOR(3)
              ELSE IF(NOCROSS) THEN
                SP(nk,nb,np)=0.0d0
              ELSE IF(nu.EQ.6) THEN !dXi1dXi2
                SP(nk,nb,np)=SCALEFACTOR(1)*SCALEFACTOR(2)
              ELSE IF(nu.EQ.9) THEN !dXi1dXi3
                SP(nk,nb,np)=SCALEFACTOR(1)*SCALEFACTOR(3)
              ELSE IF(nu.EQ.10) THEN !dXi2dXi3
                SP(nk,nb,np)=SCALEFACTOR(2)*SCALEFACTOR(3)
              ELSE IF(nu.EQ.11) THEN !dXi1dXi2dXi3
                SP(nk,nb,np)=
     '            SCALEFACTOR(1)*SCALEFACTOR(2)*SCALEFACTOR(3)
              ELSE
                WRITE(OP_STRING,'('' >>Warning: Derivative '
     '            //'not updated. nu='',I2,'' is unknown'')') nu
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO !nk

          ELSE
C          write(*,*) "node ",np," was in",NWP(np,1)," now in ",ne

          ENDIF

          ENDIF
        ENDDO ! noelem
        IF(.NOT.FOUND) THEN
          NWP(np,1)=0
        ENDIF

      ENDIF

      ENDDO ! nonode

      CALL EXITS('HANGING_NODE_DETECT')
      RETURN
 9999 CALL ERRORS('HANGING_NODE_DETECT',ERROR)
      CALL EXITS('HANGING_NODE_DETECT')
      RETURN 1
      END


