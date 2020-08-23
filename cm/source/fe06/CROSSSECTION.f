      SUBROUTINE CROSSSECTION(INDEX,IBT,IDO,INP,IW,
     '  NBJ,NEELEM,NENP,nj,NKJE,NNB,NPF,NPNE,NRE,NVJE,NXI,
     '  SE,XA,XE,XP,ZVAL,ERROR,*)

C#### Subroutine: CROSSSECTION
C###  Description:
C###    CROSSECTION draws sagittal crosssection of heart when sheet
C###    angles corrected.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INDEX,
     '  INP(NNM,NIM,NBFM),IW,NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NENP(NPM,0:NEPM,0:NRM),nj,
     '  NKJE(NKM,NNM,NJM,NEM),NNB(4,4,4,NBFM),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),NRE(NEM),NVJE(NNM,NBFM,NJM,NEM),
     '  NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 SE(NSM,NBFM,NEM),XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XP(NKM,NVM,NJM,NPM),ZVAL
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,j,n,n1elem,nb,ne,NE_ADJACENT,noelem,NPOINTS
      REAL*8 DTHETA,DY,DZ,POINTS(3,21),PXI,THETA1,THETA2,X(3),XI(3),
     '  XI2,XII(2,21),YZ(3,21),Z(3)

      CALL ENTERS('CROSSSECTION',*9999)

      CALL NENXI(NBJ,NEELEM,NENP,NNB,NPNE,NXI,ERROR,*9999)
      DO noelem=1,NEELEM(0,1)
        ne=NEELEM(noelem,1)
        IF(NXI(-3,1,ne).EQ.0.AND.
     '    (NET(1).LE.30.OR.(NET(1).gt.30.and.ne.GT.30))) THEN
!                              !restrict to endo temporary!!!!!
          CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '      NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '      SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)

          XI(3)=0.d0 !inner edge of endocardial element
          NPOINTS=0
C         DO XI2=0.d0,1.d0,0.05D0
          DO i=0,20
            XI2=DBLE(i)/20.0d0
            XI(1)=0.d0
            XI(2)=XI2
            nb=NBJ(nj,ne)
            THETA1=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
     '        XE(1,nj))
            XI(1)=1.d0
            THETA2=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
     '        XE(1,nj))
            IF(THETA1.GT.ZVAL.AND.ZVAL.GT.THETA2) THEN
              NE_ADJACENT=NXI(3,1,ne)
              IF(DOP) THEN
                WRITE(OP_STRING,'('' Element '',I4,'
     '            //''' adjacent element is '',I4)')ne,NE_ADJACENT
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              NPOINTS=NPOINTS+1
              IF(DOP) THEN
                WRITE(OP_STRING,'('' NPOINTS= '',I2)') NPOINTS
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              DTHETA=THETA1-THETA2
              IF(DABS(DTHETA).GT.1.d-6) THEN
                XI(1)=(THETA1-ZVAL)/DTHETA
              ELSE
                XI(1)=0.d0
              ENDIF
              IF(DOP) THEN
                WRITE(OP_STRING,'('' inner edge of endo element'',i4,'
     '            //''' xi(1)='',F5.3,'' xi(2)='',F5.3)') ne,XI(1),XI(2)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              DO j=1,3
                nb=NBJ(j,ne)
                X(j)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
     '            XE(1,j))
              ENDDO
              IF(DOP) THEN
                WRITE(OP_STRING,'('' Interpolated coords:'',3E12.3)')
     '            (X(j),j=1,3)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              CALL XZ(ITYP10(1),X,Z)
              POINTS(1,NPOINTS)=Z(1)
              POINTS(2,NPOINTS)=DSQRT(Z(2)**2+Z(3)**2)
              XII(1,NPOINTS)=XI(1) !to store Xi(1) to draw outer edge
              XII(2,NPOINTS)=XI(2) !to store Xi(2) to draw outer edge
              YZ(2,NPOINTS) =Z(2)  !retains y coord from inner edge
              YZ(3,NPOINTS) =Z(3)  !retains z coord from inner edge
            ENDIF
          ENDDO
          CALL POLYLINE(INDEX,IW,NPOINTS,POINTS,ERROR,*9999)

          XI(3)=1.d0 !outer edge of endocardial element
          DO n=1,NPOINTS
            XI(1)=XII(1,n) !from inner edge where theta is set
            XI(2)=XII(2,n)
            IF(DOP) THEN
              WRITE(OP_STRING,'('' outer edge of endo element'',i4,'
     '          //''' xi(1)='',F5.3,'' xi(2)='',F5.3)') ne,XI(1),XI(2)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            DO j=1,3
              nb=NBJ(j,ne)
              X(j)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
     '          XE(1,j))
            ENDDO
            IF(DOP) THEN
              WRITE(OP_STRING,'('' Interpolated coords:'',3E12.3)')
     '          (X(j),j=1,3)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            CALL XZ(ITYP10(1),X,Z)
            DY=Z(2)-YZ(2,n) !is change in y,z across element
            DZ=Z(3)-YZ(3,n) !is change in y,z across element
            IF(DOP) THEN
              WRITE(OP_STRING,'('' dY,dZ:'',2E12.3)') DY,DZ
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            POINTS(1,n)=Z(1)
            POINTS(2,n)=POINTS(2,N)+DSQRT(DY**2+DZ**2)
            YZ(2,n)    =Z(2)  !retains y coord from outer edge
            YZ(3,n)    =Z(3)  !retains z coord from outer edge
          ENDDO
          CALL POLYLINE(INDEX,IW,NPOINTS,POINTS,ERROR,*9999)

          IF(NE_ADJACENT.GT.0) THEN !draw outer bdry of adjacent element
            ne=NE_ADJACENT
            CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '        NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '        SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
            XI(3)=1.d0 !outer edge of epicardial element
            DO n=1,NPOINTS
              XI(1)=XII(1,n) !from inner edge of endo element
              XI(2)=XII(2,n) !..where theta is set
              IF(DOP) THEN
                WRITE(OP_STRING,'('' outer edge of epi element'',i4,'
     '            //''' xi(1)='',F5.3,'' xi(2)='',F5.3)') ne,XI(1),XI(2)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              DO j=1,3
                nb=NBJ(j,ne)
                X(j)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
     '            XE(1,j))
              ENDDO
              IF(DOP) THEN
                WRITE(OP_STRING,'('' Interpolated coords:'',3E12.3)')
     '            (X(j),j=1,3)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              CALL XZ(ITYP10(1),X,Z)
              DY=Z(2)-YZ(2,n) !is change in y,z across element
              DZ=Z(3)-YZ(3,n) !is change in y,z across element
              IF(DOP) THEN
                WRITE(OP_STRING,'('' dY,dZ:'',2E12.3)') DY,DZ
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              POINTS(1,n)=Z(1)
              POINTS(2,n)=POINTS(2,n)+DSQRT(DY**2+DZ**2)
            ENDDO
            CALL POLYLINE(INDEX,IW,NPOINTS,POINTS,ERROR,*9999)

          ELSE IF(NE_ADJACENT.EQ.0) THEN !ne is septum and need rv free wall
            THETA1=XE(1,3) !is theta at 1st node of endo element
            THETA2=XE(2,3) !is theta at 2nd node of endo element
            DO n1elem=1,NEELEM(0,1)
              ne=NEELEM(n1elem,1)
              IF(NET(1).gt.30.and.ne.LE.30) THEN !restrict to epi TEMPORARY!!!!
                CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '            NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '            SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
                IF(DABS(XE(1,3)-THETA1).LT.1.D-2.AND
     '            .DABS(XE(2,3)-THETA2).LT.1.D-2) THEN
                  WRITE(OP_STRING,'('' Found rv free wall element '','
     '              //'I4)') ne
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  GO TO 20
                ENDIF
              ENDIF
            ENDDO
 20         CONTINUE

            XI(3)=0.d0 !inner edge of epicardial element
            DO n=1,NPOINTS
              XI(1)=XII(1,n) !from inner edge of endo element
              XI(2)=XII(2,n) !..where theta is set
              IF(DOP) THEN
                WRITE(OP_STRING,'('' inner edge of epi element'',i4,'
     '            //''' xi(1)='',F5.3,'' xi(2)='',F5.3)') ne,XI(1),XI(2)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              DO j=1,3
                nb=NBJ(j,ne)
                X(j)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
     '            XE(1,j))
              ENDDO
              IF(DOP) THEN
                WRITE(OP_STRING,'('' Interpolated coords:'',3E12.3)')
     '            (X(j),j=1,3)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              CALL XZ(ITYP10(1),X,Z)
              DY=Z(2)-YZ(2,n) !is change in y,z across element
              DZ=Z(3)-YZ(3,n) !is change in y,z across element
              IF(DOP) THEN
                WRITE(OP_STRING,'('' dY,dZ:'',2E12.3)') DY,DZ
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              POINTS(1,n)=Z(1)
              POINTS(2,n)=POINTS(2,N)+DSQRT(DY**2+DZ**2)
              YZ(2,n)    =Z(2)  !retains y coord from inner edge
              YZ(3,n)    =Z(3)  !retains z coord from inner edge
            ENDDO
            CALL POLYLINE(INDEX,IW,NPOINTS,POINTS,ERROR,*9999)

            XI(3)=1.d0 !outer edge of epicardial element
            DO N=1,NPOINTS
              XI(1)=XII(1,n) !from inner edge of endo element
              XI(2)=XII(2,n) !..where theta is set
              IF(DOP) THEN
                WRITE(OP_STRING,'('' outer edge of epi element'',i4,'
     '            //''' xi(1)='',F5.3,'' xi(2)='',F5.3)') ne,XI(1),XI(2)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              DO j=1,3
                nb=NBJ(j,ne)
                X(j)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
     '            XE(1,j))
              ENDDO
              IF(DOP) THEN
                WRITE(OP_STRING,'('' Interpolated coords:'',3E12.3)')
     '            (X(j),j=1,3)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              CALL XZ(ITYP10(1),X,Z)
              DY=Z(2)-YZ(2,n) !is change in y,z across element
              DZ=Z(3)-YZ(3,n) !is change in y,z across element
              IF(DOP) THEN
                WRITE(OP_STRING,'('' dY,dZ:'',2E12.3)') DY,DZ
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              POINTS(1,n)=Z(1)
              POINTS(2,n)=POINTS(2,n)+DSQRT(DY**2+DZ**2)
            ENDDO
            CALL POLYLINE(INDEX,IW,NPOINTS,POINTS,ERROR,*9999)

          ENDIF
        ENDIF
      ENDDO

      CALL EXITS('CROSSSECTION')
      RETURN
 9999 CALL ERRORS('CROSSSECTION',ERROR)
      CALL EXITS('CROSSSECTION')
      RETURN 1
      END


