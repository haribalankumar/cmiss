      SUBROUTINE LDAT(IBT,IDO,INP,IW,NBH,NBJ,NDLT,nr,
     '  STATIC,XE,XIDL,ZDL,ZDLMIN,ZE,ERROR,*)

C#### Subroutine: LDAT
C###  Description:
C###    LDAT draws lines connecting the data point ZDL in the element ne
C###    with its corresponding point on the element.
C###    The points on the ensemble are identified by a cross
C###    corresponding to small increments in the element coords XIDL.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'data00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp50.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),IW,NBH(NHM),NBJ(NJM),NDLT,nr
      REAL*8 XE(NSM,NJM),XIDL(NIM,NDEM),ZDL(NHM,NDEM),ZDLMIN,ZE(NSM,NHM)
      CHARACTER ERROR*(*)
      LOGICAL STATIC
!     Local Variables
      INTEGER n,nb,NBO,nde,ni,nj
      REAL*8 DXI,PXI,X0(4),X1(4),X2(4),XI1(4),XI2(4),Z0(4),Z1(4),Z2(4),
     '  Z3(4),Z4(4),ZL(3,2),LENGTH

      CALL ENTERS('LDAT',*9999)
      DO nde=1,NDLT
c AAY 29Aug95 allow data pt projections outside element
c       IF(XIDL(1,nde).GE.0.d0.AND.XIDL(1,nde).LE.1.d0.AND.
c    '     XIDL(2,nde).GE.0.d0.AND.XIDL(2,nde).LE.1.d0) THEN
          ZL(1,1)=ZDL(1,nde)
          ZL(2,1)=ZDL(2,nde)
          ZL(3,1)=ZDL(3,nde)
          DO n=1,3
            X0(n)=0.d0
            X1(n)=0.d0
            X2(n)=0.d0
          ENDDO
          DO nj=1,NJT
            nb=NBJ(nj)
            X0(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '        XIDL(1,nde),XE(1,nj))
            Z2(nj)=0.d0
            IF(.NOT.STATIC) THEN !motion
              NBO=NBH(nj)
              Z2(nj)=PXI(IBT(1,1,NBO),IDO(1,1,0,NBO),INP(1,1,NBO),
     '          NBO,1,XIDL(1,nde),ZE(1,nj))
            ENDIF
          ENDDO
          IF(DOP) THEN
            WRITE(OP_STRING,'(A)') ' ******LDAT'
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' X0(nj)='',4(D12.6,4X))')
     '        (X0(nj),nj=1,NJT)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          CALL XZ(ITYP10(1),X0,Z0)
          LENGTH=0.d0
          DO nj=1,NJT
            ZL(nj,2)=Z0(nj)
            IF(KTYP58(nr).EQ.2) THEN
              ZL(nj,2)=Z0(nj)+Z2(nj)
            ENDIF
            LENGTH=LENGTH+(ZL(nj,2)-ZL(nj,1))**2
          ENDDO
          LENGTH=DSQRT(LENGTH)

C CPB 12/1/93 if the length of the projection is greater than the specified
C minimum length then draw the projection

          IF(LENGTH.GE.ZDLMIN) THEN

            CALL POLYLINE(3,IW,2,ZL,ERROR,*9999)

            DO ni=1,NIT(NBJ(1))
              XI1(ni)=XIDL(ni,nde)
              XI2(ni)=XIDL(ni,nde)
            ENDDO
            IF(.NOT.STATIC) THEN !last xi is time
              XI1(NXIDEF)=XIDL(ni,nde)
              XI2(NXIDEF)=XIDL(ni,nde)
            ENDIF
            DXI=0.05D0
            DO ni=1,NIT(NBJ(1))
              XI1(ni)=XI1(ni)+DXI
              XI2(ni)=XI2(ni)-DXI
              DO nj=1,NJT
                nb=NBJ(nj)
                X1(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '            XI1,XE(1,nj))
                X2(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '            XI2,XE(1,nj))
                Z3(nj)=0.d0
                Z4(nj)=0.d0
                IF(.NOT.STATIC) THEN !motion
                  NBO=NBH(nj)
                  Z3(nj)=PXI(IBT(1,1,NBO),IDO(1,1,0,NBO),INP(1,1,NBO),
     '              NBO,1,XI1,ZE(1,nj))
                  Z4(nj)=PXI(IBT(1,1,NBO),IDO(1,1,0,NBO),INP(1,1,NBO),
     '              NBO,1,XI2,ZE(1,nj))
                ENDIF
              ENDDO
              CALL XZ(ITYP10(1),X1,Z1)
              CALL XZ(ITYP10(1),X2,Z2)
              DO nj=1,NJT
                ZL(nj,1)=Z1(nj)
                ZL(nj,2)=Z2(nj)
                IF(KTYP58(nr).EQ.2) THEN
                  ZL(nj,1)=Z1(nj)+Z3(nj)
                  ZL(nj,2)=Z2(nj)+Z4(nj)
                ENDIF
              ENDDO
              CALL POLYLINE(3,IW,2,ZL,ERROR,*9999)
              XI1(ni)=XIDL(ni,nde)
              XI2(ni)=XIDL(ni,nde)
            ENDDO
          ENDIF
c       ENDIF
      ENDDO

      CALL EXITS('LDAT')
      RETURN
 9999 CALL ERRORS('LDAT',ERROR)
      CALL EXITS('LDAT')
      RETURN 1
      END


