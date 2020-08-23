      SUBROUTINE SGISOC(INDEX,ISEG,ISISOC,iw,NBJ,
     '  NEELEM,NGAP,NKJE,NPF,NPNE,NRE,NVJE,
     '  SE,THRES,XA,XE,XP,CSEG,ERROR,*)

C#### Subroutine: SGISOC
C###  Description:
C###    SGISOC creates isochrone segment ISISOC(iw,ng,ne).
C###    A fill area is drawn around the Gauss point.
C###    The fill area index IFILL is ng+(ne-1)*ngt(nb).

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
!     Parameter List
      INTEGER INDEX,ISEG(*),ISISOC,iw,NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NGAP(NIM,NBM),
     '  NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NRE(NEM),NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 SE(NSM,NBFM,NEM),THRES(3,NGM,NEM),XA(NAM,NJM,NEM),
     '  XE(NSM,NJM),XP(NKM,NVM,NJM,NPM)
      CHARACTER CSEG(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER i,IFILL,INDEX_OLD,nb,ne,ng,ng1,ng2,NGG,noelem
      REAL*8 DXI1,DXI2,PTS(3,5),X(3),XI1,XI2,Z(3)
      CHARACTER CHAR*1

      CALL ENTERS('SGISOC',*9999)
      CALL OPEN_SEGMENT(ISISOC,ISEG,iw,'ISOC',INDEX,INDEX_OLD,
     '  1,1,CSEG,ERROR,*9999)

      DO i=1,3 !MHT 24-3-00 to avoid used before set, etc. errors
        X(i)=0.d0
      ENDDO !i
      DO noelem=1,NEELEM(0,1)
        ne=NEELEM(noelem,1)
        CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '    NRE(ne),NVJE(1,1,1,ne),
     '    SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
        nb=NBJ(1,ne)
        DXI1=1.0D0/DBLE(NGAP(1,nb))
        DXI2=1.0D0/DBLE(NGAP(2,nb))
        XI2=0.0D0
        ng=0
        DO ng2=1,NGAP(2,nb)
          XI1=0.0D0
          DO ng1=1,NGAP(1,nb)
            PTS(1,1)=((1.0D0-XI1)*XE(1,1)+XI1*XE(2,1))*(1.0D0-XI2)
     '             +((1.0D0-XI1)*XE(3,1)+XI1*XE(4,1))*XI2
            PTS(2,1)=((1.0D0-XI1)*XE(1,2)+XI1*XE(2,2))*(1.0D0-XI2)
     '             +((1.0D0-XI1)*XE(3,2)+XI1*XE(4,2))*XI2
            XI1=XI1+DXI1
            PTS(1,2)=((1.0D0-XI1)*XE(1,1)+XI1*XE(2,1))*(1.0D0-XI2)
     '             +((1.0D0-XI1)*XE(3,1)+XI1*XE(4,1))*XI2
            PTS(2,2)=((1.0D0-XI1)*XE(1,2)+XI1*XE(2,2))*(1.0D0-XI2)
     '             +((1.0D0-XI1)*XE(3,2)+XI1*XE(4,2))*XI2
            XI2=XI2+DXI2
            PTS(1,3)=((1.0D0-XI1)*XE(1,1)+XI1*XE(2,1))*(1.0D0-XI2)
     '             +((1.0D0-XI1)*XE(3,1)+XI1*XE(4,1))*XI2
            PTS(2,3)=((1.0D0-XI1)*XE(1,2)+XI1*XE(2,2))*(1.0D0-XI2)
     '             +((1.0D0-XI1)*XE(3,2)+XI1*XE(4,2))*XI2
            XI1=XI1-DXI1
            PTS(1,4)=((1.0D0-XI1)*XE(1,1)+XI1*XE(2,1))*(1.0D0-XI2)
     '             +((1.0D0-XI1)*XE(3,1)+XI1*XE(4,1))*XI2
            PTS(2,4)=((1.0D0-XI1)*XE(1,2)+XI1*XE(2,2))*(1.0D0-XI2)
     '             +((1.0D0-XI1)*XE(3,2)+XI1*XE(4,2))*XI2
            XI2=XI2-DXI2
            PTS(1,5)=PTS(1,4)
            PTS(2,5)=PTS(2,4)
            ng=ng+1
            IFILL=ng+(ne-1)*NGT(nb)
            IF(DOP) THEN
              WRITE(OP_STRING,'('' ne='',I3,'' ng='',I2,'
     '        //''' IFILL='',I3,'' u='',E12.3)')
     '        ne,ng,IFILL,THRES(2,ng,ne)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            IF(NJT.EQ.2) THEN
              CALL FILL_AREA(IFILL,iw,5,PTS,ERROR,*9999)
            ELSE IF(NJT.EQ.3) THEN
              IF(iw.NE.4) THEN
                CALL XZ(ITYP10(1),X,Z)
                CALL POLYMARKER(INDEX,iw,1,Z,ERROR,*9999)
              ELSE IF(iw.EQ.4) THEN
                NGG=INT((ng-1)/9)+1
                WRITE(CHAR,'(I1)') NGG
                CALL TEXT(1,iw,CHAR(1:1),X,ERROR,*9999)
              ENDIF
            ENDIF
            XI1=XI1+DXI1
          ENDDO
          XI2=XI2+DXI2
        ENDDO
      ENDDO

      CALL CLOSE_SEGMENT(ISISOC,iw,ERROR,*9999)
      CALL EXITS('SGISOC')
      RETURN
 9999 CALL ERRORS('SGISOC',ERROR)
      CALL EXITS('SGISOC')
      RETURN 1
      END


