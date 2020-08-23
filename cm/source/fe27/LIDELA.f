      SUBROUTINE LIDELA(NBJ,NEELEM,NPNE,NVJE,NXI,XP,ZA,STRING,ERROR,*)

C#### Subroutine: LIDELA
C###  Description:
C###    LIDELA lists Delaunay simplices

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NPNE(NNM,NBFM,NEM),nr,
     '  NVJE(NNM,NBFM,NJM,NEM),NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 ZA(NAM,NHM,NCM,NEM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,N3CO,noelem,ne,nn,nj,nvi,npi,
     '  np1,nv1,nb,njj
      REAL*8 A(3,3),VOL,TOTVOL,DET
      LOGICAL CBBREV

      CALL ENTERS('LIDELA',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM list delaunay
C###  Description:
C###    Lists delaunay triangulation
C###  Parameter:      <region #[1]>
C###    total lists the total Voronoi statistics
C###    region lists Delaunay region #

        OP_STRING(1)=BLANK(1:15)//'<region #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE
        IF(CBBREV(CO,'REGION',1,noco+1,NTCO,N3CO)) THEN
          nr=IFROMC(CO(N3CO+1))
        ELSE
          nr=1
        ENDIF
        WRITE(OP_STRING,'('' Listing Delaunay Simplices:'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        TOTVOL=0.d0
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          nb=NBJ(1,ne)
          np1=NPNE(1,nb,ne)
          nv1=NVJE(1,nb,1,ne)
          DO nn=2,NNT(nb)
            npi=NPNE(nn,nb,ne)
            nvi=NVJE(nn,nb,1,ne)
            DO njj=1,NJ_LOC(NJL_GEOM,0,nr)
              nj=NJ_LOC(NJL_GEOM,njj,nr)
              A(nn-1,njj)=XP(1,nvi,nj,npi)-XP(1,nv1,nj,np1)
            ENDDO
          ENDDO
          IF(NJ_LOC(NJL_GEOM,0,nr).EQ.3) THEN
            VOL=DABS(0.166666666666667d0*DET(A))
            TOTVOL=TOTVOL+VOL
          ELSEIF(NJ_LOC(NJL_GEOM,0,nr).EQ.2) THEN
            VOL=DABS(0.5d0*(A(1,1)*A(2,2)-A(1,2)*A(2,1)))
            TOTVOL=TOTVOL+VOL
          ENDIF
          WRITE(OP_STRING,'(/,''Delaunay simplex:  '',I12)') ne
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''Circumcentre:      '',3D12.4)')
     '      (ZA(1,NJ_LOC(NJL_GEOM,nj,nr),1,ne),nj=1,NJ_LOC(NJL_GEOM,0,
     '      nr))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''Volume:            '',D12.4)') VOL
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''Connections:'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''============'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(''  Node  Vers   Opp      Coordinates'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          DO nn=1,NNT(nb)
            WRITE(OP_STRING,'(I6,I6,I6,''    '',3D12.4)')
     '        NPNE(nn,nb,ne),NVJE(nn,nb,1,ne),
     '        NXI(0,nn,ne),
     '        (XP(1,NVJE(nn,nb,1,ne),NJ_LOC(NJL_GEOM,nj,nr),
     '        NPNE(nn,nb,ne)),nj=1,NJ_LOC(NJL_GEOM,0,nr))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDDO !noelem
        WRITE(OP_STRING,'(/,''Total Volume:    '',D12.4)') TOTVOL
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('LIDELA')
      RETURN
 9999 CALL ERRORS('LIDELA',ERROR)
      CALL EXITS('LIDELA')
      RETURN 1
      END


