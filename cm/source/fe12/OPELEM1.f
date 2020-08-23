      SUBROUTINE OPELEM1(ELEM_VOL,IBT,NBJ,ne,NFF,nj_volume,NKJE,NLL,
     '  NPNE,nr,NRE,NVJE,VOLTC,PG,RG,SE,VOL,VOLT,WG,XAB,XE,XG,
     '  OUTPUT_TOT_ONLY,ERROR,*)

C#### Subroutine: OPELEM1
C###  Description:
C###    OPELEM1 outputs undeformed element topology for element ne.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      REAL*8 ONESIXTH
      PARAMETER(ONESIXTH = 0.16666666666666666667D0)
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),
     '  NBJ(NJM),ne,NFF(6),nj_volume,NKJE(NKM,NNM,NJM),NLL(12),
     '  NPNE(NNM,NBFM),nr,NRE(NEM),NVJE(NNM,NBFM,NJM),VOLTC(NBFM)
      REAL*8 PG(NSM,NUM,NGM,NBM),RG(NGM),SE(NSM,NBFM),
     '  VOL(NBFM),VOLT(NBFM),WG(NGM,NBM),XAB(NORM,NEM),XE(NSM,NJM),
     '  XG(NJM,NUM)
      CHARACTER ERROR*(*)
      LOGICAL ELEM_VOL,OUTPUT_TOT_ONLY
!     Local Variables
      INTEGER IEND,mk,nb,nf,ng,NITB,nj,njj,njl,nl,
     '  nn,np,ns
      REAL*8 A(3,3),DET,DXRCXI(3,3),GL(3,3),GU(3,3),RAD,RMULT,RWG
      CHARACTER NAME1(3)*6
      LOGICAL CONSISTENT,POSITIVE

      DATA NAME1/'Length','  Area','Volume'/

      CALL ENTERS('OPELEM1',*9999)

      nb=NBJ(1)
      NITB=NIT(nb)

      CONSISTENT=.TRUE.
      IF(IBT(1,1,nb).EQ.3.AND.NITB.EQ.3) THEN !3D simplex
        DO nn = 2,4
          DO nj = 1,3
            A(nn-1,nj) = XE(nn,nj) - XE(1,nj)
          ENDDO
        ENDDO
        VOL(nb)  = DABS(ONESIXTH*DET(A))
        VOLT(nb) = VOLT(nb) + VOL(nb)
        VOLTC(nb) = VOLTC(nb)+1
      ELSE
        VOL(nb)=0.0d0
        DO ng=1,NGT(nb)
          CALL XEXG(NBJ(1),ng,nr,PG,XE,XG,ERROR,*9999)
          IF(NJ_LOC(NJL_GEOM,0,nr).NE.NITB) THEN
            CALL XGMG(-1,NITB,nb,nr,A,GL,GU,RG(ng),XG,ERROR,*9999)
          ELSE
C KAT 3Dec98: Checking consistent Jacobian of integration
            CALL DXRCDXI(NITB,nr,DXRCXI,XG,ERROR,*9999)
            IF(NITB.EQ.1) THEN
              RG(ng)=DXRCXI(1,1)
            ELSE IF(NITB.EQ.2) THEN
              RG(ng)=DXRCXI(1,1)*DXRCXI(2,2)-DXRCXI(1,2)*DXRCXI(2,1)
            ELSE !NITB=3
              RG(ng)=DET(DXRCXI)
            ENDIF
            IF(ng.EQ.1) THEN
              POSITIVE=RG(ng).GT.0.0D0
            ELSE IF(RG(ng).GT.0.0D0.NEQV.POSITIVE) THEN
              CONSISTENT=.FALSE.
            ENDIF
            IF(.NOT.POSITIVE) RG(ng)=-RG(ng)
          ENDIF
          RWG=RG(ng)*WG(ng,nb)
C LKC 3-FEB-1999 Shallow water ?
C          IF(ni.EQ.NITB) THEN
            RMULT=1.0d0
            IF(JTYP4.EQ.2) THEN
              IF(ITYP10(nr).EQ.1) THEN
                RAD=XG(2,1)
              ELSE IF(ITYP10(nr).EQ.2.OR.ITYP10(nr).EQ.3) THEN
                RAD=XG(1,1)
              ELSE IF(ITYP10(nr).EQ.4) THEN
                RAD=FOCUS*DSINH(XG(1,1))*DSIN(XG(2,1))
              ELSE IF(ITYP10(nr).EQ.5) THEN
              ENDIF
              RMULT=2.0d0*PI*RAD
            ELSE IF(JTYP4.EQ.3) THEN
              IF(ITYP10(nr).EQ.1) THEN
                RAD=XG(1,1)
              ELSE IF(ITYP10(nr).EQ.2.OR.ITYP10(nr).EQ.3) THEN
                RAD=XG(1,1)
              ELSE IF(ITYP10(nr).EQ.4) THEN
                RAD=FOCUS*DSINH(XG(1,1))*DSIN(XG(2,1))
              ELSE IF(ITYP10(nr).EQ.5) THEN
              ENDIF
              RMULT=2.0d0*PI*RAD
            ELSE IF(JTYP4.EQ.4) THEN
              RMULT=4.0d0*PI*XG(1,1)**2
            ENDIF
            VOL(nb)=VOL(nb)+RWG*RMULT
C LKC 3-FEB-1999 Shallow water ?
C          ELSE
C            VOL(nb)=VOL(nb)+RWG*XG(NJ_LOC(NJL_GEOM,0,nr),1)
C          ENDIF
        ENDDO !ng
        VOLT(nb)=VOLT(nb)+VOL(nb)
        VOLTC(nb)=VOLTC(nb)+1
      ENDIF !3d simplex

      IF(ELEM_VOL)THEN
        XAB(nj_volume,ne)=VOL(nb)
      ENDIF
        

      IF(.NOT.OUTPUT_TOT_ONLY) THEN
        WRITE(OP_STRING,'(1P/'' Element'',I9'
     '    //'/'' ================   (Region'',I3,'')  '','
     '    //'A6,''='',E12.5)') ne,NRE(ne),NAME1(NIT(nb)+JTYP4-1),VOL(nb)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(.NOT.CONSISTENT) THEN
          WRITE(OP_STRING,'('' >>WARNING: Integration Jacobian '
     '      //'passes through zero'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
        DO njl=1,3  !geom/fibres/field
          DO njj=1,NJ_LOC(njl,0,nr)
            nj=NJ_LOC(njl,njj,nr)
            nb=NBJ(nj)
            IEND=0
            CALL WRITE_LINE(IOFI,' ',ERROR,*9999)
            CALL APPENDC(IEND,'   ',OP_STRING(1))
            IF(njl.EQ.NJL_GEOM) THEN
              CALL APPENDC(IEND,'Geometric',OP_STRING(1))
            ELSEIF(njl.EQ.NJL_FIBR) THEN
              CALL APPENDC(IEND,'Fibre',OP_STRING(1))
            ELSE
              CALL APPENDC(IEND,'Field',OP_STRING(1))
            ENDIF
            CALL APPENDC(IEND,' variable ',OP_STRING(1))
            CALL APPENDI(IEND,njj,OP_STRING(1))
            CALL APPENDC(IEND,':  (Basis function type: ',OP_STRING(1))
            CALL APPENDI(IEND,nb,OP_STRING(1))
            CALL APPENDC(IEND,')',OP_STRING(1))
            CALL WRITE_LINE(IOFI,OP_STRING(1)(:IEND),ERROR,*9999)
            IF(nb.EQ.0) THEN
              CALL WRITE_LINE(IOFI,'     No interpolation defined.',
     '          ERROR,*9999)
            ELSE
              CALL WRITE_LINE(IOFI,
     '          '         Node             Derivative'//
     '          '         Value     Scale factor',ERROR,*9999)
              CALL WRITE_LINE(IOFI,
     '          '  local global version   local global',ERROR,*9999)
              ns=0
              DO nn=1,NNT(nb)
                np=NPNE(nn,nb)
                WRITE(OP_STRING,
     '            '(1P,1X,I4,I9,I6,I8,I7,6X,E12.5,'' /'',E12.5:'
     '            //'/(20X,I8,I7,6X,E12.5,'' /'',E12.5))')
     '            nn,np,NVJE(nn,nb,nj),
     '            (mk,NKJE(mk,nn,nj),XE(ns+mk,nj),SE(ns+mk,nb),
     '            mk=1,NKT(nn,nb))
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                ns=ns+NKT(nn,nb)
              ENDDO !nn
            ENDIF
          ENDDO  !njj
        ENDDO  !njl

C        DO nb=1,NBFT
CC GMH 6-10-96 Being a bit more picky about which basis to output
CC          NBLOG=.FALSE.
CC          DO nj=1,NJ_LOC(0,0,0)
CC            IF(nb.EQ.NBJ(nj)) NBLOG=.TRUE.
CC          ENDDO
CC          IF(NNT(nb).GT.1.AND.NKT(0,nb).GT.1.
CC     '      AND.(NBLOG.OR.NIT(nb).EQ.NIT(NBJ(1)))) THEN
C          SAMETYPE=.FALSE.
C          SAMENUMXI=.FALSE.
C          NBLOG=.FALSE.
C          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
C            IF(nb.EQ.NBJ(nj)) NBLOG=.TRUE.
CC PJH 3July96 Reduce basis o/p
CC           IF(NBC(nb).EQ.NBC(NBJ(nj)).OR.NBC(nb).EQ.7)
CC    '        SAMETYPE=.TRUE. !Same basis type or extended basis
CC           IF(NIT(nb).EQ.NIT(NBJ(nj))) SAMENUMXI=.TRUE.
C            IF(NBC(nb).EQ.7) SAMETYPE=.TRUE. !extended basis
C          ENDDO !nj
CC         IF(NNT(nb).GT.0.AND.SAMETYPE.AND.(NBLOG.OR.SAMENUMXI.OR.
C          IF(NNT(nb).GT.0.AND.(SAMETYPE.OR.NBLOG.OR.SAMENUMXI.OR.
C     '      IBT(1,NIT(nb),nb).EQ.9)) THEN
C            WRITE(OP_STRING,'(/4X,''Basis function '',I2)') nb
C            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C            DO nn=1,NNT(nb)
C              WRITE(OP_STRING,'(4X,''Node '',I2,''  Global derivs: '','
C     '          //'7I2)') nn,(NKE(nk,nn,nb),nk=1,NKT(nn,nb))
C              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C            ENDDO
C            DO nn=1,NNT(nb)
C              IF(nn.EQ.1) THEN
C                FORMAT='(4X,''Node '',I2,''  Global no.'',I5,'
C     '            //'''  Scale factors: '',4D12.4,(/45X,4D12.4))'
C              ELSE
C                FORMAT='(9X,I2,12X,I5,17X,4D12.4,(/45X,4D12.4))'
C              ENDIF
CC AJP         NNK=(nn-1)*NKT(0,nb)
C              ns=0  !Added AJP 4-6-93
C              DO NN2=1,nn-1
C                ns=ns+NKT(NN2,nb)
C              ENDDO
C              NNK=ns
C              WRITE(OP_STRING,FORMAT) nn,NPNE(nn,nb),
C     '          (SE(nk+NNK,nb),nk=1,NKT(nn,nb))
C              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C            ENDDO
C          ENDIF
C        ENDDO

        nb=NBJ(1)
        IF(NNT(nb).GT.0) THEN
          WRITE(OP_STRING,'(/3X,''Global line nos. of element'
     '      //' edges : '',/10X,12I5)') (NLL(nl),nl=1,NLE(nb))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IF(NIT(nb).GT.2) THEN
            WRITE(OP_STRING,'(/3X,''Global face nos. of element '','
     '        //'''sides : '',/10X,12I5)') (NFF(nf),nf=1,NFE(nb))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF

      ELSE IF(.NOT.CONSISTENT) THEN
        WRITE(OP_STRING,'('' >>WARNING: Integration Jacobian '
     '    //'passes through zero in element'',I9)') ne
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('OPELEM1')
      RETURN
 9999 CALL ERRORS('OPELEM1',ERROR)
      CALL EXITS('OPELEM1')
      RETURN 1
      END


