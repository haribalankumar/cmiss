      SUBROUTINE OPELEMD(NBH,NBJ,ne,NFF,
     '  NHE,NKHE,NLL,NPNE,nr,NRE,NVHE,nx,VOLTC,
     '  PG,RG,SE,VOL,VOLT,WG,ZE,ZG,OUTPUT_TOT_ONLY,ERROR,*)

C#### Subroutine: OPELEMD
C###  Description:
C###    OPELEMD outputs deformed element topology for element ne.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBH(NHM),NBJ(NJM),ne,NFF(6),NHE,NKHE(NKM,NNM,NHM),
     '  NLL(12),NPNE(NNM,NBFM),nr,NRE(NEM),NVHE(NNM,NBFM,NHM),nx,
     '  VOLTC(NBFM)
      REAL*8 PG(NSM,NUM,NGM,NBM),RG(NGM),SE(NSM,NBFM),
     '  VOL(NBFM),VOLT(NBFM),WG(NGM,NBM),
     '  ZE(NSM,NHM),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
      LOGICAL OUTPUT_TOT_ONLY
!     Local Variables
      INTEGER i,j,nb,nf,ng,nh,nhx,nk,nl,nn,NNK,np,ns
      REAL*8 DXIX(3,3),GZ,GZL(3,3),GZU(3,3),RAD,RMULT,RWG
      CHARACTER CHAR1*10,NAME1(3)*6,NV_STRING*10,OPT3(5)*21
      LOGICAL NBLOG,SHOW_VERSION

      DATA NAME1/'Length','  Area','Volume'/
      DATA OPT3(1)/'rectangular cartesian'/,
     '     OPT3(2)/'cylindrical polar    '/,
     '     OPT3(3)/'spherical polar      '/,
     '     OPT3(4)/'prolate spheroidal   '/,
     '     OPT3(5)/'oblate  spheroidal   '/

      CALL ENTERS('OPELEMD',*9999)
      nb=NBJ(1)

      DO i=1,3
        DO j=1,3
          DXIX(i,j)=0.0d0
        ENDDO
      ENDDO

      VOL(nb)=0.0D0
      DO ng=1,NGT(nb)
        IF(DOP) THEN
          WRITE(OP_STRING,'(/'' Gauss pt '',I3)') ng
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'( '' ------------''/)')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF !DOP

        CALL ZEZG(0,NBH(1),ng,NHE,nx,DXIX,PG,ZE,ZG,ERROR,*9999)
        CALL ZGMG(NBH(1),nr,GZ,GZL,GZU,ZG,ERROR,*9999)
c          IF(NIT(nb).EQ.1) RG(ng)=DSQRT(DABS(GZL(1,1)))
c          IF(NIT(nb).EQ.2) RG(ng)=DSQRT(DABS(DET(GZL)*GZU(1,1)))
c          IF(NIT(nb).EQ.3) RG(ng)=DSQRT(DABS(DET(GZL)))
        RG(ng)=DSQRT(GZ)
        RWG=RG(ng)*WG(ng,nb)

C LKC 7-FEB-1999 shallow water ?
C If for shallow water ?
C        IF((ni-NIT(nb)+1).EQ.1) THEN
          RMULT=1.0d0
          IF(JTYP4.EQ.2) THEN
            IF(ITYP11(nr).EQ.1) THEN
              RAD=ZG(2,1)
            ELSE IF(ITYP11(nr).EQ.2.OR.ITYP11(nr).EQ.3) THEN
              RAD=ZG(1,1)
            ELSE IF(ITYP11(nr).EQ.4) THEN
              RAD=FOCUS*DSINH(ZG(1,1))*DSIN(ZG(2,1))
            ELSE IF(ITYP11(nr).EQ.5) THEN
            ENDIF
            RMULT=2.0d0*PI*RAD
          ELSE IF(JTYP4.EQ.3) THEN
            IF(ITYP11(nr).EQ.1) THEN
              RAD=ZG(1,1)
            ELSE IF(ITYP11(nr).EQ.2.OR.ITYP11(nr).EQ.3) THEN
              RAD=ZG(1,1)
            ELSE IF(ITYP11(nr).EQ.4) THEN
              RAD=FOCUS*DSINH(ZG(1,1))*DSIN(ZG(2,1))
            ELSE IF(ITYP11(nr).EQ.5) THEN
            ENDIF
            RMULT=2.0d0*PI*RAD
          ELSE IF(JTYP4.EQ.4) THEN
            RMULT=4.0d0*PI*ZG(1,1)**2
          ENDIF
          VOL(nb)=VOL(nb)+RWG*RMULT

          IF(DOP) THEN
            WRITE(OP_STRING,'(''  GZ:'',D12.4,''  RWG:'',D12.4,'
     '        //'''  RMULT:'',D12.4, )') GZ, RWG, RMULT
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF

C LKC 7-FEB-1999 shallow water ?
C        ELSE
C          VOL(nb)=VOL(nb)+RWG*ZG(NHE,1)
C        ENDIF
      ENDDO !ng
      VOLT(nb)=VOLT(nb)+VOL(nb)
      VOLTC(nb)=VOLTC(nb)+1


      IF(.NOT.OUTPUT_TOT_ONLY) THEN
        WRITE(OP_STRING,'(/''  Element '',I5,''  (Region '',I1,
     '    //'')  Coordinates are '','
     '    //'A,3(2X,A6,''='',D12.5))')
     '    ne,NRE(ne),OPT3(ITYP10(NRE(ne))),
     '    NAME1(NIT(nb)+JTYP4-1),VOL(nb)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(2X,12(''=''))')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        DO nhx=1,NJ_LOC(NJL_GEOM,0,nr)
          nh=NH_LOC(nhx,nx)
          nb=NBJ(nh)
          IF(nb.GT.0) THEN
            IF(nh.LE.NJ_LOC(NJL_GEOM,0,nr).OR.
     '        (NIT(nb).EQ.NIT(NBJ(1)))) THEN
C ***         #field Xi's = #geom Xi's
              WRITE(OP_STRING,'(4X,''Geometric variable no. '',I2,'
     '          //''' (Basis function type= '',I2,''):'')') nh,nb
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
!news MPN 17-Nov-94
              SHOW_VERSION=.FALSE.
              DO nn=1,NNT(nb)
                IF(NVHE(nn,nb,nh).GT.1) SHOW_VERSION=.TRUE.
              ENDDO !nn
              ns=0
              DO nn=1,NNT(nb)
                np=NPNE(nn,nb)
                IF(SHOW_VERSION) THEN !several versions of nodal value
                  WRITE(CHAR1,'(I2)') NVHE(nn,nb,nh)
                  NV_STRING='['//CHAR1(1:2)//']'
                  IF(nn.EQ.1) THEN
                    FORMAT='(4X,''Node '',I2,''  Global no.'',I5,'''
     '                //NV_STRING(1:4)//'  Coord. derivs: '',4D12.4,'
     '                //'(/49X,7D12.4))'
                  ELSE
                    FORMAT='(9X,I2,12X,I5,'''//NV_STRING(1:4)
     '                //''',17X,4D12.4,(/49X,4D12.4))'
                  ENDIF
                ELSE
                  IF(nn.EQ.1) THEN
                    FORMAT='(4X,''Node '',I2,''  Global no.'',I5,'
     '                //'''  Coord. derivs: '',4D12.4,(/45X,4D12.4))'
                  ELSE
                    FORMAT='(9X,I2,12X,I5,17X,4D12.4,(/45X,4D12.4))'
                  ENDIF
                ENDIF
C AJP 13-6-93   NNK=(nn-1)*NKT(0,nb)
                IF(nn.GT.1) ns=ns+NKT(nn-1,nb) !AJP 13-6-93
                WRITE(OP_STRING,FORMAT)
     '            nn,np,(ZE(nk+ns,nhx),nk=1,NKT(nn,nb))
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDDO
              DO nn=1,NNT(nb)
                WRITE(OP_STRING,
     '            '(4X,''Node '',I2,''  Global derivs: '','//'7I2)') nn,
     '            (NKHE(nk,nn,nh),nk=1,NKT(0,nb))
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDDO
            ELSE IF(nh.GT.NJ_LOC(NJL_GEOM,0,nr).AND.
     '          (NIT(nb).ne.NIT(NBJ(1)))) THEN
              WRITE(OP_STRING,'('' Fibre or field not defined'')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF
        ENDDO !nhx

        DO nb=1,NBFT
          NBLOG=.FALSE.
          DO nhx=1,NH_LOC(0,nx)
            nh=NH_LOC(nhx,nx)
            IF(nb.EQ.NBH(nh)) NBLOG=.TRUE.
          ENDDO
          IF(NNT(nb).GT.1.AND.NKT(0,nb).GT.1.
     '      AND.(NBLOG.OR.NIT(nb).EQ.NIT(NBJ(1)))) THEN
            WRITE(OP_STRING,'(/4X,''Basis function '',I2)') nb
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            DO nn=1,NNT(nb)
              IF(nn.EQ.1) THEN
                FORMAT='(4X,''Node '',I2,''  Global no.'',I5,'
     '            //'''  Scale factors: '',4D12.4,(/45X,4D12.4))'
              ELSE
                FORMAT='(9X,I2,12X,I5,17X,4D12.4,(/45X,4D12.4))'
              ENDIF
              NNK=(nn-1)*NKT(0,nb)
              WRITE(OP_STRING,FORMAT) nn,NPNE(nn,nb),
     '          (SE(nk+NNK,nb),nk=1,NKT(0,nb))
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDDO
          ENDIF
        ENDDO

        nb=NBJ(1)
        IF(NNT(nb).GT.0) THEN
          WRITE(OP_STRING,'(/4X,''Global line nos. of element'
     '      //' edges : '',/10X,12I5)') (NLL(nl),nl=1,NLE(nb))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IF(NIT(nb).GT.2) THEN
            WRITE(OP_STRING,'(/4X,''Global face nos. of element '','
     '        //'''sides : '',/10X,12I5)') (NFF(nf),nf=1,NFE(nb))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
      ENDIF

      CALL EXITS('OPELEMD')
      RETURN
 9999 CALL ERRORS('OPELEMD',ERROR)
      CALL EXITS('OPELEMD')
      RETURN 1
      END


