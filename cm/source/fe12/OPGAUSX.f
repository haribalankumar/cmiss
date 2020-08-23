      SUBROUTINE OPGAUSX(nb,NBH,NBJ,NELIST,NGLIST,
     '  NHE,NKHE,NKJE,NPF,NPNE,nr,NUMYGCMPTS,NVHE,NVJE,NW,nx,
     '  CURVCORRECT,PG,RG,SE,XA,XE,XG,XIG,XP,YG,ZA,ZE,ZG,ZP,
     '  TYPE,ERROR,*)

C#### Subroutine: OPGAUSX
C###  Description:
C###    OPGAUSX outputs element Gauss point array XG(nj,nu)
C###    and YG(niyg,ng,ne).

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'press00.cmn'
      INCLUDE 'inout00.cmn'

!     Parameter List
      INTEGER nx,nb,NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NELIST(0:NEM),
     '  NGLIST(0:NGM),NHE(NEM,NXM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),nr,NUMYGCMPTS,NVHE(NNM,NBFM,NHM,NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3,NXM)
      REAL*8 CURVCORRECT(2,2,NNM,NEM),PG(NSM,NUM,NGM,NBM),RG(NGM),
     '  SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),XIG(NIM,NGM,NBM),
     '  XP(NKM,NVM,NJM,NPM),YG(NIYGM,NGM,NEM),
     '  ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),ZG(NHM,NUM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),TYPE*(*)
!     Local Variables
      INTEGER nb_type,NBJ_temp(20),ne,ng,nh,nhx,niyg,nj,njj1,njj2,
     '  nolist,no_nglist,nu,ni
      REAL*8 DXIX(3,3),GL(3,3),GU(3,3),PRESSURE_value,Z(3)
      CHARACTER CHAR1*2

      CALL ENTERS('OPGAUSX',*9999)

      DO nolist=1,NELIST(0)
        ne=NELIST(nolist)
        WRITE(OP_STRING,'(/'' Element '',I4)') ne
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        IF(TYPE(1:8).EQ.'PRESSURE') THEN
          IF(ITYP1(nr,nx).EQ.3) THEN !potential flow from CMISS soln
            CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '        nr,NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
            DO no_nglist=1,NGLIST(0)
              ng=NGLIST(no_nglist)
              WRITE(OP_STRING,'('' Gauss point '',I3)') ng
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              CALL XEXG(NBJ(1,ne),ng,nr,PG,XE,XG,ERROR,*9999)
C              CALL XGMG(1,0,NBJ(1,ne),nr,DXIX,GL,GU,RG,XG,
C     '          ERROR,*9999)
              CALL XGMG(1,0,NBJ(1,ne),nr,DXIX,GL,GU,RG(ng),XG,
     '          ERROR,*9999)
              CALL ZEZG(1,NBH(1,1,ne),ng,NHE(ne,nx),nx,DXIX,PG,ZE,
     '          ZG,ERROR,*9999)
              PRESSURE_value=0.5D0*(ZG(1,2)**2+ZG(1,4)**2)
              WRITE(OP_STRING,'('' Pressure= '',D12.4)')
     '          PRESSURE_value
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDDO
          ELSE !potential flow for VSaero solution
            WRITE(OP_STRING,'('' Pressure: '',6D12.4,:/(11X,6D12.4))')
     '        (PRESS(NGLIST(no_nglist),ne),no_nglist=1,NGLIST(0))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF

        ELSE IF(TYPE(1:8).EQ.'SOLUTION') THEN
          CALL ZPZE(NBH(1,1,ne),1,NHE(ne,nx),NKHE(1,1,1,ne),NPF(1,1),
     '      NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1,nx),nx,
     '      CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),ZE,
     '      ZP,ERROR,*9999)
          DO no_nglist=1,NGLIST(0)
            ng=NGLIST(no_nglist)
            WRITE(OP_STRING,'('' Gauss point '',I3)') ng
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C GR 7Jul2009: add xi location to output if basis set
            IF(nb.NE.0) THEN
              FORMAT='( ''          '', ''  Xi(i): '',3D12.4)'
              WRITE(OP_STRING,FMT=FORMAT) (XIG(ni,ng,nb),ni=1,NIT(nb))
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
            CALL ZEZG(1,NBH(1,1,ne),ng,NHE(ne,nx),nx,DXIX,PG,ZE,
     '        ZG,ERROR,*9999)
            DO nhx=1,NH_LOC(0,nx)
              nh=NH_LOC(nhx,nx)
              WRITE(OP_STRING,'('' ZG(nh='',I2,'',nu=1..): '',6D12.4,'
     '          //':/(19X,6D12.4))')
     '          nh,(ZG(nh,nu),nu=1,NUT(NBH(nh,1,ne)))
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDDO !nhx (nh)
          ENDDO !no_nglist (ng)

        ELSE
          DO njj1=1,3   !geom/fibres/field
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              CALL ASSERT(nj.LE.20,'>>ERROR: increase size of NBJ_temp '
     '          //'to NJM',ERROR,*9999)
              IF(nb.EQ.0) THEN
                NBJ_temp(nj)=NBJ(nj,ne)
              ELSE
                NBJ_temp(nj)=nb
              ENDIF
            ENDDO !njj2
          ENDDO !njj1
          CALL XPXE(NBJ_temp,NKJE(1,1,1,ne),NPF(1,1),
     '      NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '      SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
          DO no_nglist=1,NGLIST(0)
            ng=NGLIST(no_nglist)
            IF(TYPE(1:2).EQ.'XG'.OR.TYPE(1:5).EQ.'XG&YG') THEN
              WRITE(OP_STRING,'('' Gauss point '',I3)') ng
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
            CALL XEXG(NBJ_temp,ng,nr,PG,XE,XG,ERROR,*9999)
            DO njj1=1,3   !geom/fibres/field
              DO njj2=1,NJ_LOC(njj1,0,nr)
                nj=NJ_LOC(njj1,njj2,nr)
                nb_type=NBJ_temp(nj)
                IF(TYPE(1:2).EQ.'XG'.OR.TYPE(1:5).EQ.'XG&YG') THEN
                  IF(nb_type.GT.0) THEN
                    WRITE(OP_STRING,'(1P,'' XG('',I2,'',nu):'',6E12.4,'
     '                //':/(11X,6E12.4))')
     '                nj,(XG(nj,nu),nu=1,NUT(nb_type))
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  ENDIF
                ENDIF
              ENDDO !njj2
            ENDDO !njj1
            IF(TYPE(1:2).EQ.'XG'.OR.TYPE(1:5).EQ.'XG&YG') THEN
              CALL XZ(ITYP10(nr),XG(1,1),Z)
              WRITE(OP_STRING,
     '          '(1P,'' Rect. cart. coords:    '',3E17.9)')
     '          (Z(nj),nj=1,NJT)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
            WRITE(CHAR1,'(I2)') NUMYGCMPTS
            WRITE(OP_STRING,'('' YG(1..,ng='',I4,'',ne='',I3,'')='','
     '        //CHAR1(1:2)//'D12.4)')
     '        ng,ne,(YG(niyg,ng,ne),niyg=1,NUMYGCMPTS)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO !no_nglist
        ENDIF !type
      ENDDO !nolist

      CALL EXITS('OPGAUSX')
      RETURN
 9999 CALL ERRORS('OPGAUSX',ERROR)
      CALL EXITS('OPGAUSX')
      RETURN 1
      END


