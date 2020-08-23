      SUBROUTINE CM_RANGE(id,nb_type,NBH,NBJ,NEELEM,NHE,NHP,NKH,NKHE,
     '  NKJE,NPF,NPNE,NPNODE,NRLIST,NVHE,NVHP,NVJE,NW,nx,NYNE,NYNP,
     '  CURVCORRECT,PG,SE,TYPE,XA,XE,XG,XP,YG,YP,ZA,ZE,ZG,ZP,ERROR,*)

C#### Subroutine: CM_RANGE
C###  Description:
C###    CM_RANGE finds maximum range of field or dependent variables.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'scal00.cmn'
!     Parameter List
      INTEGER id,nb_type,NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NHE(NEM),NHP(NPM,0:NRM),
     '  NKJE(NKM,NNM,NJM,NEM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKHE(NKM,NNM,NHM,NEM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NRLIST(0:NRM),NVHE(NNM,NBFM,NHM,NEM),
     '  NVHP(NHM,NPM,NCM,0:NRM),NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),nx,
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CURVCORRECT(2,2,NNM,NEM),PG(NSM,NUM,NGM,NBM),
     '  SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),XP(NKM,NVM,NJM,NPM),
     '  YG(NIYGM,NGM,NEM),YP(NYM,NIYM),ZA(NAM,NHM,NCM,NEM),
     '  ZE(NSM,NHM),ZG(NHM,NUM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),TYPE*(*)
!     Local Variables
      INTEGER i,j,nb,NBJ_temp(12),nc,ne,ng,nj,
     '  nj_field,njj1,njj2,noelem,np,
     '  npp,nr,nrr
      REAL*8 DXIX(3,3),ZVAL

      CALL ENTERS('CM_RANGE',*9999)
      nc=1 !Temporary AJP 19-12-91
      nr=NRLIST(1)

      DO i=1,3
        DO j=1,3
          DXIX(i,j)=0.0d0
        ENDDO
      ENDDO

      IF(TYPE(1:5).EQ.'FIELD') THEN          !field variable problem
!news MPN 13-Apr-95: now correctly selects field var from NJ_LOC
!       find range with nodes
        nj_field=NJ_LOC(NJL_FIEL,id,nr)
        ZMINI=XP(1,1,nj_field,NPNODE(1,nr))
        ZMAXI=XP(1,1,nj_field,NPNODE(1,nr))
        DO np=2,NPT(nr)
          ZVAL=XP(1,1,nj_field,np)
          IF(ZVAL.LT.ZMINI) ZMINI=ZVAL
          IF(ZVAL.GT.ZMAXI) ZMAXI=ZVAL
        ENDDO
!       extend range with Gauss points
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          DO njj1=1,3   !Loop over geometry, fibres and field
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              CALL ASSERT(nj.LE.12,'>>ERROR: increase size of NBJ_temp '
     '          //'to NJM',ERROR,*9999)
              IF(nb_type.EQ.0) THEN
                NBJ_temp(nj)=NBJ(nj,ne)
              ELSE
                NBJ_temp(nj)=nb_type
              ENDIF
            ENDDO !njj2
          ENDDO !njj1
          CALL XPXE(NBJ_temp,NKJE(1,1,1,ne),NPF(1,1),
     '      NPNE(1,1,ne),nr,NVJE(1,1,1,ne),SE(1,1,ne),
     '      XA(1,1,ne),XE,XP,ERROR,*9999)
          nb=NBJ_temp(nj_field)
          DO ng=1,NGT(nb)
            CALL XEXG(NBJ_temp,ng,nr,PG,XE,XG,ERROR,*9999)
            ZVAL=XG(nj_field,1)
            IF(ZVAL.LT.ZMINI) ZMINI=ZVAL
            IF(ZVAL.GT.ZMAXI) ZMAXI=ZVAL
          ENDDO !ng
        ENDDO !noelem

      ELSE IF(TYPE(1:9).EQ.'DEPENDENT') THEN !dependent variable problem
!       find range with nodes
        ZMINI=RMAX
        ZMAXI=-RMAX

        DO nrr=1,NRLIST(0)
          nr=NRLIST(nrr)
          CALL YPZP(1,NBH,NEELEM,NHE,NHP(1,nr),NKH(1,1,1,nr),NPNODE,
     '      nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
          DO npp=1,NPNODE(0,nr)
            np=NPNODE(npp,nr)
            ZVAL=ZP(1,1,id,np,nc)
            IF(ZVAL.LT.ZMINI) ZMINI=ZVAL
            IF(ZVAL.GT.ZMAXI) ZMAXI=ZVAL
          ENDDO
!         extend range with Gauss points
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            CALL ZPZE(NBH(1,1,ne),1,NHE(ne),NKHE(1,1,1,ne),NPF(1,1),
     '        NPNE(1,1,ne),nr,NVHE(1,1,1,ne),NW(ne,1),nx,
     '        CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),ZE,ZP,
     '        ERROR,*9999)
            nb=NBH(id,1,ne)
            DO ng=1,NGT(nb)
              CALL ZEZG(1,NBH(1,1,ne),ng,NHE(ne),nx,DXIX,PG,ZE,ZG,
     '          ERROR,*9999)
              ZVAL=ZG(id,1)
              IF(ZVAL.LT.ZMINI) ZMINI=ZVAL
              IF(ZVAL.GT.ZMAXI) ZMAXI=ZVAL
            ENDDO
          ENDDO
        ENDDO

      ELSE IF(TYPE(1:5).EQ.'GAUSS') THEN !Gauss variable problem
        ZMINI=YG(id,1,NEELEM(1,1))
        ZMAXI=YG(id,1,NEELEM(1,1))
        DO noelem=2,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          DO ng=1,NGT(NBJ(1,ne))
            ZVAL=YG(id,ng,ne)
            IF(ZVAL.LT.ZMINI) ZMINI=ZVAL
            IF(ZVAL.GT.ZMAXI) ZMAXI=ZVAL
          ENDDO
        ENDDO

      ENDIF

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' ZMINI='',E12.3,'' ZMAXI='',E12.3)')
     '    ZMINI,ZMAXI
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF
      ZDIFF=ZMAXI-ZMINI
      IF(DABS(ZDIFF).LT.1.0D-6) ZDIFF=1.0D0

      CALL EXITS('CM_RANGE')
      RETURN
 9999 CALL ERRORS('CM_RANGE',ERROR)
      CALL EXITS('CM_RANGE')
      RETURN 1
      END


