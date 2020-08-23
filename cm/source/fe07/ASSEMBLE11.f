      SUBROUTINE ASSEMBLE11(IBT,IDO,INP,ISC_GD,ISC_GK,ISR_GD,ISR_GK,
     '  nb,NBH,NBJ,NEELEM,NENP,NHE,NKJE,NO_NE,NORD,NPF,NPLIST2,NPNE,nr,
     &  NVJE,nx,NYNE,NYNP,NYNR,NZ_ESED,nzr,CE,CG,CGE,CP,ED,EM,ER,ES,GD,
     '  GK,GK2,GR,PG,RG,SE,STACK_ED,STACK_EM,STACK_ES,WG,XA,XAB,XG,
     &  XP,YG,ZE,ZG,UPDATE_MATRIX,UPDATE_VECTOR,UPDATE_VELOCITY,
     &  RET_ERROR,*)

C#### Subroutine: ASSEMBLE11
C###  Description:
C###    ASSEMBLE11 assembles the global unreduced matrices GK, GD,
C###    for time dependent pulmonary transport problems.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'error0.inc'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'pulm00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISC_GD(NISC_GDM),ISC_GK(NISC_GKM),ISR_GD(NISR_GDM),
     '  ISR_GK(NISR_GKM),nb,NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M),NENP(NPM,0:NEPM),NHE(NEM),
     '  NKJE(NKM,NNM,NJM,NEM),NO_NE(NEM),NORD(5,NE_R_M),NPF(9,NFM),
     '  NPLIST2(0:NPM),NPNE(NNM,NBFM,NEM),nr,nx,NVJE(NNM,NBFM,NJM,NEM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM),NYNR(0:NY_R_M,0:NRCM,NCM),
     '  NZ_ESED(0:12,NE_R_M),nzr
      REAL*8 CE(NMM,NEM),CG(NMM,NGM),CGE(NMM,NGM,NEM),CP(NMM,NPM),
     &  GD(NZ_GD_M),GK(NZ_GK_M),GK2(NZ_GK_M),GR(NYROWM),
     &  PG(NSM,NUM,NGM,NBM),SE(NSM,NBFM,NEM),
     &  STACK_ED(12,NE_R_M),STACK_EM(12,NE_R_M),STACK_ES(12,NE_R_M),
     &  WG(NGM,NBM),XA(NAM,NJM,NEM),XAB(NORM,NEM),XP(NKM,NVM,NJM,NPM),
     &  YG(NIYGM,NGM,NEM),ZE(NSM,NHM),ZG(NHM,NUM),
     &  ED(NHM*NSM,NHM*NSM),EM(NHM*NSM,NHM*NSM),ER(NHM*NSM),
     &  ES(NHM*NSM,NHM*NSM),RG(NGM),XG(NJM,NUM)
      LOGICAL UPDATE_MATRIX,UPDATE_VECTOR,UPDATE_VELOCITY
      CHARACTER RET_ERROR*(*)
      INTEGER i,nb_radius,ne,nhs_cap,nn,noelem,nonode,no_nynr1,np,nv,
     &  ny1,nz
      CHARACTER CHAR*1,ERROR*(ERRSTRLEN)
!     External functions
      INTEGER IDIGITS

      CALL ENTERS('ASSEMBLE11',*9999)

      UPDATE_VELOCITY=.FALSE.
      IF(UPDATE_MATRIX)THEN !previous: IF(.NOT.UPDATE_VELOCITY)THEN
c        CALL INIT_SPARSE_MATRIX(NYNR(0,2,1),NYNR(0,2,1),ISC_GK,ISR_GK,
c     '    0,0,1,NYT(1,1,nx),NZ_GK_M,NZT(1,nx),NYNR(0,1,1),NYNR(0,1,1),
c     '    KTYP24,GK,ERROR,*9999)

        DO nz=1,NZ_GK_M
          GK(nz)=0.d0
        ENDDO !nz
        
        NYT(1,3,nx)=NYT(1,1,nx)
        NYT(2,3,nx)=NYT(2,1,nx)
        NZT(3,nx)=NZT(1,nx)

        IF(NISC_GKM.GT.NISC_GDM) THEN
          WRITE(CHAR,'(I1)') IDIGITS(NISC_GKM)
          WRITE(ERROR,'(''>>Increase NISC_GDM to '',I'//CHAR//')')
     '      NISC_GKM
          GO TO 9999
        ENDIF
        IF(NISR_GKM.GT.NISR_GDM) THEN
          WRITE(CHAR,'(I1)') IDIGITS(NISR_GKM)
          WRITE(ERROR,'(''>>Increase NISR_GDM to '',I'//CHAR//')')
     '      NISR_GKM
          GO TO 9999
        ENDIF
        DO nz=1,NISC_GKM
          ISC_GD(nz)=ISC_GK(nz)
        ENDDO !nz
        DO nz=1,NISR_GKM
          ISR_GD(nz)=ISR_GK(nz)
        ENDDO !nz
        IF(ITYP3(nr,nx).LE.2) THEN !gas transport or water-heat
          CALL INIT_SPARSE_MATRIX(NYNR(0,2,1),NYNR(0,2,1),ISC_GD,ISR_GD,
     '      0,0,1,NYT(1,3,nx),NZ_GD_M,NZT(3,nx),NYNR(0,1,1),NYNR(0,1,1),
     '      KTYP24,GD,ERROR,*9999)
          CALL INIT_SPARSE_MATRIX(NYNR(0,2,1),NYNR(0,2,1),ISC_GK,ISR_GK,
     '      0,0,1,NYT(1,3,nx),NZ_GK_M,NZT(3,nx),NYNR(0,1,1),NYNR(0,1,1),
     '      KTYP24,GK2,ERROR,*9999)
        ENDIF
        DO nonode=1,NPM !initialising
          NPLIST2(nonode)=0
        ENDDO !nonode
        nhs_cap=0 !initialising counter used in assemble11_dynam
        NPLIST2(0)=0
      ENDIF
      IF(UPDATE_VECTOR)THEN
        DO no_nynr1=1,NYNR(0,1,1)
          ny1=NYNR(no_nynr1,1,1)
          GR(ny1)=0.0d0
        ENDDO !no_nynr (ny1)
      ENDIF
      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        NO_NE(ne)=noelem !records the local noelem for global ne
        IF(ITYP3(nr,nx).LT.3.OR.ITYP3(nr,nx).EQ.5) THEN
          nb=NBJ(1,ne)
          nb_radius=NBJ(nj_radius,ne)
          DO nn=1,2
            np=NPNE(nn,nb_radius,ne)
            nv=NVJE(nn,nb_radius,nj_radius,ne)
            IF(DOP)THEN
              WRITE(OP_STRING,'('' np '',I6,'' radius '',F12.4)') np,
     &          XP(1,nv,nj_radius,np)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            CALL ASSERT(XP(1,nv,nj_radius,np).GT.0.d0,'>>Radius is'
     '        // ' zero',ERROR,*9999)
          ENDDO !nn
        ELSE IF(ITYP3(nr,nx).EQ.3)THEN
          IF(CE(nm_Dh,ne).LT.0.d0) THEN
            WRITE(OP_STRING,
     &        '('' Diam.LT.0.d0 for ne: '',I6,'' Dh '',D12.4)') 
     &        ne,CE(nm_Dh,ne)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
        
        CALL ASSEMBLE11_DYNAM(IBT,IDO,INP,ISC_GK,ISR_GK,nb,
     &    NBH(1,1,ne),NBJ(1,ne),ne,NENP,NHE(ne),nhs_cap,NKJE(1,1,1,ne),
     &    NORD,NPF(1,1),NPLIST2,NPNE,nr,NVJE(1,1,1,ne),nx,NYNE,NYNP,
     &    NZ_ESED(0,noelem),nzr,CE(1,ne),CG,CGE(1,1,ne),CP,ED,EM,ER,ES,
     &    GK,GR,PG,rg,SE(1,1,ne),STACK_ED(1,noelem),STACK_EM(1,noelem),
     &    STACK_ES(1,noelem),WG,XA,XAB,XG,XP,YG(1,1,ne),ZE,ZG,
     &    UPDATE_MATRIX,UPDATE_VECTOR,ERROR,*9999)
        
      ENDDO !noelem (ne)
      IF(ITYP3(nr,nx).LE.2.OR.ITYP3(nr,nx).EQ.5)THEN
        DO noelem=1,NEELEM(0)
          DO i=1,NZ_ESED(0,noelem)
            nz=NZ_ESED(i,noelem)
            GK(nz)=GK(nz)+STACK_ES(i,noelem)
            GD(nz)=GD(nz)+STACK_ED(i,noelem)
            GK2(nz)=GK2(nz)+STACK_EM(i,noelem)
          ENDDO !i
        ENDDO !noelem
      ENDIF !ITYP3

      CALL EXITS('ASSEMBLE11')
      RETURN

 9999 CALL ERRORS('ASSEMBLE11',ERROR)
      RET_ERROR=ERROR
      CALL EXITS('ASSEMBLE11')
      RETURN 1
      END



