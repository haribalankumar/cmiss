      SUBROUTINE UPSOUR(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,NBH,NDIPOLES,
     '  NENP,NKH,NP_INTERFACE,NPNODE,
     '  NRLIST,NW,NXLIST,NYNP,CE,DIPOLE_CEN,DIPOLE_DIR,
     '  GD,XG,XP,STRING,ERROR,*)

C#### Subroutine: UPSOUR
C###  Description:
C###    UPSOUR updates the GD (dipole) vector from the sources
C###    defined in ipsour.
C**** Created by Martin Buist

      IMPLICIT NONE

      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'

!     Parameter List
      INTEGER DIPOLE_CEN_NTIME(NDIPOLEM,NRM,NXM),
     '  DIPOLE_DIR_NTIME(NDIPOLEM,NRM,NXM),NBH(NHM,NCM,NEM),
     '  NDIPOLES(NRM,NXM),NENP(NPM,0:NEPM,0:NRM),
     '  NKH(NHM,NPM,NCM,0:NRM),NP_INTERFACE(0:NPM,0:3),
     '  NPNODE(0:NP_R_M,0:NRM),NRLIST(0:NRM),
     '  NW(NEM,3,NXM),NXLIST(0:NXM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CE(NMM,NEM,NXM),DIPOLE_CEN(4,0:NDIPTIMM,NDIPOLEM,NRM,NXM),
     '  DIPOLE_DIR(4,0:NDIPTIMM,NDIPOLEM,NRM,NXM),GD(NZ_GD_M),
     '  XG(NJM,NUM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*),STRING*(MXCH)

!     Local Variables
      INTEGER IBEG,IEND,N3CO,nj,njj2,nonode,no_nrlist,np,nr,
     '  ny1,nx,nxc
      REAL*8 RFROMC,TIME,XPFP(3)
      LOGICAL ALL_REGIONS,CBBREV,ERROR_FLAG

      CALL ENTERS('UPSOUR',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM update source
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Parameter:      <time #[0.0]>
C###    Specify the time variable
C###  Description:
C###    This command will update the GD vector from the dipole
C###    sources defined in ipsour

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(4)=BLANK(1:15)//'<time #[0.0]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe23','doc','UPSOUR',ERROR,*9999)
      ELSE

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)

        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.GT.0,'>>No nx defined for this class',
     '    ERROR,*9999)

        IF(CBBREV(CO,'TIME',2,noco+1,NTCO,N3CO)) THEN
          TIME=RFROMC(CO(N3CO+1))
        ELSE
          TIME=0.0d0
        ENDIF

C MLB  Now done in define sour
C!!!
C LKC 15-APR-2002  Putting the initialisation back which was moved to
C DESOUR on (2001-06-11). I think Martin was doing a coupled FEM/BEM
C problem with no source update, but needed GD initialised.
C My problem is a "static" laplace problem with muliple time steps
C and therefore multiple UPSOUR calls. XPGD assumes that
C GD is initialised to 0.
C NOTE: parrallel directives are not used as there is not enough work
C to justify their use.

C C$OMP   PARALLEL DO
C C$OMP&  PRIVATE(ny1)
C C$OMP&  SHARED(GD)
        DO ny1=1,NZ_GD_M
          GD(ny1)=0.0d0
        ENDDO
C C$OMP   END PARALLEL DO

        ERROR_FLAG=.FALSE.
        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)
          IF(NDIPOLES(nr,nx).GT.0) THEN !Poisson or dipole

C$OMP       PARALLEL DO
C$OMP&      PRIVATE(nj,njj2,nonode,np,XG,XPFP)
C$OMP&      SHARED(CE,DIPOLE_CEN,DIPOLE_CEN_NTIME,DIPOLE_DIR,
C$OMP&        DIPOLE_DIR_NTIME,GD,NBH,NDIPOLES,NENP,NKH,NP_INTERFACE,
C$OMP&        NPNODE,nr,NW,nx,NYNP,TIME,XP)
            DO nonode=1,NPNODE(0,nr)
              IF(.NOT.ERROR_FLAG) THEN
                np=NPNODE(nonode,nr)
                DO njj2=1,NJ_LOC(NJL_GEOM,0,nr)
                  nj=NJ_LOC(NJL_GEOM,njj2,nr)
                  XPFP(njj2)=XP(1,1,nj,np)
                ENDDO !njj2
                CALL XPGD(DIPOLE_CEN_NTIME(1,1,nx),
     '            DIPOLE_DIR_NTIME(1,1,nx),NBH,NDIPOLES(1,nx),NENP,
     '            NKH(1,1,1,nr),NP_INTERFACE,np,nr,NW(1,1,nx),
     '            nx,NYNP,CE(1,1,nx),DIPOLE_CEN(1,0,1,1,nx),
     '            DIPOLE_DIR(1,0,1,1,nx),GD,TIME,XG,XP,XPFP,ERROR,*100)

                GO TO 102
                !this statement is designed to be skipped if no error
                !occurs. However if a error occurs within a subroutine
                !the alternate return jumps to line 100 to set the flag
 100            CONTINUE
C$OMP           CRITICAL(UPSOUR_1)
                ERROR_FLAG=.TRUE.
                WRITE(OP_STRING,'(/'' >>ERROR: An error occurred - '
     '            //'results may be unreliable'')')
                CALL WRITES(IODI,OP_STRING,ERROR,*101)
 101            CONTINUE
C$OMP           END CRITICAL(UPSOUR_1)
 102            CONTINUE
              ENDIF !not error_flag
            ENDDO !nonode
C$OMP       END PARALLEL DO
          ELSE
            WRITE(OP_STRING,'('' >>WARNING: No dipoles to update'
     '        //' in region '',I2)') nr
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF !dipole or Poisson
        ENDDO !nr
      ENDIF

      CALL EXITS('UPSOUR')
      RETURN
 9999 CALL ERRORS('UPSOUR',ERROR)
      CALL EXITS('UPSOUR')
      RETURN 1
      END



