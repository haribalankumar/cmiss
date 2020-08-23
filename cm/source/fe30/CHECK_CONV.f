      SUBROUTINE CHECK_CONV(NENQ,NPNODE,NP_INTERFACE,NQNP,NQS,NQXI,NXQ,
     &  NYNP,AQ,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,YP,YQ,STRING,ERROR,*)

C#### Subroutine: CHECK_CONV
C###  Description:
C###    This routine does a subset of 'list variables grid'.
C###    It checks for differences between YQ and YP on endocardial
C###    or epicardial surfaces and writes the differences to a
C###    common block to be accessed by the ITERATE command.
C***  Created by Martin Buist, February 1999

      IMPLICIT NONE

      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'nqloc00.inc'

!     Parameter List
      INTEGER NENQ(0:8,NQM),NPNODE(0:NP_R_M,0:NRM),
     &  NP_INTERFACE(0:NPM,0:3),NQNP(NPM),NQS(NEQM),NQXI(0:NIM,NQSCM),
     &  NXQ(-NIM:NIM,0:4,0:NQM,NAM),NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 AQ(NMAQM,NQM),CQ(NMM,NQM,NXM),DNUDXQ(3,3,NQM),
     &  DXDXIQ(3,3,NQM),DXDXIQ2(3,3,NQM),YP(NYM,NIYM,NXM),
     &  YQ(NYQM,NIQM,NAM,NXM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER ERR,i,IBEG,IBEG1,IBEG2,IEND,IEND1,IEND2,IFROMC,N3CO,nh,
     '  niqV,np,npp,nq,nrb,nrg,nxb,nxc,nxg,nx_upd,nyf,nyp
      REAL*8 FLUX,FLUXDIFF,POTEDIFF
      CHARACTER CHAR1*255,CHAR2*255
      LOGICAL CBBREV,ERROR_FLAG,FOUND
C MHT 24-03-00 removed, not used
C      LOGICAL FIXQ(NYQM,NIYFIXM,NXM)
C      CHARACTER CFROMR*12

      CALL ENTERS('CHECK_CONV',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM check convergence
C###  Description:
C###    Check the differences in the solutions of an extracellular
C###    potential problem and a boundary element problem when
C###    iterating to match potentials and fluxes on epicardial
C###    and endocardial boundaries.
C###  Parameter: <grid_region #>[2]
C###    Specify the region which contains the extracellular
C###    solution.
C###  Parameter: <bem_region #>[1]
C###    Specify the region which contains the boundary element
C###    solution.
C###  Parameter: <grid_class #s>[2]
C###    Specify the class of the extracellular potential problem
C###  Parameter: <grid_update_class #s>[3]
C###    Specify the class of the extracellular update problem
C###  Parameter: <bem_class #>[1]
C###    Specify the class of the boundary element problem.
C###  Parameter: <return (name,name)>[POTE,FLUX]
C###    Specify the return of potential differences to a user defined
C###    variable called POTE and the return of flux differences to a
C###    user defined variable called FLUX.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<grid_region #>[2]'
        OP_STRING(3)=BLANK(1:15)//'<bem_region #>[1]'
        OP_STRING(4)=BLANK(1:15)//'<grid_class #>[2]'
        OP_STRING(5)=BLANK(1:15)//'<grid_update_class #>[3]'
        OP_STRING(6)=BLANK(1:15)//'<bem_class #>[1]'
        OP_STRING(7)=BLANK(1:15)//'<return (name name)>[POTE FLUX]'
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe30','doc','CHECK_CONV',ERROR,*9999)
      ELSE
        !parse class information
        IF(CBBREV(CO,'GRID_CLASS',6,noco+1,NTCO,N3CO)) THEN
          nxc=IFROMC(CO(N3CO+1))
        ELSE
          nxc=2
        ENDIF
        CALL NX_LOC(NX_INQUIRE,nxc,nxg,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nxg.GT.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)

        IF(CBBREV(CO,'BEM_CLASS',5,noco+1,NTCO,N3CO)) THEN
          nxc=IFROMC(CO(N3CO+1))
        ELSE
          nxc=1
        ENDIF
        CALL NX_LOC(NX_INQUIRE,nxc,nxb,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nxb.GT.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)

        IF(CBBREV(CO,'GRID_UPDATE_CLASS',5,noco+1,NTCO,N3CO)) THEN
          nxc=IFROMC(CO(N3CO+1))
        ELSE
          nxc=3
        ENDIF
        CALL NX_LOC(NX_INQUIRE,nxc,nx_upd,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx_upd.GT.0,
     '   '>>No nx defined for this solve class',ERROR,*9999)

        !parse region information
        IF(CBBREV(CO,'GRID_REGION',6,noco+1,NTCO,N3CO)) THEN
          nrg=IFROMC(CO(N3CO+1))
        ELSE
          nrg=2
        ENDIF
        CALL ASSERT(nrg.GT.0,'>>Invalid grid region number',
     '    ERROR,*9999)

        IF(CBBREV(CO,'BEM_REGION',5,noco+1,NTCO,N3CO)) THEN
          nrb=IFROMC(CO(N3CO+1))
        ELSE
          nrb=1
        ENDIF
        CALL ASSERT(nrb.GT.0,'>>Invalid bem region number',
     '    ERROR,*9999)

        IF(CBBREV(CO,'RETURN',3,noco+1,NTCO,N3CO)) THEN
          CHAR1=CO(N3CO+1)
          CHAR2=CO(N3CO+2)
        ELSE
C LKC 6-DEC-2000 Zero length strings no allowed
C          CHAR1=''
C          CHAR2=''
          CHAR1='-'
          CHAR2='-'
        ENDIF
        CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
        CALL STRING_TRIM(CHAR2,IBEG2,IEND2)

        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqV,NIQ_V,ERROR,*9999)
        CALL ASSERT(.NOT.UP_NQNP,'>>nq to np mapping not yet set up',
     '    ERROR,*9999)
        CALL ASSERT(.NOT.UP_NENQ,'>>nq to ne mapping not yet set up',
     '    ERROR,*9999)

        !initialise variables
        YPYQFLUX=0.0d0
        YPYQPOTE=0.0d0
        ERROR_FLAG=.FALSE.

        !calculate differences
C$OMP PARALLEL DO
C$OMP&PRIVATE(FLUX,FLUXDIFF,FOUND,i,nh,np,npp,nq,nyf,nyp,POTEDIFF)
C$OMP&SHARED(CQ,DNUDXQ,DXDXIQ,DXDXIQ2,ERROR_FLAG,NENQ,
C$OMP&  NH_LOC,niqV,NPNODE,
C$OMP&  NP_INTERFACE,NQNP,NQS,NQXI,nrb,nrg,nxb,nxg,NXQ,nx_upd,NYNP,YP,
C$OMP&  YPYQFLUX,YPYQPOTE,YQ)
        DO npp=1,NPNODE(0,nrg)
          IF(.NOT.ERROR_FLAG) THEN
            np=NPNODE(npp,nrg)
            FOUND=.FALSE.
            DO i=1,NP_INTERFACE(np,0)
              IF(nrb.EQ.NP_INTERFACE(np,i)) FOUND=.TRUE.
            ENDDO

            !Remove excluded nodes from the calculation
            DO i=1,CPLST(0,1)
              IF(np.EQ.CPLST(i,1)) FOUND=.FALSE.
            ENDDO

            IF(FOUND) THEN
              nh=NH_LOC(1,nxb)
              nyf=NYNP(1,1,nh,np,0,2,nrb)
              nyp=NYNP(1,1,nh,np,0,1,nrb)
              nq=NQNP(np)
C              BLOOD=.FALSE.
C MHT 24-03-00 if using FIXQ, include in parameter list and multiproc.
C              IF(FIXQ(nq,3,nx_upd)) BLOOD=.TRUE.
            
              CALL GGRADPHIQDN(NENQ,niqV,nq,NQS,NQXI,NXQ,AQ,
     &          CQ(6,nq,nxg),DNUDXQ,DXDXIQ,DXDXIQ2,FLUX,YQ(1,1,1,nxg),
     &          ERROR,*100)

              FLUXDIFF=FLUX+YP(nyf,1,nxb)
              POTEDIFF=YQ(nq,niqV,1,nxg)-YP(nyp,1,nxb)

C$OMP CRITICAL(CHECK_CONV_1)
              IF(DABS(FLUXDIFF).GE.DABS(YPYQFLUX)) YPYQFLUX=
     '          DABS(FLUXDIFF)
              IF(DABS(POTEDIFF).GE.DABS(YPYQPOTE)) YPYQPOTE=
     '          DABS(POTEDIFF)
C$OMP END CRITICAL(CHECK_CONV_1)
            ENDIF
            GO TO 102
            !this statement is designed to be skipped if no error
            !occurs. However if a error occurs within a subroutine
            !the alternate return jumps to line 100 to set the flag
 100        CONTINUE
C$OMP CRITICAL(CHECK_CONV_2)
            ERROR_FLAG=.TRUE.
            WRITE(OP_STRING,'(/'' >>ERROR: An error occurred - '
     '        //'results may be unreliable'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*101)
 101        CONTINUE
C$OMP END CRITICAL(CHECK_CONV_2)
 102        CONTINUE
          ENDIF !not error_flag
        ENDDO
C$OMP END PARALLEL DO
      ENDIF

      IF(IBEG1.LT.IEND1) THEN
        CALL SET_USER_DOUBLE(CHAR1(IBEG1:IEND1),YPYQPOTE,ERR)
        IF(ERR.NE.0) THEN
          ERROR='Unable to set user var "'//CHAR1(IBEG1:IEND1)//'"'
          GO TO 9999
        ENDIF
      ENDIF
      IF(IBEG2.LT.IEND2) THEN
        CALL SET_USER_DOUBLE(CHAR2(IBEG2:IEND2),YPYQFLUX,ERR)
        IF(ERR.NE.0) THEN
          ERROR='Unable to set user var "'//CHAR2(IBEG2:IEND2)//'"'
          GO TO 9999
        ENDIF
      ENDIF

      IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(CHECK_CONV_3)
        WRITE(OP_STRING,'('' Maximum potential difference '',F12.6)')
     '    YPYQPOTE
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Maximum flux difference '',F12.6)')
     '    YPYQFLUX
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(CHECK_CONV_3)
      ENDIF

      CALL EXITS('CHECK_CONV')
      RETURN
 9999 CALL ERRORS('CHECK_CONV',ERROR)
      CALL EXITS('CHECK_CONV')
      RETURN 1
      END
