      SUBROUTINE UPAERO(NPL,NPNODE,NYNP,NYNR,DL,XP,YP,STRING,FIX,
     '  ERROR,*)

C#### Subroutine: UPAERO
C###  Description:
C###    UPAERO updates finite element variables from aero paramaters
C###    calculated by UPAERO.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'aero00.cmn'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NPL(5,0:3,NLM),NPNODE(0:NP_R_M,0:NRM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM)
      REAL*8 DL(3,NLM),XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM)
      CHARACTER ERROR*(*),STRING*(MXCH)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
!     Local Variables
      INTEGER IBEG,IEND,nl,no_exit,nonode,no_nynr,np,NP1,NP2,
     '  nr,nx,ny,ny1,NY2
      REAL*8 THETA

      CALL ENTERS('UPAERO',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM update aerofoil
C###  Description:
C###    Updates nodal values for aerofoil calculations. i.e.
C###    Update nodal phi values on upstream inflow face.
C###    Update nodal fluxes on downstream outflow face.

        OP_STRING(1)=STRING(1:IEND)
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe23','doc','UPAERO',ERROR,*9999)
      ELSE

        nx=1 !Temporary
        nr=1 !temporary
        CALL ASSERT(CALL_AERO,'aerofoil not defined',ERROR,*9999)

!       Update nodal phi values on upstream inflow face
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
C GMH 8/1/97 Update cmgui link
          CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
          ny=NYNP(1,1,NH_LOC(1,nx),np,0,1,nr)
          IF(FIX(ny,1,nx)) THEN !b.c applied to phi
            IF(DABS(XP(1,1,2,np)).GT.1.0D-6.OR.DABS(XP(1,1,1,np)).
     '        GT.1.0D-6) THEN
              THETA=DATAN2(XP(1,1,2,np),-XP(1,1,1,np))
            ELSE
              THETA=0.d0
            ENDIF
            IF(DOP) write(*,'('' theta='',e13.4)') theta
            YP(ny,1,nx)=CIRCULATION*THETA/(2.d0*PI)
          ENDIF
        ENDDO

!     Update nodal fluxes on downstream outflow face
        DO no_nynr=1,NYNR(0,0,1,nr,nx) !loop over global vars
          ny=NYNR(no_nynr,0,1,nr,nx) !is global var number
          IF(FIX(ny,2,nx)) THEN !b.c applied to phi derivative
            YP(ny,2,nx)=0.0d0
          ENDIF
        ENDDO
        DO no_exit=1,NL_EXIT(0,1)
          nl=NL_EXIT(no_exit,1) !is line# on exit face
          NP1=NPL(2,1,nl) !is 1st node#
          NP2=NPL(3,1,nl) !is 2nd node#
          ny1=NYNP(1,1,NH_LOC(1,nx),np1,0,1,nr)
          ny2=NYNP(1,1,NH_LOC(1,nx),np2,0,1,nr)
          IF(FIX(ny1,2,nx)) THEN !flux is -grad(phi)*0.5d0*Length
C GMH 8/1/97 Update cmgui link
            CALL NODE_CHANGE(NP1,.FALSE.,ERROR,*9999)
            YP(ny1,2,nx)=YP(ny1,2,nx)-0.5d0*DL(3,nl)*REF_VELOC
          ENDIF
          IF(FIX(ny2,2,nx)) THEN !flux is -grad(phi)*0.5d0*Length
C GMH 8/1/97 Update cmgui link
            CALL NODE_CHANGE(NP2,.FALSE.,ERROR,*9999)
            YP(ny2,2,nx)=YP(ny2,2,nx)-0.5d0*DL(3,nl)*REF_VELOC
          ENDIF
        ENDDO

      ENDIF

      CALL EXITS('UPAERO')
      RETURN
 9999 CALL ERRORS('UPAERO',ERROR)
      CALL EXITS('UPAERO')
      RETURN 1
      END


