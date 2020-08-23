      SUBROUTINE OPHIST(NBH,NEELEM,NHE,NHP,NKH,
     '  nolist,np,NPNODE,nr,NVHP,NYNE,NYNP,
     '  YP,ZA,ZP,OPFILE,TYPE,ERROR,*)

C#### Subroutine: OPHIST
C###  Description:
C###    OPHIST outputs history of solution at a node.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NBH(NHM,NCM,NEM),NEELEM(0:NE_R_M,0:NRM),NHE(NEM),NHP(NPM),
     '  NKH(NHM,NPM,NCM),nolist,np,NPNODE(0:NP_R_M,0:NRM),nr,
     '  NVHP(NHM,NPM,NCM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 YP(NYM,NIYM),ZA(NAM,NHM,NCM,NEM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER TYPE*(*),ERROR*(*)
      LOGICAL OPFILE
!     Local Variables
      INTEGER nc,NOHIST,nx,ny
      REAL*8 T

      CALL ENTERS('OPHIST',*9999)

      nx=1 !Temporary

      NOHIST=0
 250  NOHIST=NOHIST+1
      READ(9,'('' YP(ny_lhs,1) at t='',E11.4,'' :'')',END=350) T
      READ(9,'(10E13.5)') (YP(ny,1),ny=1,NYT(1,1,nx))
      READ(9,'('' YP(ny_rhs,1) at t='',E11.4,'' :'')') T
      READ(9,'(10E13.5)') (YP(ny,1),ny=1,NYT(2,1,nx))
      IF(TYPE(1:5).EQ.'VALUE') THEN
!new MPN 6-Jan-95: current soln now stored in YP(ny,1) for nonlin probs
        CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
     '    nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
!old
c        IF(ITYP6(nr,nx).EQ.1) THEN      !linear
c          CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
c     '      nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
c        ELSE IF(ITYP6(nr,nx).EQ.2) THEN !nonlinear
c          CALL YPZP(4,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
c     '      nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
c        ENDIF
      ELSE IF(TYPE(1:8).EQ.'REACTION') THEN
        CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
     '    nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
      ENDIF

      nc=1 !temporary
      IF(OPFILE)THEN
        IF(nolist.EQ.1) WRITE(3,'(E13.5)') T
        WRITE(OP_STRING,'(E13.5)') ZP(1,1,1,np,nc)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ELSE
        WRITE(OP_STRING,'('' T='',E13.5,'' ZP='',E13.5)')
     '    T,ZP(1,1,1,np,nc)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF
      GO TO 250
 350  CALL CLOSEF(9,ERROR,*9999)

      CALL EXITS('OPHIST')
      RETURN
 9999 CALL ERRORS('OPHIST',ERROR)
      CALL EXITS('OPHIST')
      RETURN 1
      END


