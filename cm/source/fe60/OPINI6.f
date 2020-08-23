      SUBROUTINE OPINI6(NODENVCB,NPNODE,nr,NVCB,nx,NYNP,
     '  YP,ERROR,*)

C#### Subroutine: OPINI6
C###  Description:
C###    OPINI6 inputs initial conditions and boundary conditions for
C###    FE60 problems.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'maqloc00.inc'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'voro00.inc'
!     Parameter List
      INTEGER NODENVCB(NVCBM),
     '  NPNODE(0:NP_R_M),nr,NVCB(-1:3,NVCBM),nx,
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 YP(NYM,NIYM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER bnvc,np,i,nonode,nhx

      CALL ENTERS('OPINI6',*9999)

      DO bnvc=1,NVBT
        nonode=NODENVCB(bnvc)
        np=NPNODE(nonode)
        WRITE(OP_STRING,'(/$,'' Voronoi B node '',I5,'' (np = '',
     '    I5,'') maps onto node(s) '',3I5)') bnvc,np,
     '    (NPNODE(NVCB(i,bnvc)),i=1,NVCB(0,bnvc))
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        IF(NVCB(BCTYPE,bnvc).EQ.WALL) THEN
          WRITE(OP_STRING,'('' BC Type: Wall'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Velocity components: '',3D12.5)')
     '      (YP(NYNP(1,1,nh_loc(nhx,nx),np,0,1,nr),1),
     '      nhx=1,nh_loc(0,nx)-1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ELSEIF(NVCB(BCTYPE,bnvc).EQ.INLET) THEN
          WRITE(OP_STRING,'('' BC Type: Inlet'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Velocity components: '',3D12.5)')
          WRITE(OP_STRING,'('' Velocity components: '',3D12.5)')
     '      (YP(NYNP(1,1,nh_loc(nhx,nx),np,0,1,nr),1),
     '      nhx=1,nh_loc(0,nx)-1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ELSEIF(NVCB(BCTYPE,bnvc).EQ.OUTLET) THEN
          WRITE(OP_STRING,'('' BC Type: Outlet'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ELSEIF(NVCB(BCTYPE,bnvc).EQ.FREESLIP) THEN
          WRITE(OP_STRING,'('' BC Type: Freeslip'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ELSEIF(NVCB(BCTYPE,bnvc).EQ.DRIVING) THEN
          WRITE(OP_STRING,'('' BC Type: Driving'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

      ENDDO !bnvc

      CALL EXITS('OPINI6')
      RETURN
 9999 CALL ERRORS('OPINI6',ERROR)
      CALL EXITS('OPINI6')
      RETURN 1
      END


