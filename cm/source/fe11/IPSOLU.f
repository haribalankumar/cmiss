      SUBROUTINE IPSOLU(nr,nx,ERROR,*)

C#### Subroutine: IPSOLU
C###  Description:
C###    IPSOLU defines additional solution parameters for region nr.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'solv00.cmn'
!     Parameter List
      INTEGER nr,nx
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ICHAR,INFO,NOQUES
      LOGICAL FILEIP


      CALL ENTERS('IPSOLU',*9999)

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      IF(ITYP9(nr,nx).EQ.3) THEN !multigrid acceleration
        IDEFLT(1)=2
        FORMAT='($,'' Enter #multigrid relaxations on V-cycle '
     '    //'descent [2]: '',I2)'
        IF(IOTYPE.EQ.3) IDATA(1)=Relax1(nx)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,5,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) Relax1(nx)=IDATA(1)

        IDEFLT(1)=2
        FORMAT='($,'' Enter #multigrid relaxations on V-cycle '
     '    //' ascent [2]: '',I2)'
        IF(IOTYPE.EQ.3) IDATA(1)=Relax2(nx)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,5,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) Relax2(nx)=IDATA(1)
      ENDIF !ktyp9

      IF(ITYP5(nr,nx).NE.3) THEN !Not modal
C       Explicit finite differences
        IF((ITYP4(nr,nx).EQ.3).AND.(ITYP16(nr,nx).EQ.4)) THEN
          SPARSEGKK(nx)=2 !row colum sparseness for explicit FD
C       Not explicit finite differences
        ELSE
          ! Using libsolver
          CALL IPSOLU_SOLVER(nr,nx,NOQUES,FILEIP,ERROR,*9999)
        ENDIF
      ENDIF !not modal or explicit finite differences

      IF(ITYP2(nr,nx).EQ.9) THEN !cellular modelling
        FORMAT='('' Specify option for linear solution [0]: '''//
     '    '/''   (0) No output'''//
     '    '/''   (1) Timing and error check'''//
     '    '/''   (2) & Solver output'''//
     '    '/''  *(3) & Global solution matrices'''//
     '    '/''  *(4) & Global stiffness matrices'''//
     '    '/''  *(5) & Element matrices'''//
     '    '/$,''    '',I1)'
        IF(IOTYPE.EQ.3) THEN
           IF(IWRIT4(nr,nx).EQ.-3) IWRIT4(nr,nx)=-IWRIT4(nr,nx)
           IDATA(1)=IWRIT4(nr,nx)
        ENDIF
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,3,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
           IWRIT4(nr,nx)=IDATA(1)
           IF(IWRIT4(nr,nx).EQ.3) IWRIT4(nr,nx)=-IWRIT4(nr,nx)
        ENDIF
      ELSE
        FORMAT='('' Specify option for linear solution [0]: '''//
     '    '/''   (0) No output'''//
     '    '/''   (1) Timing and error check'''//
     '    '/''   (2) & Solver output'''//
     '    '/''   (3) & Global solution matrices'''//
     '    '/''   (4) & Global stiffness matrices'''//
     '    '/''   (5) & Element matrices'''//
     '    '/$,''    '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=IWRIT4(nr,nx)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,5,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) IWRIT4(nr,nx)=IDATA(1)
      ENDIF

      CALL EXITS('IPSOLU')
      RETURN
 9999 CALL ERRORS('IPSOLU',ERROR)
      CALL EXITS('IPSOLU')
      RETURN 1
      END


