      SUBROUTINE EVLAPL(IBT,NBH,NENP,NPLIST3,NPNE,NXI,NXLIST,
     '  LAPL,LAPLSQR,XP,STRING,ERROR,*)

C#### Subroutine: EVLAPL
C###  Description:
C##       Evaluates the surface Laplacian matrix for regularisation
C##     of inverse problems.

C*** Created by Greg Sands 28-NOV-2001

      IMPLICIT NONE

      INCLUDE 'mxch.inc'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),
     '  NBH(NHM,NCM,NEM),NENP(NPM,0:NEPM,0:NRM),NPLIST3(0:NP_R_M),
     '  NPNE(NNM,NBFM,NEM),NXI(-NIM:NIM,0:NEIM,0:NEM),NXLIST(0:NXM)
      REAL*8 LAPL(NY_TRANSFER_M,NY_TRANSFER_M),
     '  LAPLSQR(NY_TRANSFER_M,NY_TRANSFER_M),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*),STRING*(MXCH)

!     Local Variables
      INTEGER IBEG,IEND,nx,nxc

      CALL ENTERS('EVLAPL',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
C#### Command: FEM evaluate laplacian
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe29','doc','EVLAPL',ERROR,*9999)
      ELSE

        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)

        CALL CALC_LAPL(IBT,NBH,NENP,NPLIST3,NPNE,nx,NXI,
     '    LAPL,LAPLSQR,XP,ERROR,*9999)

C LKC 15-JUL-2002 new check
        CALL_CALC_LAPL=.TRUE.

      ENDIF


      CALL EXITS('EVLAPL')
      RETURN
 9999 CALL ERRORS('EVLAPL',ERROR)
      CALL EXITS('EVLAPL')

      RETURN 1
      END


