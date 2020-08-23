      SUBROUTINE IPSOLU_SUPERLU(SLU_PARAM,NOQUES,FILEIP,ERROR,*)

C#### Subroutine: IPSOLU_SUPERLU
C###  Description:
C###    IPSOLU_SUPERLU gets solution parameters for the SuperLU solvers

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'solv00.cmn'
      INCLUDE 'solver.inc'
!     Parameter List
      INTEGER NOQUES,SLU_PARAM(10)
      CHARACTER ERROR*(*)
      LOGICAL FILEIP
!     Local Variables
      INTEGER INFO,ID,IS,SLU_DEFLT(10)
      LOGICAL MULTI


      CALL ENTERS('IPSOLU_SUPERLU',*9999)

C     Create the default array, and optionally reset the settings to
C     their default values.
      CALL SUPERLU_RESET(SLU_DEFLT,ERROR,*9999)
      IF(IOTYPE.NE.3) CALL SUPERLU_RESET(SLU_PARAM,ERROR,*9999)

C     Multi flags if we have an OpenMP executable
      MULTI=.FALSE.
C$    MULTI=.TRUE.

C     Start reading in the values ....
      CALL SUPERLU_GETPARAM('COLUMN ORDERING',ID,SLU_DEFLT,ERROR,*9999)
      CALL SUPERLU_GETPARAM('COLUMN ORDERING',IS,SLU_PARAM,ERROR,*9999)
      IDEFLT(1)=ID+1
      IF(IOTYPE.EQ.3) IDATA(1)=IS+1
      WRITE(FORMAT,'(A,I1,A)')
     '  '('' Specify the column ordering [',IDEFLT(1),']: '''//
     '  '/''   (1) Natural ordering'''//
     '  '/''   (2) Minimum degree ordering from A^T.A'''//
     '  '/''   (3) Minimum degree ordering from A+A^T'''//
     '  '/''   (4) Column Minimum degree ordering'''//
     '  '/$,''    '',I1)'
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '  1,ADATA,ADEFLT,CDATA,CDEFLT,0,IDATA,IDEFLT,1,4,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) CALL SUPERLU_SETPARAM('COLUMN ORDERING',
     '  IDATA(1)-1,SLU_PARAM,ERROR,*9999)


      CALL SUPERLU_GETPARAM('PANEL SIZE',ID,SLU_DEFLT,ERROR,*9999)
      CALL SUPERLU_GETPARAM('PANEL SIZE',IS,SLU_PARAM,ERROR,*9999)
      IDEFLT(1)=ID
      IF(IOTYPE.EQ.3) IDATA(1)=IS
      FORMAT='($,'' Enter the panel size [20] '',I6)'
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '  1,ADATA,ADEFLT,CDATA,CDEFLT,0,IDATA,IDEFLT,1,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) CALL SUPERLU_SETPARAM('PANEL SIZE',IDATA(1),
     '  SLU_PARAM,ERROR,*9999)


      CALL SUPERLU_GETPARAM('RELAX PARAMETER',ID,SLU_DEFLT,ERROR,*9999)
      CALL SUPERLU_GETPARAM('RELAX PARAMETER',IS,SLU_PARAM,ERROR,*9999)
      IDEFLT(1)=ID
      IF(IOTYPE.EQ.3) IDATA(1)=IS
      FORMAT='($,'' Enter the relaxation parameter [5]: '',I6)'
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '  1,ADATA,ADEFLT,CDATA,CDEFLT,0,IDATA,IDEFLT,1,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) CALL SUPERLU_SETPARAM('RELAX PARAMETER',IDATA(1),
     '  SLU_PARAM,ERROR,*9999)


      CALL SUPERLU_GETPARAM('SUPERNODE SIZE',ID,SLU_DEFLT,ERROR,*9999)
      CALL SUPERLU_GETPARAM('SUPERNODE SIZE',IS,SLU_PARAM,ERROR,*9999)
      IDEFLT(1)=ID
      IF(IOTYPE.EQ.3) IDATA(1)=IS
      FORMAT='($,'' Enter the maximum supernode size [100]: '',I6)'
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '  1,ADATA,ADEFLT,CDATA,CDEFLT,0,IDATA,IDEFLT,1,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) CALL SUPERLU_SETPARAM('SUPERNODE SIZE',IDATA(1),
     '  SLU_PARAM,ERROR,*9999)


      CALL SUPERLU_GETPARAM('MIN ROW DIM',ID,SLU_DEFLT,ERROR,*9999)
      CALL SUPERLU_GETPARAM('MIN ROW DIM',IS,SLU_PARAM,ERROR,*9999)
      IDEFLT(1)=ID
      IF(IOTYPE.EQ.3) IDATA(1)=IS
      FORMAT='($,'' Enter the minimum row dimension [400]: '',I6)'
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '  1,ADATA,ADEFLT,CDATA,CDEFLT,0,IDATA,IDEFLT,1,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) CALL SUPERLU_SETPARAM('MIN ROW DIM',IDATA(1),
     '  SLU_PARAM,ERROR,*9999)


      CALL SUPERLU_GETPARAM('MIN COL DIM',ID,SLU_DEFLT,ERROR,*9999)
      CALL SUPERLU_GETPARAM('MIN COL DIM',IS,SLU_PARAM,ERROR,*9999)
      IDEFLT(1)=ID
      IF(IOTYPE.EQ.3) IDATA(1)=IS
      FORMAT='($,'' Enter the minimum column dimension [50]: '',I6)'
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '  1,ADATA,ADEFLT,CDATA,CDEFLT,0,IDATA,IDEFLT,1,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) CALL SUPERLU_SETPARAM('MIN COL DIM',IDATA(1),
     '  SLU_PARAM,ERROR,*9999)


      CALL SUPERLU_GETPARAM('ESTIMATED FILL',ID,SLU_DEFLT,ERROR,*9999)
      CALL SUPERLU_GETPARAM('ESTIMATED FILL',IS,SLU_PARAM,ERROR,*9999)
      IDEFLT(1)=ID
      IF(IOTYPE.EQ.3) IDATA(1)=IS
      FORMAT='($,'' Enter the estimated fill [10]: '',I6)'
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '  FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,0,IDATA,
     '  IDEFLT,-IMAX,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '  INFO,ERROR,*9999)
C     For single processor code we force a positive estimated fill
      IF(.NOT.MULTI.AND.IOTYPE.NE.3) THEN
        CALL SUPERLU_SETPARAM('ESTIMATED FILL',ABS(IDATA(1)),SLU_PARAM,
     '    ERROR,*9999)

C     For multiprocessor code a positive estimated fill gives a fill
C     ratio (the estimated size of L compared to NZA), whilst a negative
C     fill sets the absolute size of L. Note that this is opposite to
C     the internals of the superLU library (and the superLU documentation),
C     but it is done for compatability between single a multi-proc versions
C     of the code. The switching of sign is handled in the libsolver library.

C     Note: for multiproc SuperLU, we really should read and set the
C     estimated L fill, U Fill, and L Sub. In the interests of .ipsolv
C     file compatability, we just set the estimated fill, that in turn
C     sets all three.
      ELSE IF(MULTI.AND.IOTYPE.NE.3) THEN
        CALL SUPERLU_SETPARAM('ESTIMATED FILL',IDATA(1),SLU_PARAM,
     '    ERROR,*9999)
      ENDIF

      CALL EXITS('IPSOLU_SUPERLU')
      RETURN
 9999 CALL ERRORS('IPSOLU_SUPERLU',ERROR)
      CALL EXITS('IPSOLU_SUPERLU')
      RETURN 1
      END


