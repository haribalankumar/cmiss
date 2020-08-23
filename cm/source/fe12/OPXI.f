      SUBROUTINE OPXI(NEP,NPNODE,XI_TYPE,LD,XID,XIP,ERROR,*)

C#### Subroutine: OPXI
C###  Description:
C###    OPXI outputs Xi coordinates of nodes or data, where XI_TYPE
C###    is 'NODE' for nodes and 'DATA' for data

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'data00.cmn'
!     Parameter List
      INTEGER LD(NDM),NEP(NPM),NPNODE(0:NP_R_M)
      REAL*8 XID(NIM,NDM),XIP(NIM,NPM)
      CHARACTER ERROR*(*),XI_TYPE*4
!     Local Variables
      INTEGER nd,nj,nonode,np

      CALL ENTERS('OPXI',*9999)

C LKC 25-MAY-1998 This code has nothing to do with xi's
C
C      DO nonode=1,NPNODE(0,nr)
C        np=NPNODE(nonode,nr)
C        WRITE(OP_STRING,'('' XP(2,1,ni,'',I5,''): '',3E12.3)')
C     '    np,(XP(2,1,ni,np),ni=1,3)
C        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C      ENDDO

C  LKC 25-MAY-1998 New xi output
      IF(XI_TYPE(1:4).EQ.'DATA') THEN
        DO nd=1,NDT
          WRITE(OP_STRING,
     '      '(''nd ='',I5,''  LD ='',I5,''  XID = '',3E16.7)')
     '      nd,LD(nd),(XID(nj,nd),nj=1,NXIDEF)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
      ELSEIF(XI_TYPE.EQ.'NODE') THEN
        DO nonode=1,NPNODE(0)
          np=NPNODE(nonode)
          WRITE(OP_STRING,
     '      '(''np ='',I5,''  NEP ='',I5,''  XIP = '',3E16.7)')
     '      np,NEP(np),(XIP(nj,np),nj=1,NJT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

      CALL EXITS('OPXI')
      RETURN
 9999 CALL ERRORS('OPXI',ERROR)
      CALL EXITS('OPXI')
      RETURN 1
      END

