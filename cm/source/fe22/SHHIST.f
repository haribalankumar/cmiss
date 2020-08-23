      SUBROUTINE SHHIST(ISEG,ISHIST,NPLIST,NPNODE,STRING,ERROR,*)

C#### Subroutine: SHHIST
C###  Description:
C###    SHHIST shows history segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISHIST(0:NPM),NPLIST(0:NPM),
     '  NPNODE(0:NP_R_M,0:NRM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,N3CO,nolist,nonode,np,nr
      LOGICAL CBBREV

      CALL ENTERS('SHHIST',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM show history
C###  Parameter:    <at (all/NODE#s)[all]>
C###    specify nodes
C###  Description:
C###    Make the history segments at the specified nodes visible.

        OP_STRING(1)=STRING(1:IEND) //' <at (all/NODE#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','SHHIST',ERROR,*9999)
      ELSE
        IF(CBBREV(CO,'AT',1,noco+1,noco+1,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),NP_R_M,NPLIST(0),NPLIST(1),ERROR,*9999)
        ELSE
          NPLIST(0)=0
          DO nr=1,NRT
            DO nonode=NPLIST(0)+1,NPLIST(0)+NPNODE(0,nr)
              NPLIST(nonode)=NPNODE(nonode,nr)
            ENDDO
            NPLIST(0)=NPLIST(0)+NPNODE(0,nr)
          ENDDO
        ENDIF
C LKC 2-MAY-1998 use iw
C        IF(IWKS(10).GT.0) THEN
C          CALL ACWK(10,1,ERROR,*9999)
        iw=10
        IF(IWKS(iw).GT.0) THEN
          CALL ACWK(iw,1,ERROR,*9999)
          DO nolist=1,NPLIST(0)
            np=NPLIST(nolist)
            IF(DOP) THEN
              WRITE(OP_STRING,'('' np='',I2,'' ISHIST(np)='',I3,'
     '          //''' ISEG='',I3)') np,ISHIST(np),ISEG(ISHIST(np))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            IF(ISEG(ISHIST(np)).EQ.1) THEN
              CALL VISIB(iw,ISEG,ISHIST(np),'VISIBLE',ERROR,*9999)
            ELSE IF(ISEG(ISHIST(np)).EQ.0) THEN
              WRITE(OP_STRING,
     '          '('' >>History is not defined at node '',I4)') np
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO
C LKC 2-MAY-1998 use iw
C          CALL DAWK(10,1,ERROR,*9999)
          CALL DAWK(iw,1,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('SHHIST')
      RETURN
 9999 CALL ERRORS('SHHIST',ERROR)
      CALL EXITS('SHHIST')
      RETURN 1
      END


