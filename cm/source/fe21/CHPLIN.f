      SUBROUTINE CHPLIN(ISEG,ISPLIN,CSEG,STRING,ERROR,*)

C#### Subroutine: CHPLIN
C###  Description:
C###    CHPLIN changes polyline segment.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'plin00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISPLIN(NWM,NGRSEGM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER NTLM
      PARAMETER(NTLM=20)
      INTEGER IBEG,IEND,IFROMC,INDEX,INDEX_PLIN,iw,N3CO,nj,nolines,
     '  nopts,NO_SECTIONS,NTPTS
C DPN 17-12-97 - add new arrays to pass through to BEZIER
      INTEGER ISL2BE_LOC(NTLM),ISL3BE_LOC(NTLM),ISN2BE_LOC(NTLM),
     '  ISN3BE_LOC(NTLM)
      REAL*8 XPOINTS(20),YPOINTS(20)
      LOGICAL ABBREV,CBBREV,MOUSE

      CALL ENTERS('CHPLIN',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM change polylines;m
C###  Parameter:      <index ID#[1]>
C###    Specify the identification number of the polyline
C###  Parameter:      <on (WS#s/all)[all]>
C###    Specify the workstation (GX window) to draw the
C###    points on.
C###  Description:
C###    Change polylines using the mouse on the specified workstation.

        OP_STRING(1)=STRING(1:IEND)//';m'
        OP_STRING(2)=BLANK(1:15)//'<index ID[1]>'
        OP_STRING(3)=BLANK(1:15)//'<on (WSS/all)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','CHPLIN',ERROR,*9999)
      ELSE
        CALL CHECKQ('M',noco,1,CO,COQU,STRING,*1)
        MOUSE=.FALSE.
        IF(ABBREV(COQU(noco,1),'M',1)) THEN
          MOUSE=.TRUE.
          IOTYPE=3
        ENDIF

        IF(CBBREV(CO,'INDEX',1,noco+1,NTCO,N3CO)) THEN
          INDEX_PLIN=IFROMC(CO(N3CO+1))
        ELSE
          INDEX_PLIN=1
        ENDIF

        IF(CBBREV(CO,'ON',1,noco+1,NTCO,N3CO)) THEN
          iw=IFROMC(CO(N3CO+1))
        ELSE
          iw=1
        ENDIF

        IF(MOUSE) THEN
          IF(INDEX_PLIN_TYPE(INDEX_PLIN).EQ.1) THEN      !Piecewise linear
          ELSE IF(INDEX_PLIN_TYPE(INDEX_PLIN).EQ.2) THEN !Bezier cubic
            INDEX=1
            NO_SECTIONS=NT_PLIN_SECTIONS(INDEX_PLIN)
            NTPTS=NO_SECTIONS+1
            CALL ASSERT(NTPTS.LE.20,'>>Increase dim of X/YPOINTS',
     '        ERROR,*9999)
            DO nopts=1,NTPTS
              XPOINTS(nopts)=DBLE(PLIN_DATA(1,nopts,INDEX_PLIN))
              YPOINTS(nopts)=DBLE(PLIN_DATA(2,nopts,INDEX_PLIN))
              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
                WRITE(OP_STRING,'('' PLIN_DATA(nj,'',I2,'
     '            //''','',I2,''):'',3E12.3)') nopts,INDEX_PLIN,
     '            (PLIN_DATA(nj,nopts,INDEX_PLIN),nj=1,NJT)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
              ENDIF
            ENDDO
C DPN 17-12-97 - initialise arrays to pass through to BEZIER(), fixing
C                ftnchk error
            DO nolines = 1,NTLM
              ISL2BE_LOC(nolines)=0
              ISL3BE_LOC(nolines)=0
              ISN2BE_LOC(nolines)=0
              ISN3BE_LOC(nolines)=0
            ENDDO
C DPN 17-12-97 - update call to BEZIER()
            CALL BEZIER(INDEX,INDEX_PLIN,ISPLIN(iw,INDEX_PLIN),ISEG,
     '        ISL2BE_LOC,ISL3BE_LOC,ISN2BE_LOC,ISN3BE_LOC,iw,
     '        NO_SECTIONS,XPOINTS,YPOINTS,CSEG,ERROR,*9999)
          ENDIF
        ENDIF
      ENDIF

      CALL EXITS('CHPLIN')
      RETURN
 9999 CALL ERRORS('CHPLIN',ERROR)
      CALL EXITS('CHPLIN')
      RETURN 1
      END


