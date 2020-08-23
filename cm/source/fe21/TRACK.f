      SUBROUTINE TRACK(NEELEM,NLL,NRLIST,NXLIST,DL,STRING,ERROR,*)

C#### Subroutine: TRACK
C###  Description:
C###    TRACK tracks current lines through solution domain.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NEELEM(0:NE_R_M,0:NRM),NLL(12,NEM),
     '  NRLIST(0:NRM),NXLIST(0:NXM)
      REAL*8 DL(3,NLM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,N3CO,NPSTART,nr,NTRL,nx,nxc
      REAL*8 XPFP(3)
      CHARACTER FILE*(MXCH),TYPE*12
      LOGICAL ABBREV,ALL_REGIONS,CBBREV,FIRST,OPFILE

      CALL ENTERS('TRACK',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

        OP_STRING(1)=STRING(1:IEND)//';<FILENAME> current'
        OP_STRING(2)=BLANK(1:15)//'<from NODE#>[1] '
        OP_STRING(3)=BLANK(1:15)//'<region #/ALL>[1]'
        OP_STRING(4)=BLANK(1:15)//'<class #>[1]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

        OP_STRING(1)=STRING(1:IEND)//';<FILENAME> current'
        OP_STRING(2)=BLANK(1:15)//'<from domain point X,Y,Z>'
        OP_STRING(3)=BLANK(1:15)//'<region #/ALL>[1]'
        OP_STRING(4)=BLANK(1:15)//'<class #>[1]'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','TRACK',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN
          FILE=COQU(noco,1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          OPFILE=.TRUE.
          IO4=IOFILE1
          CALL OPENF(IO4,'DISK',FILE(IBEG:IEND)//'.track','NEW',
     '     'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          IO4=IOOP
          OPFILE=.FALSE.
        ENDIF

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)

C CPB 22/11/94 Adding NX_LOC
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)

        IF(CBBREV(CO,'DOMAIN',2,noco+1,NTCO,N3CO)) THEN
          TYPE='DOMAIN_POINT'
        ELSE
          TYPE='NODE_POINT'
        ENDIF
        IF(TYPE(1:12).EQ.'DOMAIN_POINT') THEN
          IF(ABBREV(CO(N3CO+1),'POINT',1)) THEN
            CALL PARSRL(CO(N3CO+2),6,NTRL,XPFP,ERROR,*9999)
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,*)' track from',XPFP(1),XPFP(2),XPFP(3)
             CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
          ENDIF
          FIRST=.FALSE. !The first point is not a node.

        ELSE IF(TYPE(1:10).EQ.'NODE_POINT') THEN
          IF(CBBREV(CO,'FROM',1,noco+1,NTCO,N3CO)) THEN
            NPSTART=IFROMC(CO(N3CO+1))
          ELSE
            NPSTART=1
          ENDIF
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,*)' track from node ',NPSTART
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
          FIRST=.TRUE.

        ENDIF
        nr=NRLIST(1) !Needs extending
        CALL RUNGE_KUTTA_CURRENT_LINES(NEELEM,NLL,nr,DL,FIRST,ERROR,
     '    *9999)
        IF(OPFILE) THEN
          CALL CLOSEF(IO4,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('TRACK')
      RETURN
 9999 CALL ERRORS('TRACK',ERROR)
      IF(OPFILE) CLOSE(UNIT=IO4)
      CALL EXITS('TRACK')
      RETURN 1
      END


