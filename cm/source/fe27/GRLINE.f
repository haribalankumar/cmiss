      SUBROUTINE GRLINE(NLLINE,NLLIST,NPL,NPLIST,NRLIST,STRING,ERROR,*)

C#### Subroutine: GRLINE
C###  Description:
C###    GRLINE groups lines.

C**** NTGRLI is number of line groups currently defined.
C**** LAGRLI(nogrli) is label given to group number NOGRLI.
C**** LIGRLI(0,nogrli) is number in list for group number NOGRLI.
C**** LIGRLI(1..,nogrli) is list for group number NOGRLI.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grou00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'mach00.inc'
!     Parameter List
      INTEGER NLLINE(0:NL_R_M,0:NRM),NLLIST(0:NLM),NPL(5,0:3,NLM),
     &  NPLIST(0:NPM),NRLIST(0:NRM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IBEG2,IEND,IEND2,N1GRLI,N3CO,N3CO1,nextnode,ni,
     &  nl,NLNP(-3:3),node1,node2,nogrli,XIDIRN(6)
      CHARACTER CHAR*30,CHAR2*2,LABEL*30
      LOGICAL ALIEN_BASIS,ALL_REGIONS,CBBREV,FOUND,UNSORTED,XIEND

      CALL ENTERS('GRLINE',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
C        CHAR2=CFROMI(NTGRLI+1,'(I2)')
        WRITE(CHAR2,'(I2)') NTGRLI+1
        CALL STRING_TRIM(CHAR2,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM group lines line#s/LINE_GROUP
C###  Parameter:     <as LABEL[element_1]>
C###    Specifies the name of the line group.
C###  Parameter:     <(sorted/unsorted)[sorted]>
C###    The sorted command sorts the line numbers defined
C###    in the group in ascending order and removes duplicates.
C###  Parameter:     <region (#s/all)[1]>
C###    Specify the line file region numbers to be defined.
C###    The all value specifies all currently defined regions.
C###    Any region numbers contained in the file but are not specified
C###    will be skipped.
C###  Description:
C###    Group lines under a label. By default the list is sorted
C###    by line numbers into ascending order and duplicates are
C###    removed - the UNSORTED option prevents this.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'LINE#s/LINE_GROUP'
        OP_STRING(3)=BLANK(1:15)
     '    //'<as LABEL[element_'//CHAR2(IBEG2:IEND2)//']>'
        OP_STRING(4)=BLANK(1:15)//'<(sorted/unsorted)[sorted]>'
        OP_STRING(5)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM group lines between nodes NODE1#,NODE2#
C###  Parameter:     <as LABEL[element_1]>
C###    Specifies the name of the line group.
C###  Parameter:     <(sorted/unsorted)[sorted]>
C###    The sorted command sorts the line numbers defined
C###    in the group in ascending order and removes duplicates.
C###  Parameter:     <region (#s/all)[1]>
C###    Specify the line file region numbers to be defined.
C###    The all value specifies all currently defined regions.
C###    Any region numbers contained in the file but are not specified
C###    will be skipped.
C###  Description:
C###    Group lines between two nodes under a label. By default the list
C###    is sorted by line numbers into ascending order and duplicates
C###    are removed - the UNSORTED option prevents this.

        OP_STRING(1)=STRING(1:IEND)//'between nodes NODE1#,NODE2#'
        OP_STRING(2)=BLANK(1:15)
     '    //'<as LABEL[element_'//CHAR2(IBEG2:IEND2)//']>'
        OP_STRING(3)=BLANK(1:15)//'<(sorted/unsorted)[sorted]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','GRLINE',ERROR,*9999)
      ELSE
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)

C DMAL 13 JAN 2004 - Adding "fem group lines between nodes" command
C This can be extended to "fem group lines between lines" etc...
        IF(CBBREV(CO,'BETWEEN',2,noco+1,NTCO,N3CO)) THEN
          IF(CBBREV(CO,'NODES',2,N3CO+1,N3CO+1,N3CO1)) THEN
            ! Parse nodes
            CDATA(1)='NODES'
            CALL PARSILG(NPLIST,NPM,CDATA(1),CO(N3CO1+1),ERROR,*9999)
            IF(NPLIST(0).EQ.2) THEN
              node1=NPLIST(1)
              node2=NPLIST(2)
            ELSEIF(NPLIST(0).EQ.4)THEN
              ! To be implemented for grouping nodes on a 2d face
              CALL ASSERT(.FALSE.,
     &          '>> No support for grouping nodes on faces',
     &          ERROR,*9999)
            ELSE
              CALL ASSERT(.FALSE.,
     &        '>> Incorrect number of nodes to group between',
     &        ERROR,*9999)
            ENDIF
            
            DO nl=-3,3
              NLNP(nl)=0
            ENDDO
            
            ! finding line segments coming off node1 in each xi-direction
            DO nl=1,NLT
              IF(NPL(2,1,nl).EQ.node1)THEN
                NLNP(NPL(1,0,nl))=nl
              ELSEIF(NPL(3,1,nl).EQ.node1)THEN
                NLNP(-NPL(1,0,nl))=nl
              ENDIF
            ENDDO
            
            XIDIRN(1)=1
            XIDIRN(2)=-1
            XIDIRN(3)=2
            XIDIRN(4)=-2
            XIDIRN(5)=3
            XIDIRN(6)=-3
            
            ni=1
            FOUND=.FALSE.
            ALIEN_BASIS=.FALSE.
            DO WHILE ((.NOT.FOUND).AND.(ni.LE.6))
              IF(NLNP(XIDIRN(ni)).NE.0)THEN
                XIEND=.FALSE.
                NLLIST(0)=0
                nl=NLNP(XIDIRN(ni))
                IF((NPL(1,1,nl).NE.1).AND.(NPL(1,1,nl).NE.4))THEN
                  ALIEN_BASIS=.TRUE.
                ENDIF
                IF(XIDIRN(ni).GT.0)THEN
                  nextnode=3
                ELSE
                  nextnode=2
                ENDIF
                DO WHILE ((.NOT.FOUND).AND.(.NOT.XIEND))
                  NLLIST(0)=NLLIST(0)+1
                  NLLIST(NLLIST(0))=nl
                  IF(NPL(nextnode,1,nl).EQ.node2)THEN
                    FOUND=.TRUE.
                  ELSE
                    IF((XIDIRN(ni).GT.0).AND.(NPL(3,0,nl).NE.0))THEN
                      nl=NPL(3,0,nl)
                    ELSEIF((XIDIRN(ni).LT.0).AND.(NPL(2,0,nl).NE.0))THEN
                      nl=NPL(2,0,nl)
                    ELSE
                      XIEND=.TRUE.
                    ENDIF
                  ENDIF
                ENDDO
              ENDIF
              ni=ni+1
            ENDDO
            
            IF(.NOT.FOUND)THEN
              NLLIST(0)=0
              WRITE(OP_STRING,'('' >>WARNING: No lines found between'
     '          //' nodes!'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            
            IF(ALIEN_BASIS)THEN
              WRITE(OP_STRING,'('' >>WARNING: only linear and cubic-'
     '          //'Hermite lines are implemented!'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF
            
          ELSE
            CALL ASSERT(.FALSE.,
     &        '>> Bound type not implemented',
     &        ERROR,*9999)
          ENDIF
        ELSE
          CALL PARSE_LINES(NLLINE,NLLIST,noco-1,NRLIST,NTCO,CO,
     '      ERROR,*9999)
        ENDIF

        IF(CBBREV(CO,'UNSORTED',1,noco+1,NTCO,N3CO)) THEN
          UNSORTED=.TRUE.
        ELSE
          UNSORTED=.FALSE.
        ENDIF

        IF(.NOT.UNSORTED) THEN
C         Sort and remove duplicates from the line list
          CALL ILISTRMDUP(NLLIST(0),NLLIST(1),ERROR,*9999)
        ENDIF

C        CALL ASSERT(NLLIST(0).LE.GRLI_MAXDATA,
C     '    '>>Increase array sizes in grou00.cmn',ERROR,*9999)

        IF(CBBREV(CO,'AS',1,noco+1,NTCO,N3CO)) THEN
          !Check whether group name already exists
          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
C          CHAR=CUPPER(CO(N3CO+1)(IBEG:IEND))
          CALL CUPPER(CO(N3CO+1)(IBEG:IEND),CHAR)
          N1GRLI=0
          DO nogrli=1,NTGRLI
C            LABEL=CUPPER(LAGRLI(nogrli))
            CALL CUPPER(LAGRLI(nogrli),LABEL)
            CALL STRING_TRIM(LABEL,IBEG2,IEND2)
            IF(CHAR(IBEG:IEND).EQ.LABEL(IBEG2:IEND2)) THEN
              N1GRLI=nogrli !is existing group label ID
              GO TO 100
            ENDIF
          ENDDO
 100      IF(N1GRLI.EQ.0) THEN !need new group label
            NTGRLI=NTGRLI+1 !is new total #groups
            CALL ASSERT(NTGRLI.LE.GRLI_MAXGRP,
     '        '>>Increase array sizes in grou00.cmn',ERROR,*9999)
            N1GRLI=NTGRLI
            CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
            LAGRLI(N1GRLI)=CO(N3CO+1)(IBEG:IEND) !new line group label
          ENDIF

        ELSE
          NTGRLI=NTGRLI+1 !is new total
          CALL ASSERT(NTGRLI.LE.GRLI_MAXGRP,
     '      '>>Increase array sizes in grou00.cmn',ERROR,*9999)
          N1GRLI=NTGRLI
C          CHAR2=CFROMI(N1GRLI,'(I2)')
          WRITE(CHAR2,'(I2)') N1GRLI
          CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
          LAGRLI(N1GRLI)='line_'//CHAR2(IBEG2:IEND2) !new label
        ENDIF

C TVK 28July1999: dynamic groups (entire routine is new)
          NLIGRLI(N1GRLI)=NLLIST(0)
         CALL ALLOCATE_MEMORY(NLLIST(0),0,INTTYPE,LIGRLI_PTR(N1GRLI),
     '    MEM_INIT,ERROR,*9999)
         CALL ILIST_COPY(NLLIST(0),NLLIST(1),%VAL(LIGRLI_PTR(N1GRLI)))



      ENDIF

      CALL EXITS('GRLINE')
      RETURN
 9999 CALL ERRORS('GRLINE',ERROR)
      CALL EXITS('GRLINE')
      RETURN 1
      END


