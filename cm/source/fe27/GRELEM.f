      SUBROUTINE GRELEM(NEELEM,NEL,NELIST,NELIST2,NLLINE,NLLIST,NORD,
     '  NRLIST,NXI,NFF,STRING,ERROR,*)

C#### Subroutine: GRELEM
C###  Description:
C###    GRELEM groups elements.

C**** NTGREL is number of element groups currently defined.
C**** LAGREL(nogrel) is label given to group number NOGREL.
C**** LIGREL(0,nogrel) is number in list for group number NOGREL.
C**** LIGREL(1..,nogrel) is list for group number NOGREL.

C MPN Jul2004: several changes to enable selection of external elements 
C              with options for specified (s1, s2 or s3) boundaries of mesh

C Glenn Ramsey Feb2005: bug fix to S1/2/3 options to correctly exclude
C                       collapsed elements that are not on a boundary

C A Swan Feb2007: added option to create group by specifying elements
C                       to be excluded from the group
      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grou00.cmn'
      INCLUDE 'mach00.inc'
!     Parameter List
      INTEGER NEL(0:NELM,NLM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     &  NELIST2(0:NEM),NLLINE(0:NL_R_M,0:NRM),NLLIST(0:NLM),
     &  NORD(5,NE_R_M),NRLIST(0:NRM),NXI(-NIM:NIM,0:NEIM,0:NEM),
     &  NFF(6,NEM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER ERR,GROUP_TYPE,IBEG,IBEG2,IEND,IEND2,IFROMC,incl,
     &  N3CO,N4CO,ne,ne1,ne_parent,NELIST3(0:NEM),nl,noelem,noelem1,
     &  s_dirn,s_value(3),mxi(-3:3),ae,a_dirn
      CHARACTER CHAR2*2,STRING2*(255)
      LOGICAL ALL_REGIONS,AS,CBBREV,EXCLUDE_ELEM,EXTERNAL,
     &  EXTERNAL_ELEM,FOUND,UNSORTED,s_bdy(3)

      CALL ENTERS('GRELEM',*9999)
      IF(noco.EQ.NTCO) THEN
        CO(noco+1)='?'
        NTCO=NTCO+1
        CALL STRING_TRIM(STRING,IBEG,IEND)
        STRING=STRING(IBEG:IEND)//' elements'
      ENDIF        
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
C        CHAR2=CFROMI(NTGREL+1,'(I2)')
        WRITE(CHAR2,'(I2)') NTGREL+1
        CALL STRING_TRIM(CHAR2,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM group elements ELEM#s/ELEM_GROUP/all_elements
C###  Parameter:     <as LABEL[element_1]>
C###    Specifies the name of the element group.
C###  Parameter:     <(sorted/unsorted)[sorted]>
C###    The sorted command sorts the elements numbers defined
C###    in the group in ascending order and removes duplicates.
C###  Parameter:     <by generation/order/lobe/terminal>
C###    The 'by' command allows grouping in a tree-structured mesh by
C###    generation (recursive from stem element), Horsfield order
C###    (recursive from terminal elements), lobe (using lobe groups
C###    set up in 'define mesh'), terminal branches of tree.
C###  Parameter:     <parent PARENT#[1]>
C###    Grouping by parent (all elements subtending
C###    each in a list of parent elements).
C###  Parameter:     <(external <s(1/2/3)=(0/1)>/all)[all]>
C###    Specifies that only outer elements are to be included
C###    in the group. The s(1/2/3) option restricts the selection
C###    to the specfied direction boundary.  
C###  Parameter:     <region (#s/all)[1]>
C###    Specify the element file region numbers to be defined.
C###    The all value specifies all currently defined regions.
C###    Any region numbers contained in the file but are not specified
C###    will be skipped.
C###  Parameter:     <exclude elements ELEM#s/ELEM_GROUP>
C###    The exclude command removes the specified element/s 
C###    from the group.
C###  Description:
C###    Group elements under a label. By default the list is sorted
C###    by element numbers into ascending order and duplicates are
C###    removed - the UNSORTED option prevents this.

        OP_STRING(1)=STRING(1:IEND)//' ELEM#s/ELEM_GROUP/all_elements'
        OP_STRING(2)=BLANK(1:15)
     '    //'<as LABEL[element_'//CHAR2(IBEG2:IEND2)//']>'
        OP_STRING(3)=BLANK(1:15)//'<(sorted/unsorted)[sorted]>'
        OP_STRING(4)=BLANK(1:15)//'<by generation/order/lobe/terminal>'
        OP_STRING(5)=BLANK(1:15)//'<parent PARENT#[1]>'
        OP_STRING(6)=BLANK(1:15)//'<(external '
     &    //'<s1=(0/1)>[both] <s2=(0/1)>[both] <s3=(0/1)>[both]> '
     &    //'/all)[all]'
        OP_STRING(7)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(8)=BLANK(1:15)//'<exclude elements ELEM#s/ELEM_GROUP>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM group elements with lines LINE#s/LINE_GROUP
C###  Parameter:     <as LABEL[element_1]>
C###    Specifies the name of the element group.
C###  Parameter:     <(external/all)[all]>
C###    Specifies that only outer elements are to be included
C###    in the group.
C###  Parameter:     <s1=(0/1)>
C###    Selects only the elements at the specified direction (s1)
C###     boundary to be included in the group.
C###  Parameter:     <s2=(0/1)>
C###    Selects only the elements at the specified direction (s2)
C###     boundary to be included in the group.
C###  Parameter:     <s3=(0/1)>
C###    Selects only the elements at the specified direction (s3)
C###     boundary to be included in the group.
C###  Parameter:     <region (#s/all)[1]>
C###    Specify the element file region numbers to be defined.
C###    The all value specifies all currently defined regions.
C###    Any region numbers contained in the file but are not specified
C###    will be skipped.
C###  Description:
C###    Group elements that include the lines listed under a label. 

        OP_STRING(1)=STRING(1:IEND)//' with lines LINE#s/LINE_GROUP'
        OP_STRING(2)=BLANK(1:15)
     '    //'<as LABEL[element_'//CHAR2(IBEG2:IEND2)//']>'
        OP_STRING(3)=  BLANK(1:15)//'<(external/all)[all] '
     &    //'<s1=(0/1)>[both] <s2=(0/1)>[both] <s3=(0/1)>[both]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','GRELEM',ERROR,*9999)
      ELSE
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)

        IF(CBBREV(CO,'WITH',2,noco+1,NTCO,N3CO)) THEN
          CALL PARSE_LINES(NLLINE,NLLIST,noco-1,NRLIST,NTCO,CO,
     '      ERROR,*9999)
          
          NELIST(0)=0
          DO nl=1,NLLIST(0)
            DO ne=1,NEL(0,NLLIST(nl))
              FOUND=.FALSE.
              DO ne1=1,NELIST(0)
                IF (NEL(ne,NLLIST(nl)).EQ.NELIST(ne1)) THEN
                  FOUND=.TRUE.
                ENDIF
              ENDDO
              
              IF (.NOT.FOUND) THEN
                NELIST(0)=NELIST(0)+1
                NELIST(NELIST(0))=NEL(ne,NLLIST(nl))
              ENDIF
            ENDDO
          ENDDO

        ELSEIF(CBBREV(CO,'ALL_ELEMENTS',1,noco+1,NTCO,N3CO)) THEN
          CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '      ERROR,*9999)
        ELSE
          CALL PARSE_ELEMENTS(NEELEM,NELIST,noco-1,NRLIST,NTCO,CO,
     '      ERROR,*9999)
        ENDIF

        IF(CBBREV(CO,'UNSORTED',1,noco+1,NTCO,N3CO)) THEN
          UNSORTED=.TRUE.
        ELSE
          UNSORTED=.FALSE.
        ENDIF

        IF(.NOT.UNSORTED) THEN
C         Sort and remove duplicates from the element list
          CALL ILISTRMDUP(NELIST(0),NELIST(1),ERROR,*9999)
        ENDIF

C MPN Jul2004: select external elems on specified boundary faces
        IF(CBBREV(CO,'EXTERNAL',2,noco+1,NTCO,N3CO)) THEN
C         Select only external elements
          EXTERNAL=.TRUE.
        ELSE
          EXTERNAL=.FALSE.
        ENDIF

C       Select only elements at direction (s#) boundary
        IF(CBBREV(CO,'S1',2,noco+1,NTCO,N3CO)) THEN
          s_bdy(1)=.TRUE.
          CALL INTFROMCHAR(s_value(1),CO(N3CO+1),ERR)
          IF(ERR.NE.0) GOTO 9998
        ELSE
          s_bdy(1)=.FALSE.
        ENDIF !s1

        IF(CBBREV(CO,'S2',2,noco+1,NTCO,N3CO)) THEN
          s_bdy(2)=.TRUE.
          CALL INTFROMCHAR(s_value(2),CO(N3CO+1),ERR)
          IF(ERR.NE.0) GOTO 9998
        ELSE
          s_bdy(2)=.FALSE.
        ENDIF !s2

        IF(CBBREV(CO,'S3',2,noco+1,NTCO,N3CO)) THEN
          s_bdy(3)=.TRUE.
          CALL INTFROMCHAR(s_value(3),CO(N3CO+1),ERR)
          IF(ERR.NE.0) GOTO 9998
        ELSE
          s_bdy(3)=.FALSE.
        ENDIF !s3
        
        !Map to change direction indices from NXI to NFF format
        mxi(-1)=1
        mxi(1)=2
        mxi(-2)=3
        mxi(2)=4
        mxi(-3)=5
        mxi(3)=6

        IF(EXTERNAL) THEN
C         Select only external elements
          NELIST2(0)=0
          DO noelem=1,NELIST(0)
            ne=NELIST(noelem)
            EXTERNAL_ELEM=.FALSE.
            EXCLUDE_ELEM=.FALSE.
            DO s_dirn=1,3
              IF(NXI(s_dirn,1,ne).EQ.0.OR.NXI(-s_dirn,1,ne).EQ.0) THEN
                EXTERNAL_ELEM=.TRUE.
              ENDIF
              IF(s_bdy(s_dirn).AND.
     &          ((s_value(s_dirn).EQ.0.AND.NXI(-s_dirn,1,ne).NE.0).OR.
     &           (s_value(s_dirn).EQ.1.AND.NXI(s_dirn,1,ne).NE.0)))
     &           EXCLUDE_ELEM=.TRUE.
     
              ! test for collapsed element
              IF(NFF(mxi(s_dirn),ne).EQ.0
     &          .OR.NFF(mxi(-s_dirn),ne).EQ.0) THEN ! face is collapsed
                ! check the adjacent elements
                DO a_dirn = -3,3
                  ae=0
                  IF(a_dirn.NE.0.AND.a_dirn.NE.s_dirn
     &              .AND.a_dirn.NE.-s_dirn)
     &              ae=NXI(a_dirn,1,ne) ! adjacent element
                  IF(ae.NE.0) THEN
                    IF(s_bdy(s_dirn).AND.
     &                ((s_value(s_dirn).EQ.0.AND.NXI(-s_dirn,1,ae).NE.0)
     &                .OR.
     &                (s_value(s_dirn).EQ.1.AND.NXI(s_dirn,1,ae).NE.0)))
     &                EXCLUDE_ELEM=.TRUE.
                  ENDIF
                ENDDO
              ENDIF
            ENDDO !dirn
            IF(EXTERNAL_ELEM.AND..NOT.EXCLUDE_ELEM) THEN
              NELIST2(0)=NELIST2(0)+1
              NELIST2(NELIST2(0))=ne
            ENDIF
          ENDDO !noelem (ne)
C         Copy selected elems back into NELIST
          DO noelem=1,NELIST2(0)
            NELIST(noelem)=NELIST2(noelem)
          ENDDO !noelem
          NELIST(0)=NELIST2(0)
        ENDIF

        IF(CBBREV(CO,'AS',1,noco+1,NTCO,N3CO)) THEN
          AS=.TRUE.
          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
          CALL CUPPER(CO(N3CO+1)(IBEG:IEND),STRING2)
        ELSE
          AS=.FALSE.
        ENDIF
        
        IF(CBBREV(CO,'BY',2,noco+1,NTCO,N3CO)) THEN
          IF(CBBREV(CO,'GENERATION',3,noco+1,NTCO,N3CO)) THEN
            GROUP_TYPE=1
          ELSE IF(CBBREV(CO,'ORDER',3,noco+1,NTCO,N3CO)) THEN
            GROUP_TYPE=2
          ELSE IF(CBBREV(CO,'LOBE',3,noco+1,NTCO,N3CO)) THEN
            GROUP_TYPE=3
          ELSE IF(CBBREV(CO,'TERMINAL',3,noco+1,NTCO,N3CO)) THEN
            GROUP_TYPE=4
          ELSE
            GROUP_TYPE=0
          ENDIF
          CALL GRELEM_BY(GROUP_TYPE,NELIST,NORD,0,NXI,STRING2,
     &      AS,ERROR,*9999)
        ELSE IF(CBBREV(CO,'PARENT',3,noco+1,NTCO,N3CO)) THEN
          ne_parent=IFROMC(CO(N3CO+1))
          GROUP_TYPE=5
          CALL GRELEM_BY(GROUP_TYPE,NELIST,NORD,ne_parent,NXI,
     &      STRING2,AS,ERROR,*9999)
        ELSE IF(CBBREV(CO,'EXCLUDE',1,noco+1,NTCO,N3CO)) THEN
C AJS 02/2007: Remove specified elements from NELIST 
          NELIST3(0)=0
          CALL PARSE_ELEMENTS(NEELEM,NELIST3,N3CO,NRLIST,NTCO,
     '      CO,ERROR,*9999) !NELIST3 contains elements to be excluded
C         Ensures only included elements are in NELIST by checking 
C         CO up to the position of "exclude"
          IF(CBBREV(CO,'ELEMENTS',3,noco,N3CO,N4CO))
     '      CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,N4CO+1,CO,
     '        ERROR,*9999)
          incl=0
          DO noelem=1,NELIST(0)
            ne=NELIST(noelem)
            FOUND=.FALSE.
            noelem1=1
            DO WHILE (.NOT.FOUND.AND.noelem1.LE.NELIST3(0))
              ne1=NELIST3(noelem1)
              IF(ne1.EQ.ne) THEN 
                FOUND=.TRUE.
              ENDIF
              noelem1=noelem1+1
            ENDDO 
            IF(.NOT.FOUND) incl=ne
          ENDDO !noelem
          CALL ASSERT(incl.NE.0, 
     '      '>>Element group cannot be created. No elements specified.',
     '      ERROR,*9999) !all elements to be excluded (no elements in NELIST)
          DO noelem1=1,NELIST3(0)
            ne1=NELIST3(noelem1)
            FOUND=.FALSE.
            noelem=1
            DO WHILE (.NOT.FOUND.AND.noelem.LE.NELIST(0))
              ne=NELIST(noelem)
              IF(ne1.EQ.ne) THEN 
                NELIST(noelem)=incl !set to an included element #
                FOUND=.TRUE.
              ENDIF
              noelem=noelem+1
            ENDDO 
          ENDDO !noelem
C         Sort and remove duplicates from the element list
          !CALL ILISTRMDUP(NELIST(0),NELIST(1:NELIST(0)),ERROR,*9999)   
          CALL ILISTRMDUP(NELIST(0),NELIST(1),ERROR,*9999)   
          CALL GRELEM_SUB(NELIST,STRING2,AS,ERROR,*9999)
        ELSE
          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
          CALL GRELEM_SUB(NELIST,CO(N3CO+1)(IBEG:IEND),AS,ERROR,*9999)
        ENDIF !BY

      ENDIF !CO(noco+1)

      CALL EXITS('GRELEM')
      RETURN
 9998 ERROR=' '
 9999 CALL ERRORS('GRELEM',ERROR)
      CALL EXITS('GRELEM')
      RETURN 1
      END


