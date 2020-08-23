      SUBROUTINE GRFACE(NEELEM,NELIST,NFFACE,NFLIST,NFLIST1,NPF,NFF,
     & NXI,NRLIST,STRING,ERROR,*)

C#### Subroutine: GRFACE
C###  Description:
C###    GRFACE groups faces.

C**** NTGRFA is number of face groups currently defined.
C**** LAGRFA(nogrfa) is label given to group number NOGRFA.
C**** LIGRFA(0,nogrfa) is number in list for group number NOGRFA.
C**** LIGRFA(1..,nogrfa) is list for group number NOGRFA.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grou00.cmn'
      INCLUDE 'mach00.inc'
!     Parameter List
      INTEGER NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     & NFFACE(0:NF_R_M,NRM),NFLIST(0:NFM),NFLIST1(0:NFM),
     & NPF(9,NFM),NFF(6,NEM),NXI(-NIM:NIM,0:NEIM,0:NEM),
     & NRLIST(0:NRM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER ERR,IBEG,IBEG2,IEND,IEND2,IFROMC,N1GRFA,N3CO,nedummy,
     & nf,nfdummy,nfc,nff1,nogrfa,nr,nrr,Xi_coord,s_dirn,s_value(3),
     & mxi(-3:3),nec,ne,ae
      CHARACTER CHAR*30,CHAR2*2,LABEL*30
      LOGICAL ABBREV,ALL_REGIONS,BOTHXI,CBBREV,HIGHXI,INLIST,REVERSE,
     &  UNSORTED,s_bdy(3)

      CALL ENTERS('GRFACE',*9999)
      
      s_dirn=0
      
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
C        CHAR2=CFROMI(NTGRLI+1,'(I2)')
        WRITE(CHAR2,'(I2)') NTGRLI+1
        CALL STRING_TRIM(CHAR2,IBEG2,IEND2)

C---------------------------------------------------------------------

C#### Command: FEM group faces face#s/FACE_GROUP/allfaces
C###  Parameter:     <as LABEL[face_1]>
C###    Specifies the name of the face group.
C###  Parameter:     <(sorted/unsorted)[sorted]>
C###    The sorted command sorts the face numbers defined
C###    in the group in ascending order and removes duplicates.
C###  Parameter:     <(internal/external)[all]>
C###    Select only internal or external faces (default is all faces). 
C###  Parameter:     <xi(1/2/3) (high/low)>[all]
C###    Identifies the faces with the specified xi coordinate  
C###    as a normal. xi_i=0 (low) or xi_i=1 (high) can be specified
C###    separately. Appending 'not' selects all faces except those
C###    indicated.
C###  Parameter:     <reverse_selection>
C###    Selects all faces except those indicated.
C###  Parameter:     <element (all/GROUP/#s)[all] <s(1/2/3)=(0/1)>[all]>
C###    Select the face numbers that only belong to the specified elements.
C###    The s(1/2/3)=(0/1) option further restricts the selection to include
C###    only the faces on the specified boundary of the element group.
C###    This is different to the xi(1/2/3) option which doesn't restrict to
C###    faces on a boundary.
C###  Parameter:     <region (#s/all)[1]>
C###    Specify the face file region numbers to be defined.
C###    The all value specifies all currently defined regions.
C###    Any region numbers contained in the file but are not specified
C###    will be skipped.
C###  Description:
C###    Group faces under a label. By default the list is sorted
C###    by face numbers into ascending order and duplicates are
C###    removed - the UNSORTED option prevents this.

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'FACE#s/FACE_GROUP/allfaces'
        OP_STRING(3)=BLANK(1:15)
     '    //'<as LABEL[face_'//CHAR2(IBEG2:IEND2)//']>'
        OP_STRING(4)=BLANK(1:15)//'<(sorted/unsorted)[sorted]>'
        OP_STRING(5)=BLANK(1:15)//'<(internal/external)[all]>'
        OP_STRING(6)=BLANK(1:15)//'<xi(1/2/3) (high/low)>[all]'
        OP_STRING(7)=BLANK(1:15)//'<reverse_selection>'
        OP_STRING(8)=BLANK(1:15)//'<element (all/GROUP/#s)[all]>'
        OP_STRING(9)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe27','doc','GRFACE',ERROR,*9999)
      ELSE
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)

C MPN 28Apr2003: added several options for selecting faces
        IF(CBBREV(CO,'ALLFACES',2,noco+1,noco+2,N3CO)) THEN
C         Select all faces in specified regions
          nff1=0
          DO nrr=1,NRLIST(0)
            nr=NRLIST(nrr)
            DO nfc=1,NFFACE(0,nr)
              nf=NFFACE(nfc,nr)
              nff1=nff1+1
              IF(nff1.LE.NFM) NFLIST(nff1)=nf
            ENDDO !nf
          ENDDO !nrr (nr)
          WRITE(CHAR2,'(I2)') nff1
          CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
          CALL ASSERT(nff1.LE.NFM,'>>Increase NFM to be at least '
     &      //CHAR2(IBEG2:IEND2),ERROR,*9999)
          NFLIST(0)=nff1
        ELSE
C         Parse list of specified faces
          CALL PARSE_FACES(NFFACE,NFLIST,noco-1,NRLIST,NTCO,CO,
     &      ERROR,*9999)
        ENDIF

        IF(CBBREV(CO,'UNSORTED',2,noco+1,NTCO,N3CO)) THEN
          UNSORTED=.TRUE.
        ELSE
          UNSORTED=.FALSE.
        ENDIF

        IF(CBBREV(CO,'REVERSE_SELECTION',2,noco+1,NTCO,N3CO)) THEN
          REVERSE=.TRUE.
        ELSE
          REVERSE=.FALSE.
        ENDIF
        
        nff1=0

C new MPN 21Jan2005: adding 'elements' options for selecting faces only
C                    from certain elements

        IF(CBBREV(CO,'ELEMENTS',2,noco+1,NTCO,N3CO)) THEN
          CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     &      ERROR,*9999)
     
C GDR 9Apr2005 Select only elements at direction (s(1/2/3)) boundary.
C              A lot of this is borrowed from GRELEM.f
          IF(CBBREV(CO,'S1',2,noco+1,NTCO,N3CO)) THEN
            s_bdy(1)=.TRUE.
            CALL INTFROMCHAR(s_value(1),CO(N3CO+1),ERR)
            IF(ERR.NE.0) GOTO 9999
          ELSE
            s_bdy(1)=.FALSE.
          ENDIF !s1
  
          IF(CBBREV(CO,'S2',2,noco+1,NTCO,N3CO)) THEN
            s_bdy(2)=.TRUE.
            CALL INTFROMCHAR(s_value(2),CO(N3CO+1),ERR)
            IF(ERR.NE.0) GOTO 9999
          ELSE
            s_bdy(2)=.FALSE.
          ENDIF !s2
  
          IF(CBBREV(CO,'S3',2,noco+1,NTCO,N3CO)) THEN
            s_bdy(3)=.TRUE.
            CALL INTFROMCHAR(s_value(3),CO(N3CO+1),ERR)
            IF(ERR.NE.0) GOTO 9999
          ELSE
            s_bdy(3)=.FALSE.
          ENDIF !s3
          
          IF(s_bdy(1).OR.s_bdy(2).OR.s_bdy(3)) THEN
            !Map to change direction indices from NXI to NFF format
            mxi(-1)=1
            mxi(1)=2
            mxi(-2)=3
            mxi(2)=4
            mxi(-3)=5
            mxi(3)=6
            DO nec=1,NELIST(0)
              ne=NELIST(nec)
              nf=0
                                          
              DO s_dirn=1,3
                IF(s_bdy(s_dirn)) THEN
                  IF(s_value(s_dirn).EQ.0) THEN
                    ae=NXI(-s_dirn,1,ne) ! adjacent element
                    ! If there is no adjacent element or the adjacent
                    ! element is not in the element group then we are
                    ! on a boundary so include this face.
                    IF(ae.EQ.0.OR.
     &                .NOT.INLIST(ae,NELIST(1),NELIST(0),nedummy)) THEN
                      nf=NFF(mxi(-s_dirn),ne)
                    ENDIF
                  ELSE IF(s_value(s_dirn).EQ.1) THEN
                    ae=NXI(s_dirn,1,ne) ! adjacent element
                    IF(ae.EQ.0.OR.
     &                .NOT.INLIST(ae,NELIST(1),NELIST(0),nedummy)) THEN
                      nf=NFF(mxi(s_dirn),ne)
                    ENDIF
                  ENDIF
                ENDIF
              ENDDO             
              
              IF(nf.NE.0) THEN
                nff1=nff1+1
                NFLIST1(nff1)=nf
              ENDIF                
            ENDDO
            NFLIST1(0)=nff1
            CALL ILIST_COPY(NFLIST1(0)+1,NFLIST1(0),NFLIST(0))  
          ELSE
     
            !nff1=0
            DO nfc=1,NFLIST(0)
              nf=NFLIST(nfc) !global face number
              IF(INLIST(NPF(6,nf),NELIST(1),NELIST(0),nedummy).OR.
     &        INLIST(NPF(7,nf),NELIST(1),NELIST(0),nedummy)) THEN
C             Face number is contained within one of the specified elements
                nff1=nff1+1
                NFLIST1(nff1)=nf
              ENDIF
            ENDDO !nff (nf)
            NFLIST1(0)=nff1
            CALL ILIST_COPY(NFLIST1(0)+1,NFLIST1(0),NFLIST(0))
          ENDIF
        ENDIF
C end new MPN 21Jan2005

        IF(CBBREV(CO,'INTERNAL',2,noco+1,NTCO,N3CO)) THEN
C         Select only internal Xi faces from the list
          nff1=0
          DO nfc=1,NFLIST(0)
            nf=NFLIST(nfc) !global face number
            IF(NPF(5,nf).EQ.2) THEN !nf is an internal face
              nff1=nff1+1
              NFLIST1(nff1)=nf
            ENDIF
          ENDDO !nff (nf)
          NFLIST1(0)=nff1
          CALL ILIST_COPY(NFLIST1(0)+1,NFLIST1(0),NFLIST(0))
        ELSE IF(CBBREV(CO,'EXTERNAL',2,noco+1,NTCO,N3CO)) THEN
C         Select only external Xi faces from the list
          nff1=0
          DO nfc=1,NFLIST(0)
            nf=NFLIST(nfc) !global face number
            IF(NPF(5,nf).EQ.1) THEN !nf is an external face
              nff1=nff1+1
              NFLIST1(nff1)=nf
            ENDIF
          ENDDO !nff (nf)
          NFLIST1(0)=nff1
          CALL ILIST_COPY(NFLIST1(0)+1,NFLIST1(0),NFLIST(0))
        ENDIF

        IF(CBBREV(CO,'XI1',3,noco+1,NTCO,N3CO)
     &    .OR.CBBREV(CO,'XI2',3,noco+1,NTCO,N3CO)
     &    .OR.CBBREV(CO,'XI3',3,noco+1,NTCO,N3CO)) THEN
          Xi_coord=IFROMC(CO(N3CO)(3:3))
          N3CO=N3CO+1
          IF(ABBREV(CO(N3CO),'HIGH',2)) THEN
            BOTHXI=.FALSE.
            HIGHXI=.TRUE.
          ELSE IF(ABBREV(CO(N3CO),'LOW',2)) THEN
            BOTHXI=.FALSE.
            HIGHXI=.FALSE.
          ELSE
            BOTHXI=.TRUE.
            HIGHXI=.FALSE.
          ENDIF
C         Add only specified Xi faces to the list
C!!! 9Apr2005 GDR This is broken because it only considers the first element
C             connected to the face. For it to work properly for all cases 
C             NPF(9,nf) needs to also be considered... but that isn't
C             straightforward.
          nff1=0
          DO nfc=1,NFLIST(0)
            nf=NFLIST(nfc) !global face number
            IF(Xi_coord.EQ.1.AND.NPF(8,nf).EQ.1.AND.
     &        (.NOT.HIGHXI.OR.BOTHXI).OR.
     &        Xi_coord.EQ.1.AND.NPF(8,nf).EQ.2.AND.
     &        (HIGHXI.OR.BOTHXI).OR.
     &        Xi_coord.EQ.2.AND.NPF(8,nf).EQ.3.AND.
     &        (.NOT.HIGHXI.OR.BOTHXI).OR.
     &        Xi_coord.EQ.2.AND.NPF(8,nf).EQ.4.AND.
     &        (HIGHXI.OR.BOTHXI).OR.
     &        Xi_coord.EQ.3.AND.NPF(8,nf).EQ.5.AND.
     &        (.NOT.HIGHXI.OR.BOTHXI).OR.
     &        Xi_coord.EQ.3.AND.NPF(8,nf).EQ.6.AND.
     &        (HIGHXI.OR.BOTHXI)) THEN !NPF(8,nf) is the local face number
              nff1=nff1+1
              NFLIST1(nff1)=nf
            ENDIF
          ENDDO !nff (nf)
          NFLIST1(0)=nff1
          CALL ILIST_COPY(NFLIST1(0)+1,NFLIST1(0),NFLIST(0))
        ENDIF

        IF(REVERSE) THEN
C         Select all faces except those indicated in the list
          nff1=0
          DO nrr=1,NRLIST(0)
            nr=NRLIST(nrr)
            DO nfc=1,NFFACE(0,nr)
              nf=NFFACE(nfc,nr)
              IF(.NOT.INLIST(nf,NFLIST(1),NFLIST(0),nfdummy)) THEN
                nff1=nff1+1
                NFLIST1(nff1)=nf
              ENDIF
            ENDDO !nf
          ENDDO !nrr (nr)
          NFLIST1(0)=nff1
          CALL ILIST_COPY(NFLIST1(0)+1,NFLIST1(0),NFLIST(0))
        ENDIF

        IF(.NOT.UNSORTED) THEN
C         Sort and remove duplicates from the face list
          CALL ILISTRMDUP(NFLIST(0),NFLIST(1),ERROR,*9999)
        ENDIF

C        CALL ASSERT(NFLIST(0).LE.GRFA_MAXDATA,
C     '    '>>Increase array sizes in grou00.cmn',ERROR,*9999)

        IF(CBBREV(CO,'AS',1,noco+1,NTCO,N3CO)) THEN
          !Check whether group name already exists
          CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
C          CHAR=CUPPER(CO(N3CO+1)(IBEG:IEND))
          CALL CUPPER(CO(N3CO+1)(IBEG:IEND),CHAR)
          N1GRFA=0
          DO nogrfa=1,NTGRFA
C            LABEL=CUPPER(LAGRFA(nogrfa))
            CALL CUPPER(LAGRFA(nogrfa),LABEL)
            CALL STRING_TRIM(LABEL,IBEG2,IEND2)
            IF(CHAR(IBEG:IEND).EQ.LABEL(IBEG2:IEND2)) THEN
              N1GRFA=nogrfa !is existing group label ID
              GO TO 100
            ENDIF
          ENDDO
 100      IF(N1GRFA.EQ.0) THEN !need new group label
            NTGRFA=NTGRFA+1 !is new total #groups
            CALL ASSERT(NTGRFA.LE.GRFA_MAXGRP,
     &        '>>Increase array sizes in grou00.cmn',ERROR,*9999)
            N1GRFA=NTGRFA
            CALL STRING_TRIM(CO(N3CO+1),IBEG,IEND)
            LAGRFA(N1GRFA)=CO(N3CO+1)(IBEG:IEND) !new face group label
          ENDIF

        ELSE
          NTGRFA=NTGRFA+1 !is new total
          CALL ASSERT(NTGRFA.LE.GRFA_MAXGRP,
     &      '>>Increase array sizes in grou00.cmn',ERROR,*9999)
          N1GRFA=NTGRFA
C          CHAR2=CFROMI(N1GRFA,'(I2)')
          WRITE(CHAR2,'(I2)') N1GRFA
          CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
          LAGRFA(N1GRFA)='face_'//CHAR2(IBEG2:IEND2) !new label
        ENDIF

C TVK 09July1999: dynamic groups (entire routine is new)
        NLIGRFA(N1GRFA)=NFLIST(0)
        CALL ALLOCATE_MEMORY(NFLIST(0),0,INTTYPE,LIGRFA_PTR(N1GRFA),
     &    MEM_INIT,ERROR,*9999)
        CALL ILIST_COPY(NFLIST(0),NFLIST(1),%VAL(LIGRFA_PTR(N1GRFA)))
      ENDIF

      CALL EXITS('GRFACE')
      RETURN
 9999 CALL ERRORS('GRFACE',ERROR)
      CALL EXITS('GRFACE')
      RETURN 1
      END


