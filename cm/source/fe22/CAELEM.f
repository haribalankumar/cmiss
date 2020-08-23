      SUBROUTINE CAELEM(ISEG,ISELNO,NEELEM,NELIST,NRLIST,
     '  STRING,ERROR,*)

C#### Subroutine: CAELEM
C###  Description:
C###    CAELEM cancels element segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grou00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISELNO(NWM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NELIST(0:NEM),NRLIST(0:NRM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,IWK(6),N1ELEM,N3CO,NE_PREV,
     '  ne,noelem,noiw,nolist,no_nelist,no_nrlist,no_region,nr,NTIW
      LOGICAL ABBREV,ALL_REGIONS,CBBREV,FOUND,GROUPS,NE_EXISTS,SEGME

      CALL ENTERS('CAELEM',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM cancel elements;s
C###  Description:
C###    Cancel element segment on specified workstations. This
C###    command removes the elements from the screen but not from
C###    the problem. (GR 1/3/08 actually it does seem to remove
C###    elements from the problem)
C###  Parameter:      <on (all/WS#s)[all]>
C###    Specify the window number on which to cancel the elements.
C###    The default is to cancel the elements on all windows.
C###  Parameter:      <numbers (all/ELEMENT#s)[all]>
C###    Specify the elements which are to be cancelled. The default
C###    is to cancel all elements in the current region.
C###    A list of element numbers or the 'all' option is allowed.
C###  Parameter:      <region (all/#s)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:      <groups>
C###    Cancel all segments belonging to element groups.
C###  Parameter:    <elements (GROUP/#s/all)[all]>
C###    Specify either element groups, element numbers or all elements
C###    to be cancelled

        OP_STRING(1)=STRING(1:IEND)//';s'
        OP_STRING(2)=BLANK(1:15) //'<on (all/WS#s)[all]>'
        OP_STRING(3)=BLANK(1:15) //'<numbers (all/ELEMENT#s)[all]>'
        OP_STRING(4)=BLANK(1:15) //'<region (all/#s)[1]>'
        OP_STRING(5)=BLANK(1:15) //'<groups>'
        OP_STRING(6)=BLANK(1:15) //'<elements (GROUP/#s/all)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','CAELEM',ERROR,*9999)
      ELSE
        IF(ABBREV(COQU(noco,1),'S',1)) THEN
          SEGME=.TRUE.
          CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)
        ELSE
          SEGME=.FALSE.
        ENDIF

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)

        IF(CBBREV(CO,'GROUPS',1,noco+1,NTCO,N3CO)) THEN
          GROUPS=.TRUE.
        ELSE
          GROUPS=.FALSE.
        ENDIF

        IF(CBBREV(CO,'NUMBERS',1,noco+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),NEM,NELIST(0),NELIST(1),ERROR,*9999)
        ELSE IF(.NOT.GROUPS) THEN
          IF(CBBREV(CO,'ELEMENTS',3,noco+1,NTCO,N3CO)) THEN
            CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '        ERROR,*9999)
          ELSE
            NELIST(0)=0
            DO nolist=1,NRLIST(0)
              nr=NRLIST(nolist)
              DO noelem=NELIST(0)+1,NELIST(0)+NEELEM(0,nr)
                NELIST(noelem)=NEELEM(noelem,nr)
              ENDDO
              NELIST(0)=NELIST(0)+NEELEM(0,nr)
            ENDDO
          ENDIF
        ENDIF

        IF(SEGME) THEN
          DO noiw=1,NTIW
            iw=IWK(noiw)
            IF(IWKS(iw).GT.0) THEN
              CALL ACWK(iw,1,ERROR,*9999)
              DO nolist=1,NRLIST(0)
                nr=NRLIST(nolist)
                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  IF(ISELNO(iw,ne).GT.0) THEN
                    CALL DELETE_SEGMENT(ISELNO(iw,ne),ISEG,iw,ERROR,
     '                *9999)
                  ENDIF
                ENDDO
              ENDDO
              CALL DAWK(iw,1,ERROR,*9999)
            ENDIF
          ENDDO

        ELSE IF(GROUPS) THEN
          NTGREL=0

        ELSE IF(.NOT.SEGME) THEN
          DO no_nelist=1,NELIST(0)
            ne=NELIST(no_nelist)

!        Remove elem ne from region list & shift list down
            FOUND=.FALSE.
            DO no_nrlist=1,NRLIST(0)
              nr=NRLIST(no_nrlist)
              DO noelem=1,NEELEM(0,nr)
                IF(NEELEM(noelem,nr).eq.ne) THEN
                  N1ELEM=noelem
                  FOUND=.TRUE.
                ENDIF
              ENDDO
              IF(FOUND) THEN
                DO noelem=N1ELEM,NEELEM(0,nr)-1
                  NEELEM(noelem,nr)=NEELEM(noelem+1,nr)
                ENDDO
                NEELEM(0,nr)=NEELEM(0,nr)-1
              ENDIF !found
            ENDDO !no_nrlist
          ENDDO !no_nelist

! Reestablish maximum elem number for regions
          NEELEM(0,0)=0
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            NET(nr)=0
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              IF(ne.GT.NET(nr)) NET(nr)=ne
              NE_EXISTS=.FALSE.
              no_region=1
              DO WHILE(no_region.LE.NRT.AND..NOT.NE_EXISTS)
                n1elem=1
                DO WHILE(n1elem.LE.NEELEM(0,no_region).AND.
     '            .NOT.NE_EXISTS)
                  NE_PREV=NEELEM(n1elem,no_region)
                  IF(ne.EQ.NE_PREV) THEN
                    NE_EXISTS=.TRUE.
                  ELSE
                    n1elem=n1elem+1
                  ENDIF
                ENDDO
                no_region=no_region+1
              ENDDO
              IF(NE_EXISTS) NEELEM(0,0)=NEELEM(0,0)+1
            ENDDO !noelem
          ENDDO !no_nrlist

! Reestablish maximum element number for whole mesh
          NET(0)=0
          DO nr=1,NRT
            IF(NET(nr).GT.NET(0)) NET(0)=NET(nr)
          ENDDO

        ENDIF !segme/etc
      ENDIF

      CALL EXITS('CAELEM')
      RETURN
 9999 CALL ERRORS('CAELEM',ERROR)
      CALL EXITS('CAELEM')
      RETURN 1
      END


