      SUBROUTINE CANODE(ISEG,ISNONO,NPNODE,NPLIST,NRLIST,
     '  STRING,ERROR,*)

C#### Subroutine: CANODE
C###  Description:
C###    CANODE cancels node segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grou00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISNONO(NWM,NPM),NPNODE(0:NP_R_M,0:NRM),
     '  NPLIST(0:NPM),NRLIST(0:NRM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,
     '  INSTAT,IPICK,ISEGM,iw,IWK(6),N1NODE,N3CO,noiw,nonode,
     '  no_nplist,no_nrlist,no_region,np,NP_PREV,nr,NTIW
      LOGICAL ABBREV,ALL_REGIONS,CBBREV,FOUND,GROUPS,
     '  MOUSE,NP_EXISTS,SEGME

      CALL ENTERS('CANODE',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM cancel nodes;s
C###  Description:
C###    Cancel node segment on specified workstations.
C###  Parameter:      <on (all/WS#s)[all]>
C###    Specify the window number on which to cancel the nodes.
C###    The default is to cancel the nodes on all windows.
C###  Parameter:      <region (all/#s)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:      <nodes (all/GROUP/#s)[all]>
C###    Specify either node groups, node numbers or all nodes
C###    to be cancelled

        OP_STRING(1)=STRING(1:IEND)//';s'
        OP_STRING(2)=BLANK(1:15) //'<on (all/WS#s)[all]>'
        OP_STRING(3)=BLANK(1:15) //'<region (all/#s)[1]>'
        OP_STRING(4)=BLANK(1:15) //'<nodes (all/GROUP/#s)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM cancel nodes;m
C###  Parameter:      <on WS#[1]>
C###    Specify the window number on which to cancel the nodes.
C###    The default is to cancel the nodes on window 1.
C###  Parameter:      <numbers (all/NODE#s)[all]>
C###    Specify the nodes which are to be cancelled. The default
C###    is to cancel all nodes in the current region.
C###    A list of node numbers or the 'all' option is allowed.
C###  Parameter:      <region (all/#s)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Description:
C###    Cancel node segment with mouse on specified workstation.

        OP_STRING(1)=STRING(1:IEND)//';m'
        OP_STRING(2)=BLANK(1:15) //'<on WS#[1]>'
        OP_STRING(3)=BLANK(1:15) //'<numbers (all/NODE#s)[all]>'
        OP_STRING(4)=BLANK(1:15) //'<region (all/#s)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','CANODE',ERROR,*9999)
      ELSE
        CALL CHECKQ(' SM',noco,1,CO,COQU,STRING,*1)

        SEGME=.FALSE.
        MOUSE=.FALSE.
        IF(ABBREV(COQU(noco,1),'S',1)) THEN
          SEGME=.TRUE.
        ELSE IF(ABBREV(COQU(noco,1),'M',1)) THEN
          MOUSE=.TRUE.
        ENDIF
        IF(SEGME.OR.MOUSE) THEN
          CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)
        ENDIF

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)

        IF(CBBREV(CO,'GROUPS',1,noco+1,NTCO,N3CO)) THEN
          GROUPS=.TRUE.
        ELSE
          GROUPS=.FALSE.
        ENDIF

        IF(CBBREV(CO,'NUMBERS',2,noco+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),NP_R_M,NPLIST(0),NPLIST(1),ERROR,*9999)
        ELSE IF(.NOT.GROUPS) THEN
          IF(CBBREV(CO,'NODES',2,noco+1,NTCO,N3CO)) THEN
            CALL PARSE_NODES(NPNODE,NPLIST,noco,NRLIST,NTCO,CO,
     '        ERROR,*9999)
          ELSE !cancel all nodes
            NPLIST(0)=0
            DO no_nrlist=1,NRLIST(0)
              nr=NRLIST(no_nrlist)
              DO nonode=NPLIST(0)+1,NPLIST(0)+NPNODE(0,nr)
                NPLIST(nonode)=NPNODE(nonode,nr)
              ENDDO
              NPLIST(0)=NPLIST(0)+NPNODE(0,nr)
            ENDDO
          ENDIF
        ENDIF

        IF(SEGME) THEN
          DO noiw=1,NTIW
            iw=IWK(noiw)
            IF(IWKS(iw).GT.0) THEN
              CALL ACWK(iw,1,ERROR,*9999)
              DO no_nrlist=1,NRLIST(0)
                nr=NRLIST(no_nrlist)
                DO nonode=1,NPNODE(0,nr)
                  np=NPNODE(nonode,nr)
                  IF(ISNONO(iw,np).GT.0) THEN
                    CALL DELETE_SEGMENT(ISNONO(iw,np),ISEG,iw,
     '                ERROR,*9999)
                  ENDIF
                ENDDO
              ENDDO
              CALL DAWK(iw,1,ERROR,*9999)
            ENDIF
          ENDDO

        ELSE IF(MOUSE) THEN
          iw=IWK(1)
          CALL ACWK(iw,0,ERROR,*9999)
          WRITE(OP_STRING,'('' >>Pick nodes on '',I1)') iw
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          CALL PICK(iw,'EVENT',INSTAT,ISEGM,IPICK,ERROR,*9999)
C          CONTINUE=.TRUE.
C          DO WHILE(CONTINUE)
c           CALL EVENT(ID_WS,ID_DEVICE,INPUT_STATUS,CLASS,IDATA,
c    '        R4DATA,SDATA,ERROR,*9999)
c            IF(DOP) THEN
c              WRITE(OP_STRING,*) ' INPUT_CLASS=',CLASS
c              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
c            ENDIF
c            IF(CLASS(1:4).EQ.'PICK') THEN
c              ISEGM=IDATA(1)
c            ELSE
c             CALL INPUT_MODE(iw,LD1,'PICK','REQUEST',ERROR,*9999)
c              CONTINUE=.FALSE.
c            ENDIF
C          ENDDO
          CALL DAWK(iw,0,ERROR,*9999)

        ELSE IF(GROUPS) THEN
          NTGRNO=0

        ELSE
          DO no_nplist=1,NPLIST(0)
            np=NPLIST(no_nplist)
!Removing code since GKS and PHIGS gone  AJP 5/5/95
c            IF(GKS.OR.PHIGS) THEN
c              DO noiw=1,2*NJT-3+IMAP
c                iw=IWK(noiw)
c                IF(IWKS(iw).GT.0) THEN
c                  CALL ACWK(iw,1,ERROR,*9999)
c                  IF(ISNONO(iw,np).GT.0) THEN
c                    CALL DELETE_SEGMENT(ISNONO(iw,np),ISEG,iw,
c     '                ERROR,*9999)
c                  ENDIF
c                  CALL DAWK(iw,1,ERROR,*9999)
c                ENDIF
c              ENDDO
c            ENDIF
!end

!        Remove node np from region list & shift list down
            FOUND=.FALSE.
            DO no_nrlist=1,NRLIST(0)
              nr=NRLIST(no_nrlist)
              DO nonode=1,NPNODE(0,nr)
                IF(NPNODE(nonode,nr).EQ.np) THEN
                  N1NODE=nonode
                  FOUND=.TRUE.
                ENDIF
              ENDDO
              IF(FOUND) THEN
                DO nonode=N1NODE,NPNODE(0,nr)-1
                  NPNODE(nonode,nr)=NPNODE(nonode+1,nr)
                ENDDO
                NPNODE(0,nr)=NPNODE(0,nr)-1
              ENDIF !found
            ENDDO !no_nrlist
          ENDDO !no_nplist

! Reestablish maximum node number for regions
          NPNODE(0,0)=0
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            NPT(nr)=0
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
              IF(np.GT.NPT(nr)) NPT(nr)=np
              NP_EXISTS=.FALSE.
              no_region=1
              DO WHILE(no_region.LE.NRT.AND..NOT.NP_EXISTS)
                n1node=1
                DO WHILE(n1node.LE.NPNODE(0,no_region).AND.
     '            .NOT.NP_EXISTS)
                  NP_PREV=NPNODE(n1node,no_region)
                  IF(np.EQ.NP_PREV) THEN
                    NP_EXISTS=.TRUE.
                  ELSE
                    n1node=n1node+1
                  ENDIF
                ENDDO
                no_region=no_region+1
              ENDDO
              IF(NP_EXISTS) NPNODE(0,0)=NPNODE(0,0)+1
            ENDDO !nonode
          ENDDO !no_nrlist

! Reestablish maximum node number for whole mesh
          NPT(0)=0
          DO nr=1,NRT
            IF(NPT(nr).GT.NPT(0)) NPT(0)=NPT(nr)
          ENDDO

        ENDIF !segme/mouse/etc
      ENDIF

      CALL EXITS('CANODE')
      RETURN
 9999 CALL ERRORS('CANODE',ERROR)
      CALL EXITS('CANODE')
      RETURN 1
      END


