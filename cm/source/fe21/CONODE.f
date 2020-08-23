      SUBROUTINE CONODE(NHP,NKH,NKJ,NPLIST,NPNODE,NVJP,
     '  XP,STRING,ERROR,*)

C#### Subroutine: CONODE
C###  Description:
C###    CONODE copies nodes.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NHP(NPM,0:NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKJ(NJM,NPM),NPLIST(0:NPM),NPNODE(0:NP_R_M,0:NRM),NVJP(NJM,NPM)
      REAL*8 XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IBEG4,IEND,IEND4,N3CO,nc,nh,nhx,nj,nk,nolist,nonode,
     '  NPOLD,NPNEW,nr,NSCALE,nx
      REAL*8 RSCALE(3)
      CHARACTER CHAR4*4
      LOGICAL CBBREV

      CALL ENTERS('CONODE',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        WRITE(CHAR4(1:4),'(I4)') NPT(1)
        CALL STRING_TRIM(CHAR4,IBEG4,IEND4)

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)
     '    //'<from NODE#s[1..'//CHAR4(IBEG4:IEND4)//']>'
        OP_STRING(3)=BLANK(1:15)
     '    //'<scale COORD_SCALE_FACTOR#s[1,1,1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
C---------------------------------------------------------------------
C#### Command: FEM copy nodes
C###  Description:
C###    Copies nodes to new node numbers specified
C###  Parameter: <from NODE#s[1..0]>
C###  Parameter: <scale COORD_SCALE_FACTOR#s[1,1,1]>
C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','CONODE',ERROR,*9999)
      ELSE
        nr=1 !Needs fixing
        nc=1 !ditto
        nx=1 !temporary

        IF(CBBREV(CO,'FROM',1,noco+1,NTCO,N3CO)) THEN
          CALL PARSIL(CO(N3CO+1),NPT(nr),NPLIST(0),NPLIST(1),ERROR,
     '      *9999)
        ELSE
          NPLIST(0)=0
          DO nr=1,NRT
            DO nonode=NPLIST(0)+1,NPLIST(0)+NPNODE(0,nr)
              NPLIST(nonode)=NPNODE(nonode,nr)
            ENDDO
            NPLIST(0)=NPLIST(0)+NPNODE(0,nr)
          ENDDO
        ENDIF
        DO nj=1,NJT
          RSCALE(nj)=1.d0
        ENDDO
        IF(CBBREV(CO,'SCALE',1,noco+1,NTCO,N3CO)) THEN
          CALL PARSRL(CO(N3CO+1),3,NSCALE,RSCALE,ERROR,*9999)
        ENDIF
        CALL ASSERT(NPT(1)+NPLIST(0).LE.NP_R_M,'>>NP_R_M too small',
     '    ERROR,*9999)

        DO nolist=1,NPLIST(0)
          NPOLD=NPLIST(nolist)
          NPNEW=NPT(1)+nolist
C GMH 8/1/97 Update cmgui link
          CALL NODE_CHANGE(NPNEW,.TRUE.,ERROR,*9999)
          NPNODE(NPNODE(0,nr)+nolist,1)=NPNEW
          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
            NKJ(nj,NPNEW)=NKJ(nj,NPOLD)
            NVJP(nj,NPNEW)=NVJP(nj,NPOLD)
            XP(1,1,nj,NPNEW)=RSCALE(nj)*XP(1,1,nj,NPOLD)
            DO nk=2,NKJ(nj,NPNEW)
              XP(nk,1,nj,NPNEW)=XP(nk,1,nj,NPOLD)
            ENDDO
          ENDDO
          NHP(NPNEW,nr,nx)=NHP(NPOLD,nr,nx)
          DO nhx=1,NHP(NPNEW,nr,nx)
            nh=NH_LOC(nhx,nx)
            NKH(nh,NPNEW,nc,nr)=NKH(nh,NPOLD,nc,nr)
          ENDDO
        ENDDO
        NPT(0)=NPT(0)+NPLIST(0)
        NPT(1)=NPT(1)+NPLIST(0)
        NPNODE(0,0)=NPNODE(0,0)+NPLIST(0)
        NPNODE(0,1)=NPNODE(0,1)+NPLIST(0)
      ENDIF

      CALL EXITS('CONODE')
      RETURN
 9999 CALL ERRORS('CONODE',ERROR)
      CALL EXITS('CONODE')
      RETURN 1
      END


