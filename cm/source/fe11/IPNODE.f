      SUBROUTINE IPNODE(TYPE,nc,NHP,NKH,NKJ,N_OFFSET,NP_INTERFACE,
     &  NPLIST,NPNODE,nr,NRLIST,NVHP,NVJP,XP,ZP,ERROR,*)

C#### Subroutine: IPNODE
C###  Description:
C###    IPNODE inputs global coordinates.

C#### Variable: nk
C###  Type: INTEGER
C###  Set_up: IPNODE
C###  Description:
C###    <HTML>
C###    nk is the nodal derivative loop variable.
C###    <PRE>
C###    nk=1 function
C###    nk=2 derivative wrt direction 1
C###    nk=3 derivative wrt direction 2
C###    nk=4 derivative wrt directions 1,2
C###    nk=5 derivative wrt direction 3
C###    nk=6 derivative wrt directions 1,3
C###    nk=7 derivative wrt directions 2,3
C###    nk=8 derivative wrt directions 1,2,3
C###    </PRE>
C###    </HTML>

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'defn00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'parameters.inc'
      
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi10.cmn'

!     Parameter List
      INTEGER nc,
     '  NHP(NPM),NKH(NHM,NPM,NCM),NKJ(NJM,NPM),N_OFFSET,
     '  NP_INTERFACE(0:NPM,0:3),NPLIST(0:NPM),
     '  NPNODE(0:NP_R_M,0:NRM),nr,NVHP(NHM,NPM,NCM),NVJP(NJM,NPM),
     &  NRLIST(0:NRM)
      REAL*8 XP(NKM,NVM,NJM,NPM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER TYPE*(*),ERROR*(*)
!     Local Variables
      INTEGER ICHAR,INFO,N,n1node,N2NODE,nj,njj,
     '  NJTOT,nk,NKJP(3),nonode,NOQUES,no_region,
     '  np,NP2,NP_MAX,NP_PREV,NTNODE,nu,NUK(8),nv,NV_MAX,
     &  N3CO
      CHARACTER CHAR1*6,CHAR2*1,CHAR3*3,CHAR4*12
      LOGICAL EXISTS,FILEIP,NP_EXISTS,PROMPT_NV(3),
     &  CBBREV

      DATA NUK/1,2,4,6,7,9,10,11/

      CALL ENTERS('IPNODE',*9999)

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      IF(ITYP10(nr).EQ.4) THEN
        FORMAT='($,'' Specify the focus position [1.0]: '',D12.4)'
        IF(IOTYPE.EQ.3) RDATA(1)=FOCUS
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RONE,-RMAX,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) FOCUS=RDATA(1)
      ENDIF

C KAT 2004-01-27: nb unused.
C       IF(NET(nr).GT.0) THEN
C         ELEM=.TRUE.
C       ELSE
C         ELEM=.FALSE.
C       ENDIF

C GMH 14/2/97 Destroy this region and recreate -
C             necessary for cmgui locking
      IF(NPNODE(0,nr).GT.0) THEN
        CALL REGION_DESTROY(nr,ERROR,*9999)
      ENDIF
      DO np=0,NPM
        NPLIST(np)=0
      ENDDO !np
C GMH 22/12/96 Find all nodes

      IF((IOTYPE.EQ.3).AND.
     &   (CBBREV(CO,'NODES',2,noco+1,NTCO,N3CO))) THEN
          CALL PARSE_NODES(NPNODE,NPLIST,noco,NRLIST,NTCO,CO,
     &    ERROR,*9999)
      ELSE
      DO nonode=1,NPNODE(0,nr)
        np=NPNODE(nonode,nr)
C!!! n^2 operation seems a big waste of effort just to get a default
        DO N=1,NPLIST(0)
          IF(NPLIST(N).EQ.np) GOTO 10
        ENDDO
        NPLIST(0)=NPLIST(0)+1
        NPLIST(NPLIST(0))=np
C       CMGUI link - destroy the node
        IF(.NOT.ADD.AND.IOTYPE.NE.3) THEN
          CALL NODE_DESTROY(np,ERROR,*9999)
          DO nj=1,NJT
            NVJP(nj,np)=0
          ENDDO
        ENDIF
 10     CONTINUE
      ENDDO !np
      ENDIF
C KAT 2004-01-21: This could be but wasn't used for writing all the
C nodes associated with the region even if they weren't read into this
C region (into NPNODE).
C       DO noelem=1,NEELEM(0,nr)
C         ne=NEELEM(noelem,nr)
C         nb=NBJ(1,ne)
C         DO nn=1,NNT(nb)
C           np=NPNE(nn,nb,ne)
C           DO N=1,NPLIST(0)
C             IF(NPLIST(N).EQ.NP) GOTO 11
C           ENDDO
C           NPLIST(0)=NPLIST(0)+1
C           NPLIST(NPLIST(0))=NP
C C         CMGUI link - destroy the node
C           CALL NODE_DESTROY(np,ERROR,*9999)
C  11       CONTINUE
C         ENDDO
C       ENDDO
      CALL ISORT(NPLIST(0),NPLIST(1))

      IDEFLT(1)=NPLIST(0)
      IF(IDEFLT(1).EQ.0) IDEFLT(1)=1
      WRITE(CHAR1,'(I6)') IDEFLT(1)
      IF(IOTYPE.EQ.3) THEN
        !NTNODE=NPNODE(0,nr)
        NTNODE=NPLIST(0)
        IDATA(1)=NTNODE
      ENDIF
      FORMAT='($,'' The number of nodes is ['//CHAR1(1:6)//']: '',I6)'
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,NPM,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) THEN
C RGB 17/4/98 incorporating 'ADD' to see if the old nodes
C has a greater total in region nr.
        NTNODE=IDATA(1)
      ENDIF
      CALL ASSERT(NTNODE.LE.NP_R_M,'>>Increase NP_R_M',ERROR,*9999)

      IDEFLT(1)=NJT
      WRITE(CHAR1,'(I1)') IDEFLT(1)
      FORMAT='($,'' Number of coordinates ['//CHAR1(1:1)//']: '',I1)'
      IF(IOTYPE.EQ.3) IDATA(1)=NJT
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) THEN
        CALL ASSERT(IDATA(1).EQ.NJT,'>>Define coordinates first',
     '    ERROR,*9999)
      ENDIF

      IF(NTNODE.GT.0) THEN
!news MPN 16-Nov-94
C       ask for version prompting
        DO nj=1,NJT
          WRITE(CHAR1,'(I1)') nj
          FORMAT='($,'' Do you want prompting for different versions '
     '      //'of nj='//CHAR1(1:1)//' [N]? '',A)'
          IF(IOTYPE.EQ.3) THEN
            PROMPT_NV(nj)=.FALSE.
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
              IF(NVJP(nj,np).GT.1) PROMPT_NV(nj)=.TRUE.
            ENDDO
            IF(PROMPT_NV(nj)) THEN
              ADATA(1)='Y'
            ELSE
              ADATA(1)='N'
            ENDIF
          ENDIF
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            IF(ADATA(1).EQ.'Y') THEN
              PROMPT_NV(nj)=.TRUE.
            ELSE
              PROMPT_NV(nj)=.FALSE.
            ENDIF
          ENDIF
        ENDDO
!newe

        NOQUES=0
        DO nj=1,NJT
          WRITE(CHAR1,'(I1)') nj
          IF(IOTYPE.EQ.3) THEN
            IF(TYPE(1:9).EQ.'REFERENCE') THEN
              !Find the maximum number of derivatives in the node list
              NKJP(nj)=0
              !DO nonode=1,NTNODE
              DO nonode=1,NPNODE(0,nr)
                IF(NKJ(nj,NPNODE(nonode,nr)).GT.NKJP(nj)) THEN
                  NKJP(nj)=NKJ(nj,NPNODE(nonode,nr))
                ENDIF
              ENDDO
            ELSE IF(TYPE(1:9).EQ.'DEPENDENT') THEN
              NKJP(nj)=0
              !DO nonode=1,NTNODE
              DO nonode=1,NPNODE(0,nr)
                IF(NKH(nj,NPNODE(nonode,nr),1).GT.NKJP(nj)) THEN
                  NKJP(nj)=NKH(nj,NPNODE(nonode,nr),nc)
                ENDIF
              ENDDO
            ENDIF
            IDATA(1)=NKJP(nj)-1
          ENDIF
          FORMAT='($,'' The number of derivatives for coordinate '
     '      //CHAR1(1:1)//' is [0]: '',I1)'
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,NKM-1,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NKJP(nj)=IDATA(1)+1
        ENDDO
        N2NODE=0
        NOQUES=0
        IF(IOTYPE.NE.3.AND..NOT.ADD) THEN
          NPNODE(0,nr)=0
        ENDIF
        DO nonode=1,NTNODE
          IF(NPLIST(nonode).GT.0) THEN
            IDEFLT(1)=NPLIST(nonode)
          ELSE IF(NPLIST(nonode).EQ.0) THEN
            IF(N2NODE.EQ.0) THEN
              IF(ADD) THEN
                IDEFLT(1)=NPT(1)+1
              ELSE
                NP_MAX=0
                DO no_region=1,NRT
                  IF(no_region.NE.NR) THEN
                    DO n1node=1,NPNODE(0,no_region)
                      np=NPNODE(n1node,no_region)
                      IF(np.GT.NP_MAX) NP_MAX=NP
                    ENDDO
                  ENDIF
                ENDDO
                IDEFLT(1)=NP_MAX+1
              ENDIF
            ELSE ! IF(nonode.GT.1) THEN
              IDEFLT(1)=NPNODE(N2NODE,nr)+1
            ENDIF
          ENDIF
          WRITE(CHAR1,'(I6)') IDEFLT(1)
          IF(IOTYPE.EQ.3) THEN
            !np=NPNODE(nonode,nr)
            NP=NPLIST(nonode)
            IDATA(1)=NP
          ENDIF
          FORMAT='(/$,'' Node number ['//CHAR1(1:6)//']: '',I6)'
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NPM,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            np=IDATA(1)+N_OFFSET
C GMH 8/1/97 Update cmgui link
            CALL NODE_CREATE(np,ERROR,*9999)
          ENDIF

          CALL ASSERT(np.LE.NPM,'>>Increase NPM',ERROR,*9999)
C CS 28/6/99 At present this is never used.
C
C          IF(JTYP2.EQ.4) THEN !Number of derivatives depends on the node.
C          DO nj=1,NJT
C            WRITE(CHAR1,'(I1)') nj
C            IF(IOTYPE.EQ.3) THEN
C              IF(TYPE(1:9).EQ.'REFERENCE') THEN
C                NKJP(nj)=NKJ(nj,NPNODE(nonode,nr))
C              ELSE IF(TYPE(1:9).EQ.'DEPENDENT') THEN
C                NKJP(nj)=NKH(nj,NPNODE(nonode,nr),nc)
c              ENDIF
c              IDATA(1)=NKJP(nj)-1
c            ENDIF
c            FORMAT='($,'' The number of derivatives for coordinate '
c     '        //CHAR1(1:1)//' is [0]: '',I1)'
c            CALL GINOUT(IOTYPE,3,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
c     '        1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,NKM-1,
c     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
c            IF(IOTYPE.NE.3) NKJP(nj)=IDATA(1)+1
C            ENDDO
C        ENDIF
          IF(TYPE(1:9).EQ.'DEPENDENT') THEN
            NHP(np)=NJT
          ENDIF
          NJTOT=NJT
          DO nj=1,NJTOT

!news MPN 16-Nov-94
            IF(PROMPT_NV(nj)) THEN  !prompt for different versions of nj
              WRITE(CHAR1,'(I1)') nj
              FORMAT='($,'' The number of versions for nj='
     '          //CHAR1(1:1)//' is [1]: '',I2)'
              IF(IOTYPE.EQ.3) THEN
                IF(TYPE(1:9).EQ.'REFERENCE') THEN
                  IDATA(1)=NVJP(nj,np)
                ELSE IF(TYPE(1:9).EQ.'DEPENDENT') THEN
                  IDATA(1)=NVHP(nj,np,nc)
                ENDIF
              ENDIF
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,
     &          NVM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &          *9999)
              IF(IOTYPE.NE.3) THEN
                CALL ASSERT(IDATA(1).LE.NVM,'>>Increase NVM',
     '            ERROR,*9999)
                IF(TYPE(1:9).EQ.'REFERENCE') THEN
                  NVJP(nj,np)=IDATA(1)
                ELSE IF(TYPE(1:9).EQ.'DEPENDENT') THEN
                  NVHP(nj,np,nc)=IDATA(1)
                ENDIF
              ENDIF
            ELSE
              IF(TYPE(1:9).EQ.'REFERENCE') THEN
                NVJP(nj,np)=1
              ELSE IF(TYPE(1:9).EQ.'DEPENDENT') THEN
                NVHP(nj,np,nc)=1
              ENDIF
            ENDIF

            IF(TYPE(1:9).EQ.'REFERENCE') THEN
              NV_MAX=NVJP(nj,np)
            ELSE IF(TYPE(1:9).EQ.'DEPENDENT') THEN
              NV_MAX=NVHP(nj,np,nc)
            ENDIF
            DO nv=1,NV_MAX
              IF(NV_MAX.GT.1) THEN !ask for diff nj versions
                WRITE(CHAR1,'(I2)') nv
                FORMAT='('' For version number'//CHAR1(1:2)//':'')'
                CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,1,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '            0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '            ERROR,*9999)
c               CALL INOUT(IOTYPE,IVDU,IFILE,FORMAT,ERROR,*9999)
              ENDIF
!newe
              WRITE(CHAR1,'(I1)') nj
C KAT 2004-01-27: nb unused.
C               IF(ELEM) THEN
C                 !Note assumption that nb comes from element NEELEM(1,nr)
C                 IF(TYPE(1:9).EQ.'REFERENCE') THEN
C                   nb=NBJ(nj,NEELEM(1,nr))
C                 ELSE IF(TYPE(1:9).EQ.'DEPENDENT') THEN
C                   nb=NBH(nj,nc,NEELEM(1,nr))
C                 ENDIF
C               ELSE
C                 nb=1
C               ENDIF
              IF(IOTYPE.NE.3) THEN
                IF(TYPE(1:9).EQ.'REFERENCE') THEN
                  NKJ(nj,np)=NKJP(nj)
                ELSE IF(TYPE(1:9).EQ.'DEPENDENT') THEN
                  NKH(nj,np,nc)=NKJP(nj)
                ENDIF
              ENDIF

C KAT 2004-01-23: GINOUT should ensure this assert is true.
C C LKC 10-SEP-97 add assert
C               CALL ASSERT(NKM.GE.NKJP(nj),'>>Increase NKM',ERROR,*9999)
              DO nk=1,NKJP(nj)
C!!! n^2 operation
                NP_EXISTS=.FALSE.
                IF(IOTYPE.NE.3) THEN !For backward compatibility in defaults 
                  DO no_region=1,NRT
                    DO n1node=1,NPNODE(0,no_region)
                      NP_PREV=NPNODE(n1node,no_region)
                      IF(np.EQ.NP_PREV) THEN
                        NP_EXISTS=.TRUE.
                        GO TO 15
                      ENDIF
                    ENDDO
                  ENDDO
                ENDIF !IOTYPE.NE.3
 15             IF(NP_EXISTS) THEN
                  NP2=NP
                ELSE
                  IF(N2NODE.EQ.0) THEN
                    NP2=0
                  ELSE
                    NP2=NPNODE(N2NODE,nr)
C                   Only useful if it has enough versions and derivatives
                    IF(nv.GT.NVJP(nj,NP2).OR.nk.GT.NKJ(nj,NP2)) NP2=0
                  ENDIF
                ENDIF
                IF(NP2.EQ.0) THEN
                  RDEFLT(1)=0.0d0
                ELSE
!convert angles into radians
                  IF((ITYP10(nr).GE.2.AND.NJ.EQ.2).OR.
     '              (ITYP10(nr).GE.3.AND.NJ.EQ.3)) THEN
                    IF(TYPE(1:9).EQ.'REFERENCE') THEN
                      RDEFLT(1)=XP(nk,nv,nj,np2)*180.0d0/PI
                    ELSE IF(TYPE(1:9).EQ.'DEPENDENT') THEN
                      RDEFLT(1)=ZP(nk,nv,nj,np2,nc)*180.0d0/PI
                    ENDIF
                  ELSE
                    IF(TYPE(1:9).EQ.'REFERENCE') THEN
                      RDEFLT(1)=XP(nk,nv,nj,np2)
                    ELSE IF(TYPE(1:9).EQ.'DEPENDENT') THEN
                      RDEFLT(1)=ZP(nk,nv,nj,np2,nc)
                    ENDIF
                  ENDIF
                ENDIF
                WRITE(CHAR4,'(E12.5)') RDEFLT(1)
!                IF(JTYP2.NE.4) THEN !always (see comment on else)
                  nu=NUK(nk)
                  IF(nu.EQ.1) THEN
                    FORMAT='($,'' The Xj('//CHAR1(1:1)//') coordinate'//
     '                ' is ['//CHAR4(1:12)//']: '',G25.17)'
                  ELSE IF(nu.EQ.2.OR.nu.EQ.4.OR.nu.EQ.7) THEN
                    IF(nu.EQ.2) THEN
                      CHAR2='1'
                    ELSE IF(nu.EQ.4) THEN
                      CHAR2='2'
                    ELSE IF(nu.EQ.7) THEN
                      CHAR2='3'
                    ENDIF
                    FORMAT='($,'' The derivative wrt direction '//
     '                CHAR2(1:1)//' is ['//CHAR4(1:12)//']: '',G25.17)'
                  ELSE IF(nu.EQ.6.OR.nu.EQ.9.OR.nu.EQ.10) THEN
                    IF(nu.EQ.6) THEN
                      CHAR2='1'
                      CHAR3='2'
                    ELSE IF(nu.EQ.9) THEN
                      CHAR2='1'
                      CHAR3='3'
                    ELSE IF(nu.EQ.10) THEN
                      CHAR2='2'
                      CHAR3='3'
                    ENDIF
                    FORMAT='($,'' The derivative wrt directions '//
     '                CHAR2(1:1)//' & '//CHAR3(1:1)//
     '                ' is ['//CHAR4(1:12)//']: '',G25.17)'
                  ELSE IF(nu.EQ.11) THEN
                    FORMAT='($,'' The derivative wrt directions '//
     '                '1, 2 & 3 is ['//CHAR4(1:12)//']: '',G25.17)'
                  ENDIF
C                ELSE !At present this is never used.
C                 This is only here in case someone may want this someday.
C                  IF(nk.EQ.1) THEN
C                    FORMAT='($,'' The Xj('//CHAR1(1:1)//') coordinate'//
C     '                ' is ['//CHAR4(1:12)//']: '',G25.17)'
C                  ELSE
C                    CHAR2='1'
C                    IF(nk.EQ.2) THEN
C                      CHAR3='1st'
C                    ELSE IF(nk.EQ.3) THEN
C                      CHAR3='2nd'
C                    ELSE IF(nk.EQ.4) THEN
C                      CHAR3='3rd'
C                    ELSE IF(nk.EQ.5) THEN
C                      CHAR3='4th'
C                    ENDIF
C                    FORMAT='($,'' The '//CHAR3(1:3)//' Xj('
C     '                //CHAR1(1:1)//') derivative wrt s('//CHAR2(1:1)
C     '                //') is ['//CHAR4(1:12)//']: '',G25.17)'
C                  ENDIF
!                ENDIF

C               only increment NPNODE(0,0) if the node
C               doesn't exist in any region and nj=nv=nk=1
                IF(IOTYPE.NE.3) THEN
                  IF(.NOT.NP_EXISTS.AND.NK.EQ.1.AND.nv.EQ.1.
     '              AND.NJ.EQ.1)THEN
                    NPNODE(0,0)=NPNODE(0,0)+1
C                   CMGUI link - add the node
                    CALL NODE_CREATE(np,ERROR,*9999)
                  ENDIF
                ELSE !IF(IOTYPE.EQ.3) THEN
                  IF((nj.EQ.2.AND.ITYP10(nr).GT.1)
     '              .OR.(nj.EQ.3.AND.ITYP10(nr).GT.2)) THEN
                    IF(TYPE(1:9).EQ.'REFERENCE') THEN
                      RDATA(1)=XP(nk,nv,nj,np)*180.0d0/PI
                    ELSE IF(TYPE(1:9).EQ.'DEPENDENT') THEN
                      RDATA(1)=ZP(nk,nv,nj,np,nc)*180.0d0/PI
                    ENDIF
                  ELSE
                    IF(TYPE(1:9).EQ.'REFERENCE') THEN
                      RDATA(1)=XP(nk,nv,nj,np)
                    ELSE IF(TYPE(1:9).EQ.'DEPENDENT') THEN
                      RDATA(1)=ZP(nk,nv,nj,np,nc)
                    ENDIF
                  ENDIF
                ENDIF
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,
     '            ICHAR,IDATA,IDEFLT,0,IMAX,LDATA,LDEFLT,
     '            RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
                IF(IOTYPE.NE.3) THEN
                  IF(TYPE(1:9).EQ.'REFERENCE') THEN
                    IF(ADD.AND.(USE_VORONOI.EQ.1)) THEN
C RGB 29/4/98         Calculate node velocities at nj=NJT+1..2*NJT
                      CALL ASSERT(NJM.GE.(2*NJT),
     '                  'Increase NJM to dimension*2',ERROR,*9999)
                      XP(nk,nv,nj+NJT,np)=(RDATA(1)-XP(nk,nv,nj,np))/DT
                      XP(nk,nv,nj,np)=RDATA(1)
                    ELSE
                      XP(nk,nv,nj,np)=RDATA(1)
                    ENDIF
                  ELSE IF(TYPE(1:9).EQ.'DEPENDENT') THEN
                    ZP(nk,nv,nj,np,nc)=RDATA(1)
                  ENDIF
                ENDIF
              ENDDO !nk
            ENDDO !nv
          ENDDO !nj

          IF(IOTYPE.NE.3) THEN
            !translate all angles into radians
            IF(ITYP10(nr).GE.2) THEN
              njj=2
              IF(TYPE(1:9).EQ.'REFERENCE') THEN
                DO nv=1,NVJP(njj,np)
                  DO nk=1,NKJ(njj,np)
                    XP(nk,nv,njj,np)=XP(nk,nv,njj,np)*PI/180.0d0
                  ENDDO
                ENDDO
              ELSE IF(TYPE(1:9).EQ.'DEPENDENT') THEN
                DO nv=1,NVHP(njj,np,nc)
                  DO nk=1,NKH(njj,np,nc)
                    ZP(nk,nv,njj,np,nc)=ZP(nk,nv,njj,np,nc)*PI/180.0d0
                  ENDDO
                ENDDO
              ENDIF
            ENDIF
            IF(ITYP10(nr).GE.3) THEN
              njj=3
              IF(TYPE(1:9).EQ.'REFERENCE') THEN
                DO nv=1,NVJP(njj,np)
                  DO nk=1,NKJ(njj,np)
                    XP(nk,nv,njj,np)=XP(nk,nv,njj,np)*PI/180.0d0
                  ENDDO
                ENDDO
              ELSE IF(TYPE(1:9).EQ.'DEPENDENT') THEN
                DO nv=1,NVHP(njj,np,nc)
                  DO nk=1,NKH(njj,np,nc)
                    ZP(nk,nv,njj,np,nc)=ZP(nk,nv,njj,np,nc)*PI/180.0d0
                  ENDDO
                ENDDO
              ENDIF
            ENDIF

            !ensure all angles are non-negative
            IF(ITYP10(nr).GE.2) THEN
              IF(TYPE(1:9).EQ.'REFERENCE') THEN
                njj=2
                DO nv=1,NVJP(njj,np)
                  IF(XP(1,nv,njj,np).LT.0.0d0) THEN
                    XP(1,nv,njj,np)=XP(1,nv,njj,np)+2.0d0*PI
                  ENDIF
                ENDDO
                IF(ITYP10(nr).GE.3) THEN
                  njj=3
                  DO nv=1,NVJP(njj,np)
                    IF(XP(1,nv,njj,np).LT.0.0d0) THEN
                      XP(1,nv,njj,np)=XP(1,nv,njj,np)+2.0d0*PI
                    ENDIF
                  ENDDO
                ENDIF
              ELSE IF(TYPE(1:9).EQ.'DEPENDENT') THEN
!old MPN 9-Feb-94: OK to have -'ve anlges for deformed coords
c                njj=2
c                DO nv=1,NVHP(njj,np,nc)
c                  IF(ZP(1,nv,njj,NP,nc).LT.0.0d0) THEN
c                    ZP(1,nv,njj,NP,nc)=ZP(1,nv,njj,NP,nc)+2.0D0*PI
c                  ENDIF
c                ENDDO
c                IF(ITYP10(nr).GE.3) THEN
c                  njj=3
c                  DO nv=1,NVHP(njj,np,nc)
c                    IF(ZP(1,nv,njj,np,nc).LT.0.0d0) THEN
c                      ZP(1,nv,njj,np,nc)=ZP(1,nv,njj,np,nc)+2.0D0*PI
c     '              ENDIF
c                  ENDDO
c                ENDIF
              ENDIF
            ENDIF
          ENDIF

          IF(IOTYPE.NE.3) THEN
C!!! n^2 operation
            EXISTS=.FALSE.
            DO n1node=1,NPNODE(0,nr)
              IF(NPNODE(n1node,nr).EQ.NP) THEN
                EXISTS=.TRUE.
                N2NODE=n1node
              ENDIF
            ENDDO
            IF(.NOT.EXISTS) THEN
              NPNODE(0,nr)=NPNODE(0,nr)+1
              CALL ASSERT(NPNODE(0,nr).LE.NP_R_M,'>>Increase NP_R_M',
     '          ERROR,*9999)
              NPNODE(NPNODE(0,nr),nr)=NP
              N2NODE=NPNODE(0,nr)
            ENDIF
          ELSE ! IF(IOTYPE.EQ.3) THEN
C           Default is last element for backward compatibility
C           (but don't need a default at all).
            N2NODE=nonode
          ENDIF !IOTYPE ne/eq 3
        ENDDO !End of nonode loop

        IF(IOTYPE.NE.3) THEN
          CALL ISORT(NPNODE(0,nr),NPNODE(1,nr))
          NPT(nr)=0
          DO nonode=1,NPNODE(0,nr)
            IF(NPNODE(nonode,nr).GT.NPT(nr)) NPT(nr)=NPNODE(nonode,nr)
          ENDDO

C         Setup regions that interface nodes belong to
          CALL INTERFACE(NP_INTERFACE,NPNODE,ERROR,*9999)
        ENDIF

C KAT 18Jun98: Element based variables should be input in ipelem
C      ELSE IF(NBJ(1,1).NE.0) THEN
C        IF(NAT(NBJ(1,1)).GT.0) THEN  !Note: need separate NBH etc
C          DO nq=1,NQT
C            DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
C              WRITE(CHAR1,'(I1)') nj
C              CALL STRING_TRIM(CHAR1,IBEG,IEND)
C              FORMAT='(/$,'' the Xj('//CHAR1(IBEG:IEND)//
C     '          ') coordinate is [1]: '',I2)'
C              IF(nq.EQ.1) THEN
C                RDEFLT(1)=0.0d0
C              ELSE IF(nq.GT.1) THEN
C                RDEFLT(1)=XA(1,nj,NQ-1)
C              ENDIF
C              IF(IOTYPE.EQ.3) RDATA(1)=XA(1,nj,nq)
C              CALL GINOUT(IOTYPE,5,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
C     '          1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
C     '          LDATA,LDEFLT,RDATA,RONE,-RMAX,RMAX,INFO,ERROR,*9999)
C              IF(IOTYPE.NE.3) XA(1,nj,nq)=RDATA(1)
C            ENDDO
C          ENDDO
C        ENDIF
      ENDIF
      IF(NPNODE(0,nr).GT.0) THEN
        CALL REGION_CREATE(nr,ERROR,*9999)
      ENDIF
C KAT 2004-01-23: This doesn't seem to handle adding of old and new
C nodes, which is now handled above.
C C RGB 17/4/98 - if adding new nodes, then if the old number of nodes is
C C greater than the new number of nodes, then this should be the maximum.
C       IF(ADD) THEN
C         NPNODE(0,nr)=MAX(NTNODE,OLD_NTNODE)
C         DO nonode=1,NPNODE(0,nr)
C           IF(NPNODE(nonode,nr).GT.NPT(nr)) NPT(nr)=NPNODE(nonode,nr)
C         ENDDO
C       ENDIF

      CALL EXITS('IPNODE')
      RETURN
 9999 CALL ERRORS('IPNODE',ERROR)
      CALL EXITS('IPNODE')
      RETURN 1
      END

