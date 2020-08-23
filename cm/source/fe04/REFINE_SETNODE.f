      SUBROUTINE REFINE_SETNODE(i1,i2,i3,IBT,IDO,IDRN,INP,NBJ,ne,
     '  nenew,nj,NKJ,nn_ne,nn_nenew,NNIP,NODE,NPNE,
     '  nr,NUMI1,NUNK,NVJE,NVJP,SE,SP,
     '  XE,XI,XII,XP,EXISTS,NOCROSS,ERROR,*)

C#### Subroutine: REFINE_SETNODE
C###  Description:
C###    REFINE_SETNODE sets all properties of the refined node. If
C###    the node has been created or needs another version all values
C###    of XP are set.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grou00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER i1,i2,i3,IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),IDRN,
     '  INP(NNM,NIM,NBFM),NBJ(NJM,NEM),ne,
     '  nenew,NKJ(NJM,NPM),nn_ne,nn_nenew,NNIP(4,4,4),
     '  NODE,NPNE(NNM,NBFM,NEM),nr,
     '  NUMI1,NUNK(NKM,NJM,NPM),NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM)
      REAL*8 SE(NSM,NBFM,NEM),SP(NKM,NBFM,NPM),
     '  XE(NSM,NJM),XI(3),XII,XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
      LOGICAL EXISTS,NOCROSS
!     Local Variables
      INTEGER nb1,nb2,nb3,ni,NITE,nj,nj1,njj1,njj2,nk,nk2,nk3,NKJCHEK,
     '  NKJNEED,NKSTART,nn,nn1,nn2,nogrno,no_grlist,np1,np2,ns,nu,nv
      REAL*8 AA,DIFF,DSDXI(3),G1,G3,LLOOSE_TOL,PXI,ref_fact,R,RC,RR,RRC,
     '  SCALEFACTOR(3),SLX,SMX,SUM,X(11,60),XS(8),XTEMP
      LOGICAL CALCDSDXI,FOUND_LINEAR_BASIS,
     '  HIGHERDERIV,NEWVERSION,USED
      INTEGER*4 ILIST_PTR
!     Functions
      LOGICAL INLIST
      INTEGER*4 ILISTLOC

C      DATA IDONK/0,1,1,2,1,2,2,3/ !deriv order

C MPN 22Apr2009: increased number of components in X (for Yikun Wang)
C       and added this ASSERT to check nj does not exceed this number.
      CALL ASSERT(nj.LE.60,'>>ERROR: array X cannot handle nj>60.'
     '  //'Reduce the number of fields in your problem.',ERROR,*9999)

      CALL ENTERS('REFINE_SETNODE',*9999)

      LLOOSE_TOL=DSQRT(LOOSE_TOL)

C*** Determine if we need to create new nodal information. Need to
C*** create new node information if (a) the node didn't exist,
C*** (b) if we need to create a new version at the node or
C*** (c) the basis functions in the new element have a higher number
C*** of derivatives than the element that previously created the node.

C KAT 9Dec98: New approach to versions
C      SETPROPERTIES=.FALSE.
C      nn1=NNIP(1,i2,i3)
C      nn2=NNIP(NUMI1,i2,i3)
C      DO njj1=1,3
C        DO njj2=1,NJ_LOC(njj1,0,nr)
C          nj=NJ_LOC(njj1,njj2,nr)
C
C          nb1=NBJ(nj,ne)
C
CC*** Determine if the basis functions in the new and old elements have
CC*** higher number of derivatives than the current number of global
CC*** derivatives defined at the node (if it exists)
C          IF(EXISTS) THEN
C            HIGHERDERIV(nj)=NKT(nn1,nb1).GT.NKJ(nj,NODE).OR.
C     '        NKT(nn2,nb1).GT.NKJ(nj,NODE)
C          ELSE
C            HIGHERDERIV(nj)=.FALSE.
C          ENDIF
C
CC*** Determine if we need to create a new version at the node. If the
CC*** node is on a 'line' between two nodes in the direction of IDR1
CC*** in an elements of which any of the surrounding elements which
CC*** also contain those nodes use a different version at each end of
CC*** the line then the node will have multiple versions.
C
C          IF(nb1.GT.0) THEN
C            np1=NPNE(NNIP(i1,i2,i3),nb1,ne)
C            np2=NPNE(NNIP(i1+1,i2,i3),nb1,ne)
CC*** Find the list of adjacent elements i.e. elements that contain
CC*** both np1 and np2 or elements that contain np1 and NODE or NODE
CC*** and np2 or elements that contain at least two matching nodes
CC*** which also contain np1 or np2
C            NEELIST(0,1,nj)=0
C            DO noelem=1,NEELEM(0,nr)
C              nee=NEELEM(noelem,nr)
C              IF(nee.NE.ne) THEN
C                matchnb=NBJ(nj,nee)
C                nummatch=0
C                FOUNDNP1=.FALSE.
C                FOUNDNP2=.FALSE.
C                DO nn=1,NNT(matchnb)
C                  IF(NPNE(nn,matchnb,nee).EQ.np1.AND.
C     '              .NOT.FOUNDNP1) THEN
C                    FOUNDNP1=.TRUE.
C                  ELSE IF(NPNE(nn,matchnb,nee).EQ.np2.AND.
C     '                .NOT.FOUNDNP2) THEN
C                    FOUNDNP2=.TRUE.
C                  ENDIF
C                ENDDO !nn
C                IF(FOUNDNP1.AND.FOUNDNP2) THEN !matching element
C                  NEELIST(0,1,nj)=NEELIST(0,1,nj)+1
C                  IF(NEELIST(0,1,nj).LE.100) THEN
C                    NEELIST(NEELIST(0,1,nj),1,nj)=nee
C                    NEELIST(NEELIST(0,1,nj),2,nj)=nee
C                  ELSE
C                    ERROR='>>Increase length of NEELIST'
C                    GOTO 9999
C                  ENDIF
C                ELSE
C                  FOUNDNP1=.FALSE.
C                  FOUNDNP2=.FALSE.
C                  DO nn=1,NNT(matchnb)
C                    IF(NPNE(nn,matchnb,nee).EQ.np1.AND.
C     '                .NOT.FOUNDNP1) THEN
C                      FOUNDNP1=.TRUE.
C                    ELSE IF(NPNE(nn,matchnb,nee).EQ.NODE.AND.
C     '                  .NOT.FOUNDNP2) THEN
C                      FOUNDNP2=.TRUE.
C                    ENDIF
C                  ENDDO !nn
C                  IF(.NOT.(FOUNDNP1.AND.FOUNDNP2)) THEN
C                    nummatch=0
C                    FOUNDNP1=.FALSE.
C                    DO nn=1,NNT(matchnb)
C                      DO i=1,NUMI1
C                        IF(NPNE(nn,matchnb,nee).EQ.
C     '                    NPNE(NNIP(i,i2,i3),nb1,ne)) THEN
C                          nummatch=nummatch+1
C                        ELSE IF(NPNE(nn,matchnb,nee).EQ.np1.AND.
C     '                      .NOT.FOUNDNP1) THEN
C                          FOUNDNP1=.TRUE.
C                        ENDIF
C                      ENDDO !i
C                    ENDDO !nn
C                    FOUNDNP2=nummatch.GE.2
C                  ENDIF
C                  IF(FOUNDNP1.AND.FOUNDNP2) THEN !matching element
C                    NEELIST(0,1,nj)=NEELIST(0,1,nj)+1
C                    IF(NEELIST(0,1,nj).LE.100) THEN
C                      NEELIST(NEELIST(0,1,nj),1,nj)=nee
C                      NEELIST(NEELIST(0,1,nj),2,nj)=0
C                      DO noelem2=1,NEELEM(0,nr)
C                        IF(NEELIST(NEELIST(0,1,nj),2,nj).EQ.0) THEN
C                          nee2=NEELEM(noelem2,nr)
C                          IF(nee2.NE.ne) THEN
C                            matchnb=NBJ(nj,nee2)
C                            FOUNDNP1=.FALSE.
C                            FOUNDNP2=.FALSE.
C                            DO nn=1,NNT(matchnb)
C                              IF(NPNE(nn,matchnb,nee2).EQ.NODE.AND.
C     '                          .NOT.FOUNDNP1) THEN
C                                FOUNDNP1=.TRUE.
C                              ELSE IF(NPNE(nn,matchnb,nee2).EQ.np2.AND.
C     '                            .NOT.FOUNDNP2) THEN
C                                FOUNDNP2=.TRUE.
C                              ENDIF
C                            ENDDO !nn
C                            IF(.NOT.(FOUNDNP1.AND.FOUNDNP2)) THEN
C                              nummatch=0
C                              FOUNDNP2=.FALSE.
C                              DO nn=1,NNT(matchnb)
C                                DO i=1,NUMI1
C                                  IF(NPNE(nn,matchnb,nee2).EQ.
C     '                              NPNE(NNIP(i,i2,i3),nb1,ne)) THEN
C                                    nummatch=nummatch+1
C                                  ELSE IF(NPNE(nn,matchnb,nee2).EQ.np2
C     '                                .AND..NOT.FOUNDNP2) THEN
C                                    FOUNDNP2=.TRUE.
C                                  ENDIF
C                                ENDDO !i
C                              ENDDO !nn
C                              FOUNDNP1=nummatch.GE.2
C                            ENDIF
C                            IF(FOUNDNP1.AND.FOUNDNP2) THEN !matching element
C                              NEELIST(NEELIST(0,1,nj),2,nj)=nee2
C                            ENDIF
C                          ENDIF
C                        ENDIF
C                      ENDDO !neelem2
C                      CALL ASSERT(NEELIST(NEELIST(0,1,nj),2,nj).NE.0,
C     '                  '>>Could not find matching element. '
C     '                  //'Hanging node?',ERROR,*9999)
C                    ELSE
C                      ERROR='>>Increase length of NEELIST'
C                      GOTO 9999
C                    ENDIF
C                  ENDIF
C                ENDIF
C              ENDIF
C            ENDDO !neelem
C            NEWVERSION(nj)=.FALSE.
C            current_nv1=NVJE(NNIP(i1,i2,i3),nb1,nj,ne)
C            current_nv2=NVJE(NNIP(i1+1,i2,i3),nb1,nj,ne)
C            NVLIST(0,1,nj)=1
C            NVLIST(1,1,nj)=ne
C            NVLIST(1,2,nj)=MAX(current_nv1,current_nv2)
C            IF(DOP) THEN
CC$            call mp_setlock()
C              WRITE(OP_STRING,'('' nj='',I2,'', nb1='',I2,'
C     '          //''', #elements='',I3)') nj,nb1,NEELIST(0,1,nj)
C              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C              WRITE(OP_STRING,'('' np1='',I5,'', current_nv1='',I2,'
C     '          //''', np2='',I5'', current_nv2='',I2)')
C     '          np1,current_nv1,np2,current_nv2
C              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
C            ENDIF
C            DO noelem=1,NEELIST(0,1,nj)
C              nee=NEELIST(noelem,1,nj)
C              nee2=NEELIST(noelem,2,nj)
C              nbb=NBJ(nj,nee)
C              nbb2=NBJ(nj,nee2)
C              IF(nbb.NE.0.AND.nbb.NE.0) THEN
C                FOUND=.FALSE.
C                nn1=1
C                DO WHILE(.NOT.FOUND.AND.nn1.LE.NNT(nbb))
C                  IF(NPNE(nn1,nbb,nee).EQ.np1) THEN
C                    FOUND=.TRUE.
C                  ELSE
C                    nn1=nn1+1
C                  ENDIF
C                ENDDO !nn
C                CALL ASSERT(FOUND,'>>Could not find local node',
C     '            ERROR,*9999)
C                FOUND=.FALSE.
C                nn2=1
C                DO WHILE(.NOT.FOUND.AND.nn2.LE.NNT(nbb2))
C                  IF(np1.EQ.np2) THEN
C                    IF(NPNE(nn2,nbb2,nee2).EQ.np2.AND.nn1.NE.nn2) THEN
C                      FOUND=.TRUE.
C                    ELSE
C                      nn2=nn2+1
C                    ENDIF
C                  ELSE
C                    IF(NPNE(nn2,nbb2,nee2).EQ.np2) THEN
C                      FOUND=.TRUE.
C                    ELSE
C                      nn2=nn2+1
C                    ENDIF
C                  ENDIF
C                ENDDO !nn
C                CALL ASSERT(FOUND,'>>Could not find local node',
C     '            ERROR,*9999)
C                nv1=NVJE(nn1,nbb,nj,nee)
C                nv2=NVJE(nn2,nbb2,nj,nee2)
C                IF(nv1.NE.current_nv1.OR.nv2.NE.current_nv2) THEN
C                  NEWVERSION(nj)=.TRUE.
C                  nv=NVLIST(0,1,nj)+1
C                  DO nvv=1,NVLIST(0,1,nj)
C                    IF(NVLIST(nvv,2,nj).GT.MAX(nv1,nv2)) nv=nvv
C                  ENDDO !nvv
C                  NVLIST(0,1,nj)=NVLIST(0,1,nj)+1
C                  IF(NVLIST(0,1,nj).LE.250) THEN
C                    IF(nv.LT.NVLIST(0,1,nj)) THEN
C                      DO nvv=NVLIST(0,1,nj),nv+1,-1
C                        NVLIST(nvv,1,nj)=NVLIST(nvv-1,1,nj)
C                        NVLIST(nvv,2,nj)=NVLIST(nvv-1,2,nj)
C                      ENDDO
C                      NVLIST(nv,1,nj)=nee
C                      NVLIST(nv,2,nj)=MAX(nv1,nv2)
C                    ELSE
C                      NVLIST(NVLIST(0,1,nj),1,nj)=nee
C                      NVLIST(NVLIST(0,1,nj),2,nj)=MAX(nv1,nv2)
C                    ENDIF
C                  ELSE
C                    ERROR='>>Increase length of NVLIST'
C                    GOTO 9999
C                  ENDIF
C                ENDIF
C                IF(DOP) THEN
CC$                call mp_setlock()
C                  WRITE(OP_STRING,'('' nee ='',I5,'', nbb ='',I2,'
C     '              //''', nn1='',I2,'', nv1='',I2)') nee,nbb,nn1,nv1
C                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                  WRITE(OP_STRING,'('' nee2='',I5,'', nbb2='',I2,'
C     '              //''', nn2='',I2,'', nv2='',I2)') nee2,nbb2,nn2,nv2
C                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                  WRITE(OP_STRING,'('' NEWVERSION(nj)='',L1)')
C     '              NEWVERSION(nj)
C                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                  WRITE(OP_STRING,'('' NVLIST(0,1,nj)='',I2)')
C     '              NVLIST(0,1,nj)
C                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                  WRITE(OP_STRING,'('' NVLIST(1..,1,nj)='',20I6)')
C     '              (NVLIST(nv,1,nj),nv=1,NVLIST(0,1,nj))
C                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C                  WRITE(OP_STRING,'('' NVLIST(1..,2,nj)='',20I6)')
C     '              (NVLIST(nv,2,nj),nv=1,NVLIST(0,1,nj))
C                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
C                ENDIF
C              ENDIF
C            ENDDO !noelem
CC*** Check adjacent refined elements
C          ELSE
C            NEWVERSION(nj)=.FALSE.
C          ENDIF !nb1 > 0
C        ENDDO !njj2
C      ENDDO !njj1

      NITE=NIT(NBJ(1,ne))
      nn1=NNIP(1,i2,i3)
      nn2=NNIP(NUMI1,i2,i3)
C KAT 22Jul99: Only calculating things needed for this nj
CC     Calculate zeroth derivatives and
C     See if arc length scale factors are needed.
      CALCDSDXI=.FALSE.
C      DO njj1=1,3
C        DO njj2=1,NJ_LOC(njj1,0,nr)
C          nj=NJ_LOC(njj1,njj2,nr)
      nb1=NBJ(nj,ne)
      IF(nb1.GT.0) THEN
C            X(1,nj)=PXI(IBT(1,1,nb1),IDO(1,1,0,nb1),INP(1,1,nb1),nb1,
C     '        1,XI,XE(1,nj))
        IF(NBI(nb1).GE.5.AND.NBI(nb1).LE.7) CALCDSDXI=.TRUE.
      ENDIF !nb1>0
C        ENDDO !njj2
C      ENDDO !njj1
      IF(CALCDSDXI) THEN
        IF(ITYP10(nr).GE.2) THEN
          nb1=NBJ(1,ne)
          X(1,1)=PXI(IBT(1,1,nb1),IDO(1,1,0,nb1),INP(1,1,nb1),nb1,
     '        1,XI,XE(1,1))
C         If an angle is 2*pi set it to zero
          IF(ITYP10(nr).EQ.3) THEN
            nb1=NBJ(3,ne)
            X(1,3)=PXI(IBT(1,1,nb1),IDO(1,1,0,nb1),INP(1,1,nb1),nb1,
     '        1,XI,XE(1,3))
            IF(DABS(X(1,3)-2.0d0*PI).LT.LOOSE_TOL) X(1,3)=0.0d0
          ENDIF
          IF(ITYP10(nr).EQ.4) THEN
            nb1=NBJ(2,ne)
            X(1,2)=PXI(IBT(1,1,nb1),IDO(1,1,0,nb1),INP(1,1,nb1),nb1,
     '        1,XI,XE(1,2))
            IF(DABS(X(1,2)-2.0d0*PI).LT.LOOSE_TOL) X(1,2)=0.0d0
          ENDIF
        ENDIF
        DO ni=1,NITE
          nu=1+ni*(1+ni)/2
          DO nj1=1,NJT
            nb1=NBJ(nj1,ne)
            X(nu,nj1)=PXI(IBT(1,1,nb1),IDO(1,1,0,nb1),INP(1,1,nb1),nb1,
     '        nu,XI,XE(1,nj1))
          ENDDO !nj1
          SUM=0.0d0
          IF(ITYP10(nr).EQ.1) THEN
            SUM=X(nu,1)**2+X(nu,2)**2
            IF(NJT.GT.2) SUM=SUM+X(nu,3)**2
          ELSE IF(ITYP10(nr).EQ.2) THEN
            R=X(1,1)
            RR=R*R
            SUM=SUM+X(nu,1)**2+RR*X(nu,2)**2
            IF(NJT.GT.2) SUM=SUM+X(nu,3)**2
          ELSE IF(ITYP10(nr).EQ.3) THEN
            R=X(1,1)
            RR=R*R
            RC=R*DCOS(X(1,3))
            RRC=RC*RC
            SUM=SUM+X(nu,1)**2+RRC*X(nu,2)**2+RR*X(nu,3)**2
          ELSE IF(ITYP10(nr).EQ.4) THEN
            AA=FOCUS*FOCUS
            SLX=DSINH(X(1,1))
            SMX=DSIN(X(1,2))
            G1=AA*(SLX*SLX+SMX*SMX)
            G3=AA* SLX*SLX*SMX*SMX
            SUM=SUM+G1*(X(nu,1)**2+X(nu,2)**2)
            IF(NJT.GT.2) SUM=SUM+G3*X(nu,3)**2
          ENDIF
          DSDXI(ni)=DSQRT(SUM)
          IF(DABS(DSDXI(ni)).LT.LOOSE_TOL) THEN
            DSDXI(ni)=1.0d0
            WRITE(OP_STRING,
     '        '('' >>Warning: Arc length derivative is zero'')')
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          ENDIF
        ENDDO !ni
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' DSDXI(ni):'',3(1X,D12.5))')
     '      (DSDXI(ni),ni=1,NITE)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
      ENDIF !CALCDSDXI
C KAT 8Oct99: Initiallizing in REFINE_FINDNODE instead
CC KAT 22Jul99: Only calculating things needed for this nj
CC      DO njj1=1,3
CC        DO njj2=1,NJ_LOC(njj1,0,nr)
CC          nj=NJ_LOC(njj1,njj2,nr)
C          SETPROPERTIES=.NOT.EXISTS
C          IF(SETPROPERTIES) THEN
CC           Initialize for all nj
C            DO njj1=1,3
C              DO njj2=1,NJ_LOC(njj1,0,nr)
C                nj1=NJ_LOC(njj1,njj2,nr)
C                IF(nj1.LE.NJT) THEN
C                  NKJ(nj1,NODE)=1
C                  NUNK(1,nj1,NODE)=1
CC                 NVJP and XP(1,1,nj1,NODE) already set in REFINE_FINDNODE
C                ELSE
C                  NKJ(nj1,NODE)=0
C                  NVJP(nj1,NODE)=0
C                ENDIF
C              ENDDO !njj2
C            ENDDO !njj1
CC            IF(nj.LE.NJT) THEN
CC              NKJ(nj,NODE)=1
CC              NUNK(1,nj,NODE)=1
CCC             NVJP and XP(1,1,nj,NODE) already set in REFINE_FINDNODE
CC            ELSE
CC              NKJ(nj,NODE)=0
CC              NVJP(nj,NODE)=0
CC            ENDIF
C          ENDIF
          nb1=NBJ(nj,ne)
C         Calculate remaining required derivatives
          IF(nb1.GT.0) THEN
C           Number of derivatives needed at the new node
            IF(NKT(nn1,nb1).GT.NKT(nn2,nb1)) THEN
              nn=nn1
            ELSE
              nn=nn2
            ENDIF
            NKJNEED=NKT(nn,nb1)
C           Check how many of these derivatives already at the node are used
            nk=NKJ(nj,NODE)
            USED=.FALSE. !for now
            DO WHILE(nk.GT.0.AND..NOT.USED)
              IF(NUNK(nk,nj,NODE).NE.0) THEN
                USED=.TRUE.
              ELSE
                nk=nk-1
              ENDIF
            ENDDO
C           Number of derivatives that can be compared with existing node
            NKJCHEK=MIN(nk,NKJNEED)
C           One more then number of derivs already at node
            NKSTART=nk+1
C           Check that there are enough derivatives at the new node for
C           both the new element and the old element.
            HIGHERDERIV=NKJNEED.GE.NKSTART
            IF(HIGHERDERIV) THEN
              NKJ(nj,NODE)=NKJNEED
              DO nk=NKSTART,NKJNEED
                NUNK(nk,nj,NODE)=NUNK(nk,nj,NPNE(nn,nb1,ne))
              ENDDO !nk
            ENDIF !HIGHERDERIV
C KAT 22Jul99: Only calculating things needed for this nj
            DO nk=1,NKJNEED
C            DO nk=2,NKJNEED
              nu=NUNK(nk,nj,NODE)
C              IF(.NOT.CALCDSDXI.OR.nj.GT.NJT
C     '          .OR.IDONK(nk).GT.1) THEN !not already calc'd
              X(nu,nj)=PXI(IBT(1,1,nb1),IDO(1,1,0,nb1),INP(1,1,nb1),
     '          nb1,nu,XI,XE(1,nj))
              IF((nj.EQ.2.AND.ITYP10(nr).GE.2).OR.
     '          (nj.EQ.3.AND.ITYP10(nr).GE.3)) THEN
C KAT 18Oct00:  XPXE often modifies nodal angles usually by adding 2pi.
C               Here we make a half hearted attempt to keep them
C               within [0,2pi).
                XTEMP=X(1,nj)-2d0*PI
                IF(XTEMP.GE.0d0) X(1,nj)=XTEMP
CC***  If an angle > 2*pi set it to zero
C                IF(X(1,nj)-2.0d0*PI).LT.LOOSE_TOL) X(1,nj)=0.0d0
              ENDIF
C              ENDIF
            ENDDO !nk
C KAT 2/3/00: It seems that something has changed in
C           REFINE_SETSCALEFACTORS so this now needs to be done even for
C           unit scale factors.
C            IF(NBI(nb1).GT.1) THEN !correct for scale factors
              IF(NBI(nb1).GE.5.AND.NBI(nb1).LE.7) THEN
C               arc-length or ave. arc-length scale factors
                DO ni=1,NITE
                  SCALEFACTOR(ni)=DSDXI(ni)
                ENDDO !ni
C KAT 13Jan00: Lower dimensions now handled in REFINE_SETSCALEFACTORS
C              ELSE IF(NITE.LT.3) THEN
CC??? Not sure about this bit. Calculate derivatives of s* wrt Xi
CC??? Note: Only correct for bicubic at present
C                IF(IDRN.EQ.1) THEN
C                  IF(i3.EQ.1) THEN
C                    IF(i2.EQ.1) THEN
C                      SCALEFACTOR(1)=
C     '                  0.25d0*(SE( 2,nb1,ne)+SE( 6,nb1,ne))
C                      SCALEFACTOR(2)=
C     '                  0.50d0*(SE( 3,nb1,ne)+SE( 7,nb1,ne))
C                    ELSE IF(i2.EQ.2) THEN
C                      SCALEFACTOR(1)=
C     '                  0.25d0*(SE(10,nb1,ne)+SE(14,nb1,ne))
C                      SCALEFACTOR(2)=
C     '                  0.50d0*(SE(11,nb1,ne)+SE(15,nb1,ne))
C                    ENDIF
C                  ELSE IF(i3.EQ.2) THEN
C                    IF(i2.EQ.1) THEN
C                      SCALEFACTOR(1)=
C     '                  0.25d0*(SE(18,nb1,ne)+SE(22,nb1,ne))
C                      SCALEFACTOR(2)=
C     '                  0.50d0*(SE(19,nb1,ne)+SE(23,nb1,ne))
C                    ELSE IF(i2.EQ.2) THEN
C                      SCALEFACTOR(1)=
C     '                  0.25d0*(SE(26,nb1,ne)+SE(30,nb1,ne))
C                      SCALEFACTOR(2)=
C     '                  0.50d0*(SE(27,nb1,ne)+SE(31,nb1,ne))
C                    ENDIF
C                  ENDIF
C                ELSE IF(IDRN.EQ.2) THEN
C                  IF(i2.EQ.1) THEN
C                    IF(i3.EQ.1) THEN
C                      SCALEFACTOR(1)=
C     '                  0.50d0*(SE( 2,nb1,ne)+SE(10,nb1,ne))
C                      SCALEFACTOR(2)=
C     '                  0.25d0*(SE( 3,nb1,ne)+SE(11,nb1,ne))
C                    ELSE IF(i3.EQ.2) THEN
C                      SCALEFACTOR(1)=
C     '                  0.50d0*(SE( 6,nb1,ne)+SE(14,nb1,ne))
C                      SCALEFACTOR(2)=
C     '                  0.25d0*(SE( 7,nb1,ne)+SE(15,nb1,ne))
C                    ENDIF
C                  ELSE IF(i2.EQ.2) THEN
C                    IF(i3.EQ.1) THEN
C                      SCALEFACTOR(1)=
C     '                  0.50d0*(SE(18,nb1,ne)+SE(26,nb1,ne))
C                      SCALEFACTOR(2)=
C     '                  0.25d0*(SE(19,nb1,ne)+SE(27,nb1,ne))
C                    ELSE IF(i3.EQ.2) THEN
C                      SCALEFACTOR(1)=
C     '                  0.50d0*(SE(22,nb1,ne)+SE(30,nb1,ne))
C                      SCALEFACTOR(2)=
C     '                  0.25d0*(SE(23,nb1,ne)+SE(31,nb1,ne))
C                    ENDIF
C                  ENDIF
C                ENDIF !idrn
C                IF(DOP) THEN
CC$                call mp_setlock()
C                  WRITE(OP_STRING,'('' SCALEFACTOR(ni):'',3(1X,D12.5))')
C     '              (SCALEFACTOR(ni),ni=1,NITE)
C                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
C                ENDIF
C                DO ni=1,NITE
C                  IF(DABS(SCALEFACTOR(ni)).LT.LOOSE_TOL) THEN
C                    SCALEFACTOR(ni)=1.0d0
C                    WRITE(OP_STRING,
C     '                '('' >>Warning: Scale factor is zero'')')
C                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C                  ENDIF
C                ENDDO !ni
              ELSE
                DO ni=1,3
                  SCALEFACTOR(ni)=1.0d0
                ENDDO !ni
              ENDIF !nbi
              SP(1,nb1,NODE)=1.0d0
              DO nk=2,NKJNEED
                nu=NUNK(nk,nj,NODE)
                IF(IDO(nk,nn_nenew,IDRN,NBJ(nj,nenew)).EQ.1) THEN
                  ref_fact=1.0d0
                ELSE
                  ref_fact=1.0d0-XII
                ENDIF

C Hanging nodes - taken from HANGING_NODE_DETECT
            IF(JTYP2C.EQ.1) THEN
C For the refined direction set the scale factor from interpolation
            FOUND_LINEAR_BASIS=.FALSE.
            nb2=0
            DO WHILE((.NOT.FOUND_LINEAR_BASIS)
     '        .AND.(nb2.LE.NBFT))
              nb2=nb2+1
              IF(NIT(nb1).EQ.3) THEN
                IF((IBT(1,1,nb2).EQ.1).AND.
     '            (IBT(1,2,nb2).EQ.1).AND.(IBT(1,3,nb2).EQ.1)) THEN
                  FOUND_LINEAR_BASIS=.TRUE.
                ENDIF
              ELSE
                IF((IBT(1,1,nb2).EQ.1).AND.
     '            (IBT(1,2,nb2).EQ.1)) THEN
                  FOUND_LINEAR_BASIS=.TRUE.
                ENDIF
              ENDIF
            ENDDO
            CALL ASSERT(nb2.LE.NBFT,
     '        '>>Define linear basis',ERROR,*9999)

            DO nb3=1,NBFT
              IF(NKT(0,nb3).GT.1.AND.NBI(nb3).NE.1) THEN !scale fac not unit
                DO nk2=2,NKT(1,nb3)
                  ns=0
                  DO nn=1,NNT(nb3)
                    ns=ns+1
                    DO nk3=2,NKT(nn,nb3)
                      ns=ns+1
                      IF(nk3.EQ.nk2) XS(nn)=SE(ns,nb3,ne)
                    ENDDO !nk
                  ENDDO !nn
                  !only set in refined direction
C!!! setting using basis 1 for now
                  IF(IDRN.EQ.1)THEN
                    IF((nk2.EQ.2).AND.(nb3.EQ.1)) THEN
                       SCALEFACTOR(1)=
     '               PXI(IBT(1,1,nb2),IDO(1,1,0,nb2),
     '               INP(1,1,nb2),nb2,1,XI,XS)
                    ENDIF
                  ELSE IF(IDRN.EQ.2) THEN
                    IF((nk2.EQ.3).AND.(nb3.EQ.1)) THEN
                       SCALEFACTOR(2)=
     '               PXI(IBT(1,1,nb2),IDO(1,1,0,nb2),
     '               INP(1,1,nb2),nb2,1,XI,XS)
                    ENDIF
                  ELSE
                    IF((nk2.EQ.5).AND.(nb3.EQ.1)) THEN
                       SCALEFACTOR(3)=
     '               PXI(IBT(1,1,nb2),IDO(1,1,0,nb2),
     '               INP(1,1,nb2),nb2,1,XI,XS)
                    ENDIF
                  ENDIF
                ENDDO !nk2
              ENDIF
            ENDDO !nb3
            ENDIF

                IF(nu.EQ.2) THEN !dXi1
                  X(nu,nj)=X(nu,nj)/SCALEFACTOR(1)
C CS new July 1999 Added the set up of SE for new nodes here.
C KAT 13Jan00:    What is the relationship between SE values set here
C                 and in REFINE_SETSCALEFACTORS?
                  SE(nk+((nn_nenew-1)*NKJNEED),NBJ(nj,nenew),nenew)=
     '              ref_fact*SCALEFACTOR(1)
                  SP(nk,nb1,NODE)=SCALEFACTOR(1)
                ELSE IF(nu.EQ.4) THEN !dXi2
                  X(nu,nj)=X(nu,nj)/SCALEFACTOR(2)
                  SE(nk+((nn_nenew-1)*NKJNEED),NBJ(nj,nenew),nenew)=
     '              ref_fact*SCALEFACTOR(2)
                  SP(nk,nb1,NODE)=SCALEFACTOR(2)
                ELSE IF(nu.EQ.7) THEN !dXi3
                  X(nu,nj)=X(nu,nj)/SCALEFACTOR(3)
                  SE(nk+((nn_nenew-1)*NKJNEED),NBJ(nj,nenew),nenew)=
     '              ref_fact*SCALEFACTOR(3)
                  SP(nk,nb1,NODE)=SCALEFACTOR(3)
                ELSE IF(NOCROSS) THEN
                  X(nu,nj)=0.0d0
                  SE(nk+((nn_nenew-1)*NKJNEED),NBJ(nj,nenew),nenew)=1
                  SP(nk,nb1,NODE)=0.0d0
                ELSE IF(nu.EQ.6) THEN !dXi1dXi2
                  X(nu,nj)=X(nu,nj)/(SCALEFACTOR(1)*SCALEFACTOR(2))
                  SE(nk+((nn_nenew-1)*NKJNEED),NBJ(nj,nenew),nenew)=
     '              ref_fact*(SCALEFACTOR(1)*SCALEFACTOR(2))
                  SP(nk,nb1,NODE)=SCALEFACTOR(1)*SCALEFACTOR(2)
                ELSE IF(nu.EQ.9) THEN !dXi1dXi3
                  X(nu,nj)=X(nu,nj)/(SCALEFACTOR(1)*SCALEFACTOR(3))
                  SE(nk+((nn_nenew-1)*NKJNEED),NBJ(nj,nenew),nenew)=
     '              ref_fact*(SCALEFACTOR(1)*SCALEFACTOR(3))
                  SP(nk,nb1,NODE)=SCALEFACTOR(1)*SCALEFACTOR(3)
                ELSE IF(nu.EQ.10) THEN !dXi2dXi3
                  X(nu,nj)=X(nu,nj)/(SCALEFACTOR(2)*SCALEFACTOR(3))
                  SE(nk+((nn_nenew-1)*NKJNEED),NBJ(nj,nenew),nenew)=
     '              ref_fact*(SCALEFACTOR(2)*SCALEFACTOR(3))
                  SP(nk,nb1,NODE)=SCALEFACTOR(2)*SCALEFACTOR(3)
                ELSE IF(nu.EQ.11) THEN !dXi1dXi2dXi3
                  X(nu,nj)=X(nu,nj)/
     '              (SCALEFACTOR(1)*SCALEFACTOR(2)*SCALEFACTOR(3))
                  SE(nk+((nn_nenew-1)*NKJNEED),NBJ(nj,nenew),nenew)=
     '              ref_fact*(SCALEFACTOR(1)*SCALEFACTOR(2)*
     '              SCALEFACTOR(3))
                  SP(nk,nb1,NODE)=
     '              SCALEFACTOR(1)*SCALEFACTOR(2)*SCALEFACTOR(3)
                ELSE
                  WRITE(OP_STRING,'('' >>Warning: Derivative '
     '              //'not updated. nu='',I2,'' is unknown'')') nu
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDIF
              ENDDO !nk
C            ENDIF !scale factors not unit
C KAT 9Dec98: New approach to versions
C           Determine the version to use by checking if a version with the
C           appropriate information already exists.
            NEWVERSION=.TRUE. !assume a new version unless we find a match
            nv=0
            DO WHILE(nv.LT.NVJP(nj,NODE).AND.NEWVERSION)
              nv=nv+1
              NEWVERSION=.FALSE. !assume node matches until it doesn't
              nk=0
              DO WHILE(nk.LT.NKJCHEK.AND..NOT.NEWVERSION)
                nk=nk+1
                nu=NUNK(nk,nj,NODE)
                DIFF=DABS(X(nu,nj)-XP(nk,nv,nj,NODE))
                IF(DIFF.GT.LLOOSE_TOL*(DABS(X(nu,nj))+LLOOSE_TOL)) THEN
                  NEWVERSION=.TRUE. !doesn't match unless special cases
                  IF(nk.EQ.1) THEN
C                   angles of 0 and 2pi are the same
                    IF((nj.EQ.2.AND.ITYP10(nr).GE.2).OR.
     '                (nj.EQ.3.AND.ITYP10(nr).GE.3)) THEN
                      NEWVERSION=DABS(2d0*PI-DIFF).GT.LOOSE_TOL
                    ENDIF
                  ENDIF
                ENDIF
              ENDDO !nk
            ENDDO !nv
            IF(NEWVERSION) THEN
              nv=nv+1
              NKSTART=1
              CALL ASSERT(nv.LE.NVM,' >>ERROR: Increase NVM',ERROR,
     '          *9999)
              NVJP(nj,NODE)=nv
            ENDIF !NEWVERSION
C KAT 9Dec98: New approach to versions
CC*** Determine the version number from the list of orderer maximum
CC*** version numbers of the adjacent elements.
C                nv=0
C                DO nvv=1,NVLIST(0,1,nj)
C                  IF(NVLIST(nvv,1,nj).EQ.ne) nv=nvv
C                ENDDO !nvv
C                CALL ASSERT(nv.NE.0,
C     '            '>>Could not find element in NVLIST',ERROR,*9999)
C                IF(nv.GT.NVJP(nj,NODE)) NVJP(nj,NODE)=nv
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
              WRITE(OP_STRING,
     '          '('' HIGHERDERIV='',L1,'', NEWVERSION='',L1)')
     '          HIGHERDERIV,NEWVERSION
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' nj='',I2,'', nb1='',I2,'
     '          //''', nv='',I2)') nj,nb1,nv
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' X(nu,nj):'',6(1X,D10.3),'
     '          //'/:(10X,6(1X,D10.3)))') (X(nu,nj),nu=1,NUT(nb1))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
            ENDIF
            DO nk=NKSTART,NKJNEED
              nu=NUNK(nk,nj,NODE)
              XP(nk,nv,nj,NODE)=X(nu,nj)
            ENDDO !nk
            IF(nn_ne.NE.0) NVJE(nn_ne,nb1,nj,ne)=nv
            nb2=NBJ(nj,nenew)
            IF(nb2.NE.0.AND.nn_nenew.NE.0) THEN
              NVJE(nn_nenew,nb2,nj,nenew)=nv
            ENDIF
          ENDIF !nb1>0
C        ENDDO !njj2
C      ENDDO !njj1

C*** NEW CS 28/11/97 If np1 and np2 are in the same group then the new
C*** node NODE is made a member of the group
      IF(.NOT.EXISTS) THEN
        nb1=NBJ(1,ne)
        np1=NPNE(NNIP(i1,i2,i3),nb1,ne)
        np2=NPNE(NNIP(i1+1,i2,i3),nb1,ne)
        DO nogrno = 1,NTGRNO
C KAT 3Nov99 Updating for dynamic groups
          IF(INLIST(np1,%VAL(LIGRNO_PTR(nogrno)),
     '      NLIGRNO(nogrno),no_grlist)) THEN
            IF(INLIST(np2,%VAL(LIGRNO_PTR(nogrno)),
     '        NLIGRNO(nogrno),no_grlist)) THEN
              ILIST_PTR=0
              CALL ALLOCATE_MEMORY(NLIGRNO(nogrno)+1,0,INTTYPE,
     '          ILIST_PTR,MEM_INIT,ERROR,*9999)
              CALL ILIST_COPY(NLIGRNO(nogrno),
     '          %VAL(LIGRNO_PTR(nogrno)),%VAL(ILIST_PTR))
              NLIGRNO(nogrno)=NLIGRNO(nogrno)+1
              CALL ILIST_COPY(1,
     '          NODE,%VAL(ILISTLOC(%VAL(ILIST_PTR),NLIGRNO(nogrno))))
              CALL FREE_MEMORY(LIGRNO_PTR(nogrno),ERROR,*9999)
              LIGRNO_PTR(nogrno)=ILIST_PTR
            ENDIF
          ENDIF
C          DO no_grlist = 1,0!LIGRNO(0,nogrno)
C            IF(LIGRNO(no_grlist,nogrno).EQ.np1) np1_IN_GROUP=.TRUE.
C            IF(LIGRNO(no_grlist,nogrno).EQ.np2) np2_IN_GROUP=.TRUE.
C          ENDDO !no_grlist
C          IF(np1_IN_GROUP.AND.np2_IN_GROUP) THEN
C            LIGRNO(0,nogrno)=LIGRNO(0,nogrno)+1
C            LIGRNO(LIGRNO(0,nogrno),nogrno)=NODE
C          ENDIF
        ENDDO !nogrno
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'(/'' New node properites:'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NODE='',I5)') NODE
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          DO njj1=1,3
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj1=NJ_LOC(njj1,njj2,nr)
              nb1=NBJ(nj1,ne)
              nb2=NBJ(nj1,nenew)
              WRITE(OP_STRING,'('' njj1='',I1,'', njj2='',I1,'
     '          //''', nj='',I2,'', nb1='',I2,'', nb2='',I2)')
     '          njj1,njj2,nj1,nb1,nb2
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' NVJP(nj,NODE)='',I2)')
     '          NVJP(nj1,NODE)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              DO nv=1,NVJP(nj1,NODE)
                WRITE(OP_STRING,'('' nv='',I2,'', XP: '',8D12.3)')
     '            nv,(XP(nk,nv,nj1,NODE),nk=1,NKJ(nj1,NODE))
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDDO !nv
              WRITE(OP_STRING,'('' NVJE(1..,nb1,nj,ne)   : '',8I3)')
     '          (NVJE(nn,nb1,nj1,ne),nn=1,NNT(nb1))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' NVJE(1..,nb2,nj,nenew): '',8I3)')
     '          (NVJE(nn,nb2,nj1,nenew),nn=1,NNT(nb2))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO !njj2
          ENDDO !njj1
CC$        call mp_unsetlock()
        ENDIF
        EXISTS=.TRUE.
      ENDIF

      CALL EXITS('REFINE_SETNODE')
      RETURN
 9999 CALL ERRORS('REFINE_SETNODE',ERROR)
      CALL EXITS('REFINE_SETNODE')
      RETURN 1
      END


