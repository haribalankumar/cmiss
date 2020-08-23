      SUBROUTINE GETEQVNONY(IDO,INP,NBH,NBJ,NENP,NKB,NKHE,NNB,
     '  NONY,NPNE,NPNY,nr,NVHE,NVHP,NXI,ny1,ny2,
     '  NYNP,NYNY,CYNY,RATIO,FIX,*)

C#### Subroutine: GETEQVNONY
C###  Description:
C###    GETEQVNONY searches for a corresponding col number ny2 to ny1 for
C###    which a solution degree of freedom (no) has already been set up.
C###    RATIO is set to 1 or -1 according to whether the corresponding col
C###    number is equal or opposite (derivative direction).  If there is
C###    no corresponding col number that has already been set up then ny2
C###    is set to ny1.  If there is no corresponding solution degree of
C###    freedom then ny2 is set to 0.

      IMPLICIT NONE
      CHARACTER*(*) ROUTINENAME
      PARAMETER(ROUTINENAME='GETEQVNONY')
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
!     Parameter List
      INTEGER IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NENP(NPM,0:NEPM),NKB(2,2,2,NNM,NBFM),
     '  NKHE(NKM,NNM,NHM,NEM),NNB(4,4,4,NBFM),NONY(0:NOYM,NYM,NRCM),
     '  NPNE(NNM,NBFM,NEM),NPNY(0:6,NYM,0:NRCM),nr,
     '  NVHE(NNM,NBFM,NHM,NEM),
     '  NVHP(NHM,NPM,NCM),NXI(-NIM:NIM,0:NEIM,0:NEM),ny1,ny2,
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM),NYNY(0:NYYM,NYM)
      REAL*8 CYNY(0:NYYM,NYM),RATIO
      LOGICAL FIX(NYM,NIYFIXM)
!     Local Variables
      INTEGER i,IDNNU(11),IDOI(3),IDONK(8),INPI1(3),INPI2(3),
     '  INPSTART,INPFSTART,INPT(2),IRATIO,IRATIOF,j,nb1,nb2,nc,ND1,ne1,
     '  ne2,nep,nep1,NEPT,nh,ni,ni1,ni2,nif,NIM3,NIX1(2),NIX2(2),
     '  nk1,nk2,nk_e,NKTB,nn1,nn2,nnn1,nnn2,NNTB,NOYT,np,nv1,nv2,NVT,ny
      LOGICAL COLLAPSED,DONE,FOUND
      CHARACTER ERROR*100

      DATA IDNNU/0,1,0,2,0,3,3,0,2,1,0/ !deriv or cross dirn of 2 ders
      DATA IDONK/0,1,1,2,1,2,2,3/ !deriv order

      CALL ENTERS(ROUTINENAME,*9999)

C     Set default output if no equivalent no ny has been set up.
      ny2=ny1
      RATIO=1

      IF(NPNY(0,ny1,0).EQ.1) THEN !ny is node based
        nh=NPNY(3,ny1,0)
        np=NPNY(4,ny1,0)
        nc=NPNY(5,ny1,0)
C       There may an equivalent dof if there is >1 version at this node.
        NVT=NVHP(nh,np,nc)
        IF(NVT.GT.1) THEN
          nv1=NPNY(2,ny1,0)
          nk1=NPNY(1,ny1,0)
          ND1=IDONK(nk1) !derivative order
          DONE=.FALSE.
          IF((ND1.EQ.0).AND.(ITYP21(nr).LE.2)) THEN !not a derivative
C           Search all other versions to see if mapping is already done.
            nv2=0
            DO WHILE(nv2.LT.NVT.AND..NOT.DONE)
              nv2=nv2+1
              IF(nv2.NE.nv1) THEN
                ny=NYNP(nk1,nv2,nh,np,0,nc)
                IF(FIX(ny,1)) THEN
                  DONE=.TRUE.
                ELSE
                  NOYT=NONY(0,ny,2)
                  DONE=NOYT.NE.0 !assuming nrc=1 is done too
                ENDIF
                IF(KTYP11.EQ.2) THEN ! faces
                  IF(.NOT.FIX(ny,NIYFIXM)) THEN
                    IF(FIX(ny,1)) THEN
                       DONE=.TRUE.
                     ELSE
                       NOYT=NONY(0,ny,2)
                       DONE=NOYT.NE.0 !assuming nrc=1 is done too
                     ENDIF
                  ENDIF
                ENDIF
              ENDIF !nv2.NE.nv1
            ENDDO !nv2
          ELSE IF((ND1.LT.3).AND.(ITYP21(nr).EQ.1)) THEN !not 3rd deriv
C           Search for an element using the dof ny1
            NEPT=NENP(np,0)
            FOUND=.FALSE.
            nep1=0
            DO WHILE(nep1.LT.NEPT.AND..NOT.FOUND)
              nep1=nep1+1
              ne1=NENP(np,nep1)
C             Search test element for nv1 of np
              IF(KTYP11.EQ.2) THEN !faces
                nb1=NBJ(1,ne1) !!!!just use the first geom basis for now
              ELSE
              nb1=NBH(nh,nc,ne1)
              ENDIF
              NNTB=NNT(nb1)
              nn1=0
              DO WHILE(nn1.LT.NNTB.AND..NOT.FOUND)
                nn1=nn1+1
                IF(NPNE(nn1,nb1,ne1).EQ.np
     '            .AND.NVHE(nn1,nb1,nh,ne1).EQ.nv1) THEN
C                 Check that the same global nk is used
                  NKTB=NKT(nn1,nb1)
                  nk_e=0
                  DO WHILE(nk_e.LT.NKTB.AND..NOT.FOUND)
                    nk_e=nk_e+1
                    IF(NKHE(nk_e,nn1,nh,ne1).EQ.nk1) THEN
                      FOUND=.TRUE.
                      ni1=IDNNU(IDO(nk_e,nn1,0,nb1))
                      NIX1(1)=1+MOD(ni1,3)
                      NIX1(2)=6-ni1-NIX1(1)
                    ENDIF !nk1
                  ENDDO !nk_e
                ENDIF !np,nv1
              ENDDO !nn1
            ENDDO !nep1
            IF(.NOT.FOUND) THEN
C             dof was not found in an element so it is not
C             used and no mapping is needed.
              ny2=0
              WRITE(OP_STRING,
     '          '('' >>WARNING: Degree of freedom not used'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ELSE IF((ND1.EQ.1).AND.(ITYP21(nr).EQ.1)) THEN !1st deriv
C             Node position indices in first element found.
C             There is an assumption here that yp1 is only used on one
C             line of this element.
              NIM3=MIN(NIM,3)
              DO i=1,2
                ni=NIX1(i)
                IF(ni.LE.NIM) THEN
                  INPI1(ni)=INP(nn1,ni,nb1)
                ELSE
                  INPI1(ni)=1
                  INPI2(ni)=1
                ENDIF
              ENDDO
C             If the line corresponding to this derivative is collapsed
C             then fix the dof.
              INPI1(ni1)=3-INP(nn1,ni1,nb1)
              nn2=NNB(INPI1(1),INPI1(2),INPI1(3),nb1)
              IF(NPNE(nn2,nb1,ne1).EQ.np) THEN
                DONE=.TRUE.
                ny=0
              ENDIF !collapsed
C             Find number of nodes on this line
              INPI1(ni1)=2
              DO WHILE(NNB(INPI1(1),INPI1(2),INPI1(3),nb1).NE.0
     '          .AND.INPI1(ni1).LE.4)
                INPI1(ni1)=INPI1(ni1)+1
              ENDDO
              INPT(1)=INPI1(ni1)-1
C             Search for elements using different versions of np.
              NEPT=NENP(np,0)
              nep=0
              DO WHILE(nep.LT.NEPT.AND..NOT.DONE)
                nep=nep+1
C               Don't exclude nep1 as a collapsed face will have
C               corresponding dofs in the same element.  But perhaps
C               this could be optimized for when nep==nep1.
                ne2=NENP(np,nep)
C               Search test element for np
                nb2=NBH(nh,nc,ne2)
                NNTB=NNT(nb2)
                nn2=0
                DO WHILE(nn2.LT.NNTB.AND..NOT.DONE)
                  nn2=nn2+1
                  IF(NPNE(nn2,nb2,ne2).EQ.np) THEN !same node
                    nv2=NVHE(nn2,nb2,nh,ne2)
                    IF(nv2.NE.nv1) THEN !diff version
C                     Check all (inconsistent) xi directions
                      ni2=0
                      DO WHILE(ni2.LT.NIM3.AND..NOT.DONE)
                        ni2=ni2+1
                        NIX2(1)=1+MOD(ni2,3)
                        NIX2(2)=6-ni2-NIX2(1)
                        DO i=1,2
                          ni=NIX2(i)
                          IF(ni.LE.NIM) INPI2(ni)=INP(nn2,ni,nb2)
                        ENDDO
C                       Check + and - directions
                        IRATIO=3
                        INPSTART=-1-INPT(1)
                        DO WHILE(IRATIO.GT.0.AND..NOT.DONE)
                          IRATIO=IRATIO-2
                          INPSTART=INPSTART+INPT(1)+1
C                         Check if nodes in deriv dirn are same.
                          FOUND=.TRUE.
                          INPI1(ni1)=0
                          INPI2(ni2)=INPSTART
                          DO WHILE(INPI1(ni1).LT.INPT(1).AND.FOUND)
                            INPI1(ni1)=INPI1(ni1)+1
                            INPI2(ni2)=INPI2(ni2)+IRATIO
                            nnn1=NNB(INPI1(1),INPI1(2),INPI1(3),nb1)
                            nnn2=NNB(INPI2(1),INPI2(2),INPI2(3),nb2)
                            IF(nnn2.EQ.0) THEN !nnn1!=0
                              FOUND=.FALSE.
                            ELSEIF(NPNE(nnn2,nb2,ne2)
     '                          .NE.NPNE(nnn1,nb1,ne1)) THEN
                              FOUND=.FALSE.
                            ENDIF
                          ENDDO !inpi1(ni1)
                          IF(FOUND.AND.INPT(1).LT.4) THEN
C                           Check that there are the same number of
C                           nodes in each line.  A shorter line will
C                           have been detected by zero nnn.  Here we
C                           check that we don't have a longer line.
                            INPI2(ni2)=INPT(1)+1
                            FOUND=
     '                        NNB(INPI2(1),INPI2(2),INPI2(3),nb2).EQ.0
                          ENDIF
                          IF(FOUND) THEN
                            RATIO=DBLE(IRATIO)
                            IDOI(ni2)=2 !1st deriv
                            IDOI(NIX2(1))=1 !no deriv
                            IDOI(NIX2(2))=1 !no deriv
                            nk_e=NKB(IDOI(1),IDOI(2),IDOI(3),nn2,nb2)
                            IF(nk_e.NE.0) THEN
                              nk2=NKHE(nk_e,nn2,nh,ne2)
                              ny=NYNP(nk2,nv2,nh,np,0,nc)
                              IF(FIX(ny,1)) THEN
                                DONE=.TRUE.
                              ELSE
                                NOYT=NONY(0,ny,2)
                                DONE=NOYT.NE.0 !assuming nrc=1 is done too
                              ENDIF !FIX
                            ENDIF !nk_e
                          ENDIF !FOUND
                        ENDDO !+-
                      ENDDO !ni2
                    ENDIF !nv2.NE.nv1
                  ENDIF !np
                ENDDO !nn2
              ENDDO !nep
            ELSE IF(ITYP21(nr).EQ.1) THEN!2nd deriv
C             There is an assumption here that yp1 is only used on one
C             face of this element.
              INPI1(ni1)=INP(nn1,ni1,nb1)
              DO i=1,2
                ni=NIX1(i)
                nif=NIX1(3-i)
C               Find the number of nodes in this direction on the face.
C               Assumes that the largest number of nodes in this
C               direction can be found by searching the first line in
C               NNB.  Multidirectional mapping will probably only work on
C               square elements.
                INPI1(nif)=1! first line
                INPI1(ni)=2
                DO WHILE(NNB(INPI1(1),INPI1(2),INPI1(3),nb1).NE.0
     '            .AND.INPI1(ni).LE.4)
                  INPI1(ni)=INPI1(ni)+1
                ENDDO
                INPT(i)=INPI1(ni)-1
C               If the face corresponding to this derivative is
C               collapsed then fix the dof.
                COLLAPSED=.TRUE.
                INPI1(ni)=0
                DO WHILE(INPI1(ni).LT.INPT(i).AND.COLLAPSED)
                  INPI1(ni)=INPI1(ni)+1
                  INPI1(nif)=1
                  nnn1=NNB(INPI1(1),INPI1(2),INPI1(3),nb1)
                  INPI1(nif)=2
                  nnn2=NNB(INPI1(1),INPI1(2),INPI1(3),nb1)
                  IF(nnn2.NE.0.AND.nnn1.NE.0) THEN
                    IF(NPNE(nnn2,nb1,ne1).NE.NPNE(nnn1,nb1,ne1)) THEN
                      COLLAPSED=.FALSE.
                    ENDIF !NPNE
                  ENDIF !nnn2
                ENDDO !inpi1(ni)
                IF(COLLAPSED) THEN
                  DONE=.TRUE.
                  ny=0
                ENDIF !collapsed
              ENDDO !i
C             Find adjacent element sharing the cross deriv
              IF(INPI1(ni1).EQ.1) THEN
                ne2=NXI(-ni1,1,ne1)
              ELSE
                ne2=NXI(ni1,1,ne1)
              ENDIF
              IF(ne2.NE.0.AND..NOT.DONE) THEN
C               Search test element for np
                nb2=NBH(nh,nc,ne2)
                NNTB=NNT(nb2)
                nn2=0
                DO WHILE(nn2.LT.NNTB.AND..NOT.DONE)
                  nn2=nn2+1
                  IF(NPNE(nn2,nb2,ne2).EQ.np) THEN !same node
                    nv2=NVHE(nn2,nb2,nh,ne2)
                    IF(nv2.NE.nv1) THEN !diff version
C                     Check all (inconsistent) face-normal xi directions
                      ni2=0
                      DO WHILE(ni2.LT.3.AND..NOT.DONE)
                        ni2=ni2+1
                        INPI2(ni2)=INP(nn2,ni2,nb2)
C                       Check both (inconsistent) face xi directions
                        j=0
                        DO WHILE(j.LT.2.AND..NOT.DONE)
                          NIX2(1)=1+MOD(ni2+j,3)
                          NIX2(2)=6-ni2-NIX2(1)
                          j=j+1
C                         Check that there are the same number of nodes
C                         in each direction.  Less will be detected
C                         by inconsistent nnn.  Here we check that we
C                         don't have more.
                          FOUND=.TRUE.
                          DO i=1,2
                            IF(INPT(i).LT.4) THEN
                              ni=NIX2(i)
                              nif=NIX2(3-i)
                              INPI2(nif)=1 ! first line
                              INPI2(ni)=INPT(i)+1
                              IF(NNB(INPI2(1),INPI2(2),INPI2(3),nb2)
     '                          .NE.0) FOUND=.FALSE.
                            ENDIF
                          ENDDO !i
                          IF(FOUND) THEN
C                           Check + and - directions
                            IRATIO=3
                            INPSTART=-1-INPT(1)
                            DO WHILE(IRATIO.GT.0.AND..NOT.DONE)
                              IRATIO=IRATIO-2
                              INPSTART=INPSTART+INPT(1)+1
C                             Check that local node is in the same
C                             position on the face.  This needs to be
C                             done in addition to checking nodes in
C                             face are same as the global node may
C                             appear more than once in the element
C                             with different versions.
                              IF(INPSTART+IRATIO*INP(nn2,NIX2(1),nb2)
     '                          .EQ.INP(nn1,NIX1(1),nb1)) THEN
                                IRATIOF=3
                                INPFSTART=-1-INPT(2)
                                DO WHILE(IRATIOF.GT.0.AND..NOT.DONE)
                                  IRATIOF=IRATIOF-2
                                  INPFSTART=INPFSTART+INPT(2)+1
                                  IF(INPFSTART+
     '                              IRATIOF*INP(nn2,NIX2(2),nb2)
     '                              .EQ.INP(nn1,NIX1(2),nb1)) THEN
C                     Check if nodes in face are same.
                      INPI1(NIX1(1))=0
                      INPI2(NIX2(1))=INPSTART
                      FOUND=.TRUE.
                      DO WHILE(INPI1(NIX1(1)).LT.INPT(1).AND.FOUND)
                        INPI1(NIX1(1))=INPI1(NIX1(1))+1
                        INPI2(NIX2(1))=INPI2(NIX2(1))+IRATIO
                        INPI1(NIX1(2))=0
                        INPI2(NIX2(2))=INPFSTART
                        DO WHILE(INPI1(NIX1(2)).LT.INPT(2).AND.FOUND)
                          INPI1(NIX1(2))=INPI1(NIX1(2))+1
                          INPI2(NIX2(2))=INPI2(NIX2(2))+IRATIOF
                          nnn1=NNB(INPI1(1),INPI1(2),INPI1(3),nb1)
                          nnn2=NNB(INPI2(1),INPI2(2),INPI2(3),nb2)
C                         Can't use NPNE with nnn2 or nnn1 0
                          IF(nnn2.EQ.0.OR.nnn1.EQ.0) THEN
                            IF(nnn2.NE.nnn1) FOUND=.FALSE.
                          ELSEIF(NPNE(nnn2,nb2,ne2)
     '                        .NE.NPNE(nnn1,nb1,ne1)) THEN
                            FOUND=.FALSE.
                          ENDIF
                        ENDDO !INPI(NIX(2))
                      ENDDO !INPI(NIX(1))
                      IF(FOUND) THEN
                        RATIO=DBLE(IRATIO*IRATIOF)
                        IDOI(ni2)=1 !no deriv
                        IDOI(NIX2(1))=2 !1st deriv
                        IDOI(NIX2(2))=2 !1st deriv
                        nk_e=NKB(IDOI(1),IDOI(2),IDOI(3),nn2,nb2)
                        IF(nk_e.NE.0) THEN
                          nk2=NKHE(nk_e,nn2,nh,ne2)
                          ny=NYNP(nk2,nv2,nh,np,0,nc)
                          IF(FIX(ny,1)) THEN
                            DONE=.TRUE.
                          ELSE
                            NOYT=NONY(0,ny,2)
                            DONE=NOYT.NE.0 !assuming nrc=1 is done too
                          ENDIF !FIX
                        ENDIF !nk_e
                      ENDIF !FOUND
                                  ENDIF !INP NIX(2)
                                ENDDO !IRATIOF
                              ENDIF !INP NIX(1)
                            ENDDO !IRATIO
                          ENDIF !same number of nodes in each direction
                        ENDDO !j (NIX2)
                      ENDDO !ni2
                    ENDIF !nv2.NE.nv1
                  ENDIF !np
                ENDDO !nn2
              ENDIF !ne2.NE.0
            ENDIF !ND1
          ENDIF !ND1
          IF(DONE) THEN !corresponding solution dof is already mapped
            ny2=ny
!        print *,'GETEQVNONY: NYNY(0,',ny1,')=',NYNY(0,ny1),' ny2=',ny2

            IF(DOP) THEN
              WRITE(OP_STRING,'('' GETEQVNONY:'',I10,'' ->'',I10)')
     '          NY1,NY2
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF !DONE
        ENDIF !NVT
      ENDIF !node based

C new CS 29/9/2001 Explicit control for arbitrary ny->ny mapping
      ! This is done after the automatic mapping rather than instead
      ! of it to allow the addition of extra non-standard mappings
      ! after the automatic mapping code has created the standard ones.
      ! Note that the automatic mapping does not correctly determine
      ! mappings for all cases, namely collapsed elements that do
      ! not occur in an apex configuration.
      IF(ITYP21(nr).GE.2) THEN
        !print *,'explicit mapping: NYNY(0,',ny1,')=',NYNY(0,ny1)
        IF(NYNY(0,ny1).NE.0) THEN
          ny2=NYNY(1,ny1)
          RATIO=CYNY(1,ny1)
         !print *,'explicit mapping: NYNY(1,',ny1,')=',NYNY(1,ny1)
          IF(DOP) THEN
              WRITE(OP_STRING,'('' Expl Mapping:'',I10,'' ->'',I10)')
     '          NY1,NY2
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
      ENDIF

      CALL EXITS(ROUTINENAME)
      RETURN

 9999 CALL ERRORS(ROUTINENAME,ERROR)
      CALL EXITS(ROUTINENAME)
      RETURN 1
      END


