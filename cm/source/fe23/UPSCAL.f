      SUBROUTINE UPSCAL(IBT,IDO,INP,NBJ,NEELEM,NEL,NENP,
     '  NKJE,NKJ,NLL,NLLINE,NNL,NPL,NPNE,NPNODE,
     '  NRE,NRLIST,NSB,NVJE,NVJL,NVJP,DL,DLL,SE,XP,STRING,ERROR,*)

C#### Subroutine: UPSCAL
C###  Description:
C###    UPSCAL updates scale factors for a basis.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'iter00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NEL(0:NELM,NLM),
     '  NENP(NPM,0:NEPM,0:NRM),NKJE(NKM,NNM,NJM,NEM),NKJ(NJM,NPM),
     '  NLL(12,NEM),NLLINE(0:NL_R_M,0:NRM),NNL(0:4,12,NBFM),
     '  NPL(5,0:3,NLM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NRE(NEM),NRLIST(0:NRM),NSB(NKM,NNM,NBFM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJL(4,NJM,NLM),NVJP(NJM,NPM)
      REAL*8 DL(3,NLM),DLL(3,NLM),SE(NSM,NBFM,NEM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,JEST,nb,NBB1,NBB2,ne,nep,NEPT,nj,njj,
     '  NJL_TYPE,nk,nk_e,nk1,NKLIST(3),NKMAX,NKTB,NKTOT,nl,N3CO,
     '  nn,NNTB,no_np,no_nr,np,nr,ns,NUMBER,nv,NVTOT
      !SMAR009 22/12/98 ,nk2
      REAL*8 SCALE,SL2,SM2,SUM,XD(3)
      LOGICAL ALL_REGIONS,CBBREV,FIELD,FOUND,GEOMETRY,
     '  ITERATE,NORMALISE,SOLUTION,UNIT
      DATA NKLIST/2,3,5/

      CALL ENTERS('UPSCAL',*9999)

      IF(CO(noco+1).EQ.'?') THEN

C---------------------------------------------------------------------

C#### Command: FEM update scale_factors
C###  Parameter:    <region (all/#s)[1]>
C###    Specify the region numbers to update.
C###  Parameter:    <estimate>
c###   Specify wiether  initial estimates of scale factors are made first
C###  Parameter:    <normalise/unit>
c###   Specify wiether  nodal geometry/field/solution derivatives
C###   are normalised first
C###  Parameter:    <geometry/field/solution[geometry]>
C###  Specify wiether geometry,field or solution is to be updated
C###  Description:
C###    Recalculates the line lengths and updates the scale factors
C###    based on the scaling model and the new line lengths.
C###    If 'estimate' is specified then initial estimates of scale
C###    factors are made first.  Otherwise the present values are used
C###    as estimates.  If 'normalise' is specified then the nodal
C###    geometry/field/solution derivatives are normalised first.  Fibre
C###    derivative will also be updated if 'geometry' is specified and
C###    scale factors are node based.  If 'unit' is specified then node
C###    based scale factors are set to 1 and the nodal derivatives are
C###    updated attempting to retain the same interpolated values.  The
C###    <geometry/field/solution> option applies only to 'normalise'.
C###    For 'unit', geometry, fibre, and field values will be updated
C###    but solution value updates are not yet implemented.

        CALL STRING_TRIM(STRING,IBEG,IEND)
        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<region (all/#s)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<estimate>'
        OP_STRING(4)=BLANK(1:15)//'<normalise/unit>'
        OP_STRING(5)=BLANK(1:15)//'<geometry/field/solution[geometry]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update scale_factors iterate
C###  Parameter:      <basis (#/all)[all]>
C###    Specify the basis function numbers
C###  Description:
C###    Recalculates the line lengths and updates the scale factors
C###    based on the scaling model and the new line lengths

        CALL STRING_TRIM(STRING,IBEG,IEND)
        OP_STRING(1)=STRING(1:IEND)//' iterate'
        OP_STRING(2)=BLANK(1:15)//'<basis (#/all)[all]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe23','doc','UPSCAL',ERROR,*9999)
      ELSE
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        IF(CBBREV(CO,'ITERATE',1,noco+1,NTCO,N3CO)) THEN
          ITERATE=.TRUE.
          IF(CBBREV(CO,'BASIS',1,noco+1,NTCO,N3CO)) THEN
            nb=IFROMC(CO(N3CO+1))
            NBB1=nb
            NBB2=nb
          ELSE
            NBB1=1
            NBB2=NBFT
          ENDIF
        ELSE
          ITERATE=.FALSE.
          NORMALISE=CBBREV(CO,'NORMALISE',2,noco+1,NTCO,N3CO)
          UNIT=CBBREV(CO,'UNIT',2,noco+1,NTCO,N3CO)
          CALL ASSERT(.NOT.(UNIT.AND.NORMALISE),
     '      '>>Can''t normalise with unit scale factors',ERROR,*9999)
          IF(CBBREV(CO,'ESTIMATE',2,noco+1,NTCO,N3CO)) THEN
            JEST=0
          ELSE
            JEST=1
          ENDIF
          FIELD=.FALSE.
          SOLUTION=.FALSE.
          GEOMETRY=.FALSE.
          IF(CBBREV(CO,'FIELD',2,noco+1,NTCO,N3CO)) THEN
            FIELD=.TRUE.
          ELSE IF(CBBREV(CO,'SOLUTION',2,noco+1,NTCO,N3CO)) THEN
            SOLUTION=.TRUE.
          ELSE
            GEOMETRY=.TRUE.
          ENDIF
        ENDIF

        IF(ITERATE) THEN
          DO nl=1,NLT
            DLL(3,nl)=DL(3,nl)
          ENDDO
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
            WRITE(OP_STRING,'('' DLL(3,nl): '',6D13.5,/(12X,'
     '        //'6D13.5))') (DLL(3,nl),nl=1,NLT)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
          ENDIF
          CALL LINCAL(IBT,IDO,INP,0,NBJ,NEELEM,NEL,NENP,
     '      NKJE,NLL,NLLINE,NNL,NPL,NPNE,NPNODE,NRE,NVJE,NVJL,
     '      DL,SE,XP,ERROR,*9999)
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
            WRITE(OP_STRING,'('' DL(3,nl):  '',6D13.5,/(12X,'
     '        //'6D13.5))') (DL(3,nl),nl=1,NLT)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(/'' ITER_ERR_OLD(2)='',D12.5)')
     '        ITER_ERR_OLD(2)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
          ENDIF

          IF(KTYP20.EQ.1) THEN
          ELSE IF(KTYP20.EQ.2) THEN
            IF(KTYP21.EQ.1) THEN !Data fitting by optimistion
              ITER_ERR_OLD(2)=ITER_ERR(2)
              ITER_ERR(2)=0.0d0
              DO nl=1,NLT
                IF(DABS(DL(3,nl)-DLL(3,nl)).GT.ITER_ERR(2))
     '            ITER_ERR(2)=DABS(DL(3,nl)-DLL(3,nl))
              ENDDO
              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
                WRITE(OP_STRING,'('' ITER_ERR(2)    ='',D12.5)')
     '            ITER_ERR(2)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
              ENDIF
            ENDIF
          ENDIF

          DO nb=NBB1,NBB2
            CALL DLSE(IBT,IDO,nb,NEELEM,NLL,NNL,NPL,
     '        DL,SE,ERROR,*9999)
          ENDDO

        ELSE
          IF(NORMALISE) THEN
C           update normalised arc-length derivs and scale factors
            CALL ASSERT(.NOT.SOLUTION,'>> Solution update is not '
     '        //'implemented yet',ERROR,*9999)
            IF(GEOMETRY) NJL_TYPE=NJL_GEOM
            IF(FIELD) NJL_TYPE=NJL_FIEL

            DO no_nr=1,NRLIST(0)
              nr=NRLIST(no_nr)
              DO no_np=1,NPNODE(0,nr)
                np=NPNODE(no_np,nr)
                CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
C               find number of first derivatives at that node
                nj=NJ_LOC(NJL_TYPE,1,nr)
                NKMAX=NKJ(nj,np)
                NVTOT=NVJP(nj,np)
                DO njj=2,NJ_LOC(NJL_TYPE,0,nr)
                  nj=NJ_LOC(NJL_TYPE,njj,nr)
                  IF(NKJ(nj,np).NE.NKMAX.AND.NKMAX.NE.0) THEN
                    WRITE(OP_STRING,
     '                '('' Differing numbers of derivatives at node '''
     '                //',I3)') np
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    NKMAX=0
                  ENDIF
                  IF(NVJP(nj,np).NE.NVTOT.AND.NVTOT.NE.0) THEN
                    WRITE(OP_STRING,
     '                '('' Differing numbers of versions at node '''
     '                //',I3)') np
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    NVTOT=0
                  ENDIF
                ENDDO
                IF(NKMAX.EQ.2) THEN
                  NKTOT=1
                ELSE IF(NKMAX.EQ.4) THEN
                  NKTOT=2
                ELSE IF(NKMAX.EQ.8) THEN
                  NKTOT=3
                ELSE
                  NKTOT=0
                ENDIF
                DO nv=1,NVTOT
                  IF(ITYP10(1).EQ.4) THEN !prolate spheroidal
                    nj=NJ_LOC(NJL_TYPE,1,nr) !lambda
                    SL2=DSINH(XP(1,nv,nj,np))**2
                    nj=NJ_LOC(NJL_TYPE,2,nr) !mu
                    SM2=DSIN(XP(1,nv,nj,np))**2
                  ENDIF
                  DO nk1=1,NKTOT !loop over 1st derivs
                    nk=NKLIST(nk1)
                    DO njj=1,NJ_LOC(NJL_TYPE,0,nr)
                      nj=NJ_LOC(NJL_TYPE,njj,nr)
                      XD(njj)=XP(nk,nv,nj,np)
                    ENDDO
                    DO njj=NJ_LOC(NJL_TYPE,0,nr)+1,3
                      XD(njj)=0.0d0
                    ENDDO
                    IF(ITYP10(1).EQ.1) THEN !rectangular Cartesian
                      SCALE=XD(1)**2+XD(2)**2+XD(3)**2
                    ELSE IF(ITYP10(1).EQ.4) THEN !prolate spheroidal
                      SCALE=FOCUS*FOCUS*(
     '                  (SL2+SM2)*(XD(1)**2+XD(2)**2)+SL2*SM2*XD(3)**2)
                    ELSE
                      ERROR='>> Not implemented for these coords'
                      GO TO 9999
                    ENDIF
                    IF(SCALE.GT.ZERO_TOL) THEN
                      SCALE=1.d0/DSQRT(SCALE)
C                     Normalise the nodal derivatives
                      DO njj=1,NJ_LOC(NJL_TYPE,0,nr)
                        nj=NJ_LOC(NJL_TYPE,njj,nr)
                        XP(nk,nv,nj,np)=XP(nk,nv,nj,np)*SCALE
C                       adjust cross-derivatives
                        IF(NKTOT.GE.2) THEN
                          IF(nk1.EQ.1.OR.nk1.EQ.2) THEN
                            XP(4,nv,nj,np)=XP(4,nv,nj,np)*SCALE
                          ENDIF
                          IF(NKTOT.EQ.3) THEN
                            IF(nk1.EQ.1.OR.nk1.EQ.3) THEN
                              XP(6,nv,nj,np)=XP(6,nv,nj,np)*SCALE
                            ENDIF
                            IF(nk1.EQ.2.OR.nk1.EQ.3) THEN
                              XP(7,nv,nj,np)=XP(7,nv,nj,np)*SCALE
                            ENDIF
                            XP(8,nv,nj,np)=XP(8,nv,nj,np)*SCALE
                          ENDIF
                        ENDIF
                      ENDDO ! njj
                    ENDIF !>0
                  ENDDO ! nk1
                ENDDO !nv
                IF(GEOMETRY) THEN !update fibres also
C                 Multiply by present scale factors (divide by new later)
                  DO njj=1,NJ_LOC(NJL_FIBR,0,nr)
                    nj=NJ_LOC(NJL_FIBR,njj,nr)
                    DO nv=1,NVJP(nj,np)
                      DO nk=2,NKJ(nj,np) !skip nodal value
                        FOUND=.FALSE.
C                       Search for an element using the dof.
                        NEPT=NENP(np,0,nr)
                        nep=0
                        DO WHILE(nep.LT.NEPT.AND..NOT.FOUND)
                          nep=nep+1
                          ne=NENP(np,nep,nr)
C                         Search test element for version nv of np
                          nb=NBJ(nj,ne)
                          NNTB=NNT(nb)
                          nn=0
                          DO WHILE(nn.LT.NNTB.AND..NOT.FOUND)
                            nn=nn+1
                            IF(NPNE(nn,nb,ne).EQ.np
     '                        .AND.NVJE(nn,nb,nj,ne).EQ.nv) THEN
C                             Check that the same global nk is used
                              NKTB=NKT(nn,nb)
                              nk_e=0
                              DO WHILE(nk_e.LT.NKTB.AND..NOT.FOUND)
                                nk_e=nk_e+1
                                IF(NKJE(nk_e,nn,nj,ne).EQ.nk) THEN
                                  FOUND=.TRUE.
                                ENDIF !nk
                              ENDDO !nk_e
                            ENDIF !np,nv
                          ENDDO !nn1
                        ENDDO !nep1
                        IF(.NOT.FOUND) THEN
                          IEND=0
                          CALL APPENDC(IEND,
     '                      ' No element found for fibre angle ',
     '                      OP_STRING(1))
                          CALL APPENDI(IEND,njj,OP_STRING(1))
                          CALL APPENDC(IEND,' node ',OP_STRING(1))
                          CALL APPENDI(IEND,np,OP_STRING(1))
                          CALL APPENDC(IEND,' version ',OP_STRING(1))
                          CALL APPENDI(IEND,nv,OP_STRING(1))
                          CALL APPENDC(IEND,' derivative ',OP_STRING(1))
                          CALL APPENDI(IEND,nk,OP_STRING(1))
C                          WRITE(OP_STRING,
C     '                      '('' No element found for node '',I3,'//
C     '                      ''', version '',I3,'', derivative '',I1)')
C     '                      np,nv,nk
                          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                        ELSE
                          IF(NBI(nb).EQ.3.
     '                      OR.NBI(nb).EQ.6.OR.NBI(nb).EQ.7) THEN
C                           Nodal scale factor
                            ns=NSB(nk_e,nn,nb)
                            XP(nk,nv,nj,np)=XP(nk,nv,nj,np)*SE(ns,nb,ne)
                          ENDIF
                        ENDIF !found scale factor
                      ENDDO !nk
                    ENDDO !nv
                  ENDDO !njj
                ENDIF !GEOMETRY
              ENDDO !no_np

              CALL ASSERT(.NOT.FIELD,'>> Nodal derivatives have been'
     '          //'normalised, but arc-length derivatives and '
     '          //'scale factors have not been updated',ERROR,*9999)

C KAT 17Sep98: olds
C            DO no_nl=1,NLLINE(0,nr) !initialise DLL
C              nl=NLLINE(no_nl,nr)
C              DLL(3,nl)=DL(3,nl)
C            ENDDO !nl
C
C            jder=0
C            jest=0
C            jsca=0
C            iter=0
C            CONVERGED=.FALSE.
C            DO WHILE ((.NOT.CONVERGED).AND.(iter.LT.100))
C              MAX_DIFF=0.0d0
C              iter=iter+1
C              IF(DOP) THEN
CC$              call mp_setlock()
C                WRITE(OP_STRING,'('' Iteration '',I2)') iter
C                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
C              ENDIF
C              DO no_nl=1,NLLINE(0,nr) !loop over global lines
C                nl=NLLINE(no_nl,nr)
CC ***           Calculate DL array derivatives from arclength
C                IF(NBI(NBJ(1,NEL(1,nl))).EQ.1) THEN
CC                 Set DL array derivatives to unity
C                  DL(1,nl)=1.0D0
C                  DL(2,nl)=1.0D0
C                ELSE IF(NBI(NBJ(1,NEL(1,nl))).EQ.4) THEN
CC                 Calculate DL array derivatives from angle change
C                  CALL ANGSCA(NPL(1,0,nl),DL(1,nl),XP,ERROR,*9999)
C                ELSE IF(NBI(NBJ(1,NEL(1,nl))).EQ.5) THEN
CC                 Calculate DL array derivatives from arclength
C                  CALL ARCSCA(IDO,jder,jest,jsca,NBJ,NEL(0,nl),nl,
C     '              NPL(1,0,nl),NPNE,DL,1.0d-6,XP,ERROR,*9999)
C                ELSE IF(NBI(NBJ(1,NEL(1,nl))).EQ.6.
C     '              OR.NBI(NBJ(1,NEL(1,nl))).EQ.7) THEN
CC                 Calculate DL array derivatives from ave. arclength
C                  CALL AVE_ARCSCA(IDO,jder,jest,jsca,NBJ,NEL(0,nl),nl,
C     '              NPL(1,0,nl),NPNE,DL,1.0d-6,XP,ERROR,*9999)
C                ELSE
C                  CALL ASSERT(.FALSE.,'>>Not implemented for '
C     '              //'scale factor type',ERROR,*9999)
C                ENDIF
C                DIFF=DLL(3,nl)-DL(3,nl)
C                IF(DABS(DIFF).GT.MAX_DIFF) THEN
C                  MAX_DIFF=DABS(DIFF)
C                ENDIF
C                IF(DOP) THEN
CC$                call mp_setlock()
C                  WRITE(OP_STRING,'(''nl'',I4,'' old length'',
C     '              F12.8,'' new length'',
C     '              F12.8,'' Difference '',E12.5)')
C     '              nl,DLL(3,nl),DL(3,nl),DIFF
C                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
C                ENDIF
C                DLL(3,nl)=DL(3,nl)
C              ENDDO !nl
C              CONVERGED=.TRUE.
CC ***         Using the L-infinity norm to check convergence
C              IF(MAX_DIFF.GT.ZERO_TOL) CONVERGED=.FALSE.
C            ENDDO !while not converged
C            IF(iter.EQ.100) THEN
C              WRITE(OP_STRING,'('' Iterative calculation of arc lengths'
C     '          //' has not converged in UPSCAl'')')
C              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C            ENDIF
C
CC ***       Calculate SE from DL
C            DO nb=1,NBFT
C              CALL DLSE(IBT,IDO,nb,NEELEM,NLL,NNL,NPL,DL,SE,
C     '          ERROR,*9999)
C            ENDDO
C
C KAT 17Sep98: olde
            ENDDO ! no_nr

          ELSE IF(UNIT) THEN
C           Update the derivative values
            DO no_nr=1,NRLIST(0)
              nr=NRLIST(no_nr)
              DO no_np=1,NPNODE(0,nr)
                np=NPNODE(no_np,nr)
                DO NJL_TYPE=1,3 !geometry/fibre/field
                  DO njj=1,NJ_LOC(NJL_TYPE,0,nr)
                    nj=NJ_LOC(NJL_TYPE,njj,nr)
                    DO nv=1,NVJP(nj,np)
                      DO nk=2,NKJ(nj,np) !skip nodal value
                        SUM=0.0d0
                        NUMBER=0
C                       Find the average scale factor of elements using
C                       this dof. derivative nk at version nv of np
                        DO nep=1,NENP(np,0,nr)
                          ne=NENP(np,nep,nr)
                          nb=NBJ(nj,ne)
                          IF(nb.NE.0) THEN
                            DO nn=1,NNT(nb)
C                             Check if version nv of node np
                              IF(NPNE(nn,nb,ne).EQ.np
     '                          .AND.NVJE(nn,nb,nj,ne).EQ.nv) THEN
                                DO nk_e=1,NKT(nn,nb)
C                                 Check that the same global nk is used
                                  IF(NKJE(nk_e,nn,nj,ne).EQ.nk) THEN
                                    NUMBER=NUMBER+1
                                    ns=NSB(nk_e,nn,nb)
                                    SUM=SUM+SE(ns,nb,ne)
                                  ENDIF !nk
                                ENDDO !nk_e
                              ENDIF !np,nv
                            ENDDO !nn1
                          ENDIF !nb
                        ENDDO !nep1
                        IF(NUMBER.EQ.0) THEN
                          IEND=0
                          CALL APPENDC(IEND,' No element found for ',
     '                      OP_STRING(1))
                          IF(NJL_TYPE.EQ.NJL_FIBR) THEN
                            CALL APPENDC(IEND,'fibre angle ',
     '                        OP_STRING(1))
                          ELSE IF(NJL_TYPE.EQ.NJL_FIEL) THEN
                            CALL APPENDC(IEND,'field ',OP_STRING(1))
                          ELSE
                            CALL APPENDC(IEND,'coordinate ',
     '                        OP_STRING(1))
                          ENDIF
                          CALL APPENDI(IEND,njj,OP_STRING(1))
                          CALL APPENDC(IEND,' node ',OP_STRING(1))
                          CALL APPENDI(IEND,np,OP_STRING(1))
                          CALL APPENDC(IEND,' version ',OP_STRING(1))
                          CALL APPENDI(IEND,nv,OP_STRING(1))
                          CALL APPENDC(IEND,' derivative ',OP_STRING(1))
                          CALL APPENDI(IEND,nk,OP_STRING(1))
C                          WRITE(OP_STRING,
C     '                      '('' No element found for node '',I3,'//
C     '                      ''', version '',I3,'', derivative '',I1)')
C     '                      np,nv,nk
                          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                        ELSE
                          XP(nk,nv,nj,np)=XP(nk,nv,nj,np)*SUM/NUMBER
                        ENDIF !found scale factor
                      ENDDO !nk
                    ENDDO !nv
                  ENDDO !njj
                ENDDO !njj1
              ENDDO !nonode
            ENDDO ! no_nr
C           Change the basis function types
            DO nb=1,NBT
              NBI(nb)=1
            ENDDO !nb
          ENDIF !NORMALISE/UNIT

          IF(ALL_REGIONS) THEN
            CALL LINSCA(IBT,IDO,0,JEST,NBJ,NEELEM,NEL,NLL,NLLINE,
     '        NNL,NPL,NPNE,0,NRE,NVJL,DL,SE,XP,ERROR,*9999)
          ELSE
            DO no_nr=1,NRLIST(0)
              nr=NRLIST(no_nr)
              CALL LINSCA(IBT,IDO,0,JEST,NBJ,NEELEM,NEL,NLL,NLLINE,
     '          NNL,NPL,NPNE,nr,NRE,NVJL,DL,SE,XP,ERROR,*9999)
            ENDDO ! no_nr
          ENDIF

          IF(NORMALISE.AND.GEOMETRY) THEN !update fibres also
C           Divide by new scale factors
            DO no_nr=1,NRLIST(0)
              nr=NRLIST(no_nr)
              DO no_np=1,NPNODE(0,nr)
                np=NPNODE(no_np,nr)
                  DO njj=1,NJ_LOC(NJL_FIBR,0,nr)
                    nj=NJ_LOC(NJL_FIBR,njj,nr)
                    DO nv=1,NVJP(nj,np)
                      DO nk=2,NKJ(nj,np) !skip nodal value
                        FOUND=.FALSE.
C                       Search for an element using the dof.
                        NEPT=NENP(np,0,nr)
                        nep=0
                        DO WHILE(nep.LT.NEPT.AND..NOT.FOUND)
                          nep=nep+1
                          ne=NENP(np,nep,nr)
C                         Search test element for version nv of np
                          nb=NBJ(nj,ne)
                          NNTB=NNT(nb)
                          nn=0
                          DO WHILE(nn.LT.NNTB.AND..NOT.FOUND)
                            nn=nn+1
                            IF(NPNE(nn,nb,ne).EQ.np
     '                        .AND.NVJE(nn,nb,nj,ne).EQ.nv) THEN
C                             Check that the same global nk is used
                              NKTB=NKT(nn,nb)
                              nk_e=0
                              DO WHILE(nk_e.LT.NKTB.AND..NOT.FOUND)
                                nk_e=nk_e+1
                                IF(NKJE(nk_e,nn,nj,ne).EQ.nk) THEN
                                  FOUND=.TRUE.
                                ENDIF !nk
                              ENDDO !nk_e
                            ENDIF !np,nv
                          ENDDO !nn1
                        ENDDO !nep1
                        IF(FOUND) THEN
                          IF(NBI(nb).EQ.3.
     '                      OR.NBI(nb).EQ.6.OR.NBI(nb).EQ.7) THEN
C                           Nodal scale factor
                            ns=NSB(nk_e,nn,nb)
                            XP(nk,nv,nj,np)=XP(nk,nv,nj,np)/SE(ns,nb,ne)
                          ENDIF
                        ENDIF !found scale factor
                      ENDDO !nk
                    ENDDO !nv
                  ENDDO !njj
              ENDDO !no_np
            ENDDO ! no_nr
          ENDIF !NORMALISE.AND.GEOMETRY
        ENDIF !ITERATE

      ENDIF

      CALL EXITS('UPSCAL')
      RETURN
 9999 CALL ERRORS('UPSCAL',ERROR)
      CALL EXITS('UPSCAL')
      RETURN 1
      END


