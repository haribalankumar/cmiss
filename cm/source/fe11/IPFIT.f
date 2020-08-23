      SUBROUTINE IPFIT(IBT,IDO,INP,ISC_GK,ISR_GK,LD,LGE,LN,NBH,
     '  NBJ,NBJF,NEELEM,NELIST,NFFACE,NFLIST,
     '  NHE,NHP,NKEF,NKH,NKHE,NKJE,NMNO,NNF,NPF,NP_INTERFACE,
     '  NPLIST,NPNE,NPNF,NPNODE,NPNY,
     '  nr,NRLIST,NVHE,NVJE,NVJF,NVHP,NVJP,nx,
     '  NYNE,NYNP,NYNR,SE,SF,WU,
     '  XA,XE,XID,XP,ZD,ZP,TYPE,FIX,ERROR,*)

C#### Subroutine: IPFIT
C###  Description:
C###    IPFIT inputs parameters for fibre, field, geometry,
C###    material, sheet or signal fit for region nr.

C#### Variable: KTYP1V
C###  Type: INTEGER
C###  Set_up: IPFIT
C###  Description:
C###    KTYP1V is 0 for no volume conservation and 1 for
C###    penalty based volume constraint on host-mesh.

C#### Variable: KTYP11
C###  Type: INTEGER
C###  Set_up: IPFIT,DEXI
C###  Description:
C###    KTYP11 is 1,2 for element of face fitting

C#### Variable: KTYP12
C###  Type: INTEGER
C###  Set_up: IPFIT
C###  Description:
C###    KTYP12 is 0 for no smoothing / 1 for Sobolev smoothing on the
C###    field / 2 for Sobolev smoothing on deviation from initial field
C###    / 3 for strain energy

C*** See archive version for non-linear and Fourier fitting problems

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'fit000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'solv00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),ISC_GK(NISC_GKM),ISR_GK(NISR_GKM),LD(NDM),
     '  LGE(NHM*NSM,NRCM),LN(0:NEM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NBJF(NJM,NFM),NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NFFACE(0:NF_R_M,NRM),NFLIST(0:NFM),NHE(NEM),NHP(NPM,0:NRM),
     '  NKEF(0:4,16,6,NBFM),NKH(NHM,NPM,NCM,0:NRM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),NMNO(1:2,0:NOPM),
     '  NNF(0:17,6,NBFM),NPF(9,NFM),NP_INTERFACE(0:NPM,0:3),
     '  NPLIST(0:NPM),NPNE(NNM,NBFM,NEM),NPNF(NNM,NBFM),
     '  NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM),nr,
     '  NRLIST(0:NRM),NVHE(NNM,NBFM,NHM,NEM),NVJE(NNM,NBFM,NJM,NEM),
     '  NVJF(NNM,NBFM,NJM),NVHP(NHM,NPM,NCM,0:NRM),NVJP(NJM,NPM),
     '  nx,NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM)
      REAL*8 SE(NSM,NBFM,NEM),SF(NSM,NBFM),
     '  WU(0:NUM+1,NEM),XA(NAM,NJM,NEM),XE(NSM,NJM),XID(NIM,NDM),
     '  XP(NKM,NVM,NJM,NPM),ZD(NJM,NDM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),TYPE*(*)
      LOGICAL FIX(NYM,NIYFIXM)
!     Local Variables
      INTEGER i,IBEG1,ICHAR,IEND1,ILIM,INFO,iy,l,n,N1,nb,nbf,nd,ne,
     '  nef,nf,nh,nhj,nhx,nhx_max,ni,nj,njj,njj2,nk,
     '  NKJF(NKM,NNM,NJM),nm,nn,noelem,noface,nonode,no_nynr,NOQUES,
     '  np,np2,nu,NUTOT,nv,nv2,ny,ny_first
      INTEGER*4 WORK_PTR
      REAL*8 PXI,THETA1,THETA2,X(4),XI(3)
      CHARACTER CHAR1*30,CHAR2*10,CHAR3*10,CHAR4*10
      LOGICAL CALC_SPARSITY,CONTINUE,FACESIG,FILEIP,INLIST,ISFIXED,
     '  ALLSET,USED
      ! Functions
      INTEGER LEN_TRIM

      CALL ENTERS('IPFIT',*9999)

      FACESIG=.TRUE.
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      IF((KTYP8.NE.6).AND.(KTYP8.NE.8).AND.(KTYP8.NE.9)) THEN
C                                     !Not fitting by optimisation

        IDEFLT(1)=1
        IF(IOTYPE.EQ.3) IDATA(1)=ITYP6(nr,nx)
        FORMAT='($,'' Specify whether problem is (1) linear or '//
     '    '(2) nonlinear [1]: '',I1)'
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) ITYP6(nr,nx)=IDATA(1)
      ELSE
        ITYP6(nr,nx)=1
      ENDIF

      nh=0
      NCT(nr,nx)=1

      CALC_GLOBAL(nr,nx)=.TRUE.

C------------------- geometric fitting problem ------------------

      IF(KTYP8.EQ.1) THEN !geometric fitting problem

C LKC 26-MAY-1998 added assert
        CALL ASSERT(NDT.GT.0,'>> No data defined',ERROR,*9999)

        IF(ITYP6(nr,nx).EQ.1) THEN !...is linear

C LKC 25-JUL-2001 Changing the defaults to more appropriate ones,
C  the number of fitting problems is independent of the number of
C  geometric coordinates
C          IDEFLT(1)=NJ_LOC(NJL_FIEL,0,nr)
C          WRITE(CHAR1,'(I1)') IDEFLT(1)

          IDEFLT(1)=1
          WRITE(CHAR1,'(I1)') IDEFLT(1)
          FORMAT='($,'' Specify the number of fitting problems '
     '      //'['//CHAR1(1:1)//']: '',I1)'
          IF(IOTYPE.EQ.3) IDATA(1)=NUM_FIT(0)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,40,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NUM_FIT(0)=IDATA(1)
          CALL ASSERT(NUM_FIT(0).LE.NFITVARMX,
     '      '>>Increase NFITVARMX in fit000.cmn',ERROR,*9999)
          DO njj=1,NUM_FIT(0)
            WRITE(CHAR1,'(I1)') njj

C LKC 25-JUL-2001 Changing the default number of fitted variables
C  to equal the number of fields currently defined.
C            IDEFLT(1)=1
            IDEFLT(1)=NJ_LOC(NJL_FIEL,0,nr)
            WRITE(CHAR2,'(I1)') IDEFLT(1)

            FORMAT='($,'' Specify #geometric vars to be fitted for '
     '        //'problem '//CHAR1(1:1)//' ['//CHAR2(1:1)//']: '',I1)'
            IF(IOTYPE.EQ.3) IDATA(1)=NUM_FIT(njj)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     '        NJ_LOC(NJL_FIEL,0,nr),LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     '        RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) NUM_FIT(njj)=IDATA(1)
            CALL ASSERT(NUM_FIT(njj).LE.NHFITMX,
     '        '>>Increase NHFITMX in fit000.cmn',ERROR,*9999)
            DO nhj=1,NUM_FIT(njj)
              WRITE(CHAR2,'(I1)') nhj
              IDEFLT(1)=nhj
              WRITE(CHAR3,'(I1)') IDEFLT(1)
              FORMAT='($,'' Specify field# to store var '//CHAR2(1:1)
     '          //' of problem '//CHAR1(1:1)
     '          //' ['//CHAR3(1:1)//']: '',I1)'
              IF(IOTYPE.EQ.3) IDATA(1)=NLH_FIT(nhj,2,njj)
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     '          NJ_LOC(NJL_FIEL,0,nr),LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     '          RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) THEN
                NLH_FIT(nhj,1,njj)=NJ_LOC(NJL_FIEL,IDATA(1),nr)
                NLH_FIT(nhj,2,njj)=IDATA(1)
              ENDIF
            ENDDO !nhj
            DO nhj=1,NUM_FIT(njj)
              WRITE(CHAR2,'(I1)') nhj
              IDEFLT(1)=nhj
              WRITE(CHAR3,'(I1)') IDEFLT(1)
              FORMAT='($,'' Specify geometric# to be '
     '          //'fitted for var '//CHAR2(1:1)//' of problem '
     '          //CHAR1(1:1)//' ['//CHAR3(1:1)//']: '',I1)'
              IF(IOTYPE.EQ.3) IDATA(1)=NJ_FIT(nhj,njj)
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     '          NJ_LOC(NJL_GEOM,0,nr),LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     '          RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) THEN
C!!! cpb 7/3/95 Assumes geometric var nj is in the nj position i.e.
C!!! does not use NJ_LOC(NJL_GEOM,...).
                NJ_FIT(nhj,njj)=IDATA(1)
              ENDIF
            ENDDO !nhj
            IF(ITYP10(nr).GT.1) THEN
              DO nd=1,NDT
                CALL ZX(ITYP10(nr),ZD(1,nd),X)
                DO nhj=1,NUM_FIT(njj)
                  nj=NJ_FIT(nhj,njj) !is geo nj to be fitted for fit njj
                  ZD(nj,nd)=X(nj)
                ENDDO !nh
              ENDDO !nd
            ENDIF
            IF(ITYP10(nr).GT.1) THEN !polar coords
              DO nhj=1,NUM_FIT(njj)
                nj=NJ_FIT(nhj,njj)
                IF(nj.NE.1) THEN !not radius
                  DO nd=1,NDT
                    THETA1=ZD(nj,nd) !is theta data value
                    ne=LD(nd)
                    IF(ne.GT.0) THEN !data point projected into an elem.
                      nb=NBJ(nj,ne)
                      DO ni=1,NIT(nb)
                        XI(ni)=XID(ni,nd)
                      ENDDO !ni
                      CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '                  NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '                  SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
                      THETA2=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '                  nb,1,XI,XE(1,nj)) !is interpolated theta
                      IF(DOP) THEN
                        WRITE(OP_STRING,
     '                    '('' Orig. data angle ='',D12.3,'' deg.'','
     '                    //'/'' Interp. angle    ='',D12.3,'' deg.'','
     '                    //'/'' Difference       ='',D12.3,'' deg.'')')
     '                    THETA1*180.0d0/PI,THETA2*180.0d0/PI,
     '                    (THETA1-THETA2)*180.0d0/PI
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                      ENDIF
                      IF(DABS(THETA1-THETA2).GT.PI) THEN !2n*pi correctn
                        IF(DOP) THEN
                          WRITE(OP_STRING,'('' Correction needed..'')')
                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                        ENDIF
                        N=-3
                        CONTINUE=.TRUE.
                        DO WHILE(CONTINUE)
                          N=N+1
                          IF(DABS(THETA1-THETA2+2.0d0*DBLE(N)*PI).LT.PI)
     '                      THEN
                            THETA1=THETA1+2.0d0*DBLE(N)*PI
                            CONTINUE=.FALSE.
                          ELSE IF(N.GT.5) THEN
                            WRITE(OP_STRING,
     '                        '('' >>Warning: Data correction failed '
     '                        //'for nd='',I6)') nd
                            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                            CONTINUE=.FALSE.
                          ENDIF
                        ENDDO
                        ZD(nj,nd)=THETA1 !is new value
                        IF(DOP) THEN
                          WRITE(OP_STRING,'('' New ZD('',I1,'') for nd'
     '                      //' = '',I4,'' is '',D12.3)')
     '                      nj,nd,THETA1
                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                        ENDIF
                      ENDIF
                    ENDIF
                  ENDDO !nd
                ENDIF !not radius
              ENDDO !nhj
            ENDIF !polar coords
          ENDDO !njj

        ELSE IF(ITYP6(nr,nx).EQ.2) THEN !nonlinear geometry fit
          ERROR='>>Not implemented'
          GOTO 9999
        ENDIF

C------------------- fibre/sheet fitting problem ----------------

      ELSE IF(KTYP8.EQ.2) THEN !fibre/sheet fitting problem
        IDEFLT(1)=NJ_LOC(NJL_FIBR,0,nr)
        WRITE(CHAR1,'(I1)') IDEFLT(1)
        FORMAT='($,'' Specify #fit variables in the fit'
     '    //' ['//CHAR1(1:1)//']: '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=NUM_FIT(0)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,40,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) NUM_FIT(0)=IDATA(1)
        CALL ASSERT(NUM_FIT(0).LE.NFITVARMX,
     '    '>>Increase NFITVARMX in fit000.cmn',ERROR,*9999)
        DO njj=1,NUM_FIT(0)
          WRITE(CHAR1,'(I1)') njj
          IDEFLT(1)=1
          WRITE(CHAR2,'(I1)') IDEFLT(1)
          FORMAT='($,'' Specify #fibre vars to be fitted for fit var '
     '      //CHAR1(1:1)//' ['//CHAR2(1:1)//']: '',I1)'
          IF(IOTYPE.EQ.3) IDATA(1)=NUM_FIT(njj)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     '      NJ_LOC(NJL_FIBR,0,nr),LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     '      RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NUM_FIT(njj)=IDATA(1)
          CALL ASSERT(NUM_FIT(njj).LE.NHFITMX,
     '      '>>Increase NHFITMX in fit000.cmn',ERROR,*9999)
          DO nhj=1,NUM_FIT(njj)
            WRITE(CHAR2,'(I1)') nhj
            IDEFLT(1)=nhj
            WRITE(CHAR3,'(I1)') IDEFLT(1)
            FORMAT='($,'' Specify fibre var# to be '
     '        //'fitted for var '//CHAR2(1:1)//' of fit '
     '        //'var '//CHAR1(1:1)//' ['//CHAR3(1:1)//']: '',I1)'
            IF(IOTYPE.EQ.3) IDATA(1)=NLH_FIT(nhj,2,njj)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     '        NJ_LOC(NJL_FIBR,0,nr),LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     '        RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              NLH_FIT(nhj,1,njj)=NJ_LOC(NJL_FIBR,IDATA(1),nr)
              NLH_FIT(nhj,2,njj)=IDATA(1)
            ENDIF
C Specify where data comes from (in ZD or YG array)
            IF(TYPE(1:4).EQ.'DATA') THEN
              IDEFLT(1)=NLH_FIT(nhj,2,njj)
              WRITE(CHAR3,'(I1)') IDEFLT(1)
              FORMAT='($,'' Specify data var# to be fitted '
     '          //'for var '//CHAR2(1:1)//' of fit '
     '          //'var '//CHAR1(1:1)//' ['//CHAR3(1:1)//']: '',I2)'
              IF(IOTYPE.EQ.3) IDATA(1)=NJ_FIT(nhj,njj)-NJT
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     '          NJ_LOC(NJL_FIBR,0,nr),LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     '          RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) NJ_FIT(nhj,njj)=IDATA(1)+NJT
            ELSE IF(TYPE(1:5).EQ.'GAUSS') THEN
              IDEFLT(1)=nhj
              WRITE(CHAR3,'(I1)') IDEFLT(1)
              FORMAT='($,'' Specify Gauss pt var to be '
     '          //'fitted for var '//CHAR2(1:1)//' of fit '
     '          //'var '//CHAR1(1:1)//' ['//CHAR3(1:1)//']: '',I2)'
              IF(IOTYPE.EQ.3) IDATA(1)=NG_FIT(nhj,njj)
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &          NIYGM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &          *9999)
              IF(IOTYPE.NE.3) NG_FIT(nhj,njj)=IDATA(1)
            ENDIF !type
          ENDDO !nhj
        ENDDO !njj
        IF(TYPE(1:5).EQ.'GAUSS') THEN
          l=0
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            l=l+1
            LN(l)=ne
          ENDDO !noelem (ne)
          LN(0)=l
        ENDIF !type

C------------------- field fitting problem ----------------------

      ELSE IF(KTYP8.EQ.3.OR.KTYP8.EQ.8) THEN!field/patch fitting problem
        IDEFLT(1)=NJ_LOC(NJL_FIEL,0,nr)
        WRITE(CHAR1,'(I1)') IDEFLT(1)
        FORMAT='($,'' Specify #fit variables in the fit'
     '    //' ['//CHAR1(1:1)//']: '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=NUM_FIT(0)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,40,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) NUM_FIT(0)=IDATA(1)
        CALL ASSERT(NUM_FIT(0).LE.NFITVARMX,
     '    '>>Increase NFITVARMX in fit000.cmn',ERROR,*9999)
        DO njj=1,NUM_FIT(0)
          WRITE(CHAR1,'(I1)') njj
          IDEFLT(1)=1
          WRITE(CHAR2,'(I1)') IDEFLT(1)
          FORMAT='($,'' Specify #field vars to be fitted for fit var '
     '      //CHAR1(1:1)//' ['//CHAR2(1:1)//']: '',I1)'
          IF(IOTYPE.EQ.3) IDATA(1)=NUM_FIT(njj)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     '      NJ_LOC(NJL_FIEL,0,nr),LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '      INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NUM_FIT(njj)=IDATA(1)
          CALL ASSERT(NUM_FIT(njj).LE.NHFITMX,
     '      '>>Increase NHFITMX in fit000.cmn',ERROR,*9999)
          DO nhj=1,NUM_FIT(njj)
            WRITE(CHAR2,'(I1)') nhj
            IDEFLT(1)=njj
            WRITE(CHAR3,'(I1)') IDEFLT(1)
C Specify where result of fit is to be stored
            FORMAT='($,'' Specify field var# to store var '
     '        //CHAR2(1:1)//' of fit var '
     '        //CHAR1(1:1)//' ['//CHAR3(1:1)//']: '',I1)'
            IF(IOTYPE.EQ.3) IDATA(1)=NLH_FIT(nhj,2,njj)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     '        NJ_LOC(NJL_FIEL,0,nr),LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     '        RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              NLH_FIT(nhj,1,njj)=NJ_LOC(NJL_FIEL,IDATA(1),nr)
              NLH_FIT(nhj,2,njj)=IDATA(1)
            ENDIF
          ENDDO !nh
C Specify where data comes from (in ZD or YG or XP or YQS array)
          IF(TYPE(1:4).EQ.'DATA') THEN
            DO nhj=1,NUM_FIT(njj)
              WRITE(CHAR2,'(I1)') nhj
              IDEFLT(1)=njj
              WRITE(CHAR3,'(I1)') IDEFLT(1)
              FORMAT='($,'' Specify data var# to be fitted '
     '          //'for var '//CHAR2(1:1)//' of fit var '
     '          //CHAR1(1:1)//' ['//CHAR3(1:1)//']: '',I1)'
              IF(IOTYPE.EQ.3) IDATA(1)=NJ_FIT(nhj,njj)-NJT
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &          NJT,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &          *9999)
              IF(IOTYPE.NE.3) NJ_FIT(nhj,njj)=IDATA(1)+NJT
            ENDDO !nhj
          ELSE IF(TYPE(1:4).EQ.'GRID') THEN
            DO nhj=1,NUM_FIT(njj)
              WRITE(CHAR2,'(I1)') nhj
              IDEFLT(1)=njj
              WRITE(CHAR3,'(I1)') IDEFLT(1)
              FORMAT='($,'' Specify grid var# to be fitted '
     '          //'for var '//CHAR2(1:1)//' of fit var '
     '          //CHAR1(1:1)//' ['//CHAR3(1:1)//']: '',I1)'
              IF(IOTYPE.EQ.3) IDATA(1)=1
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &          NIQSM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &          *9999)
              IF(IOTYPE.NE.3) NQ_FIT(nhj,njj)=IDATA(1)
            ENDDO !nhj
          ELSE IF(TYPE(1:10).EQ.'NODE_GROUP') THEN
            FORMAT='($,'' Enter node group name '',I1)'
            IF(IOTYPE.EQ.3) CDATA(1)=GRNAME_NG
            CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,1,
     '          LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) GRNAME_NG=CDATA(1)(1:30)
            DO nhj=1,NUM_FIT(njj)
              WRITE(CHAR2,'(I1)') nhj
              IDEFLT(1)=njj
              WRITE(CHAR3,'(I1)') IDEFLT(1)
              FORMAT='($,'' Specify data var# to be fitted '
     '          //'for var '//CHAR2(1:1)//' of fit var '
     '          //CHAR1(1:1)//' ['//CHAR3(1:1)//']: '',I1)'
              IF(IOTYPE.EQ.3) IDATA(1)=NJ_FIT(nhj,njj)-(NJ_LOC(
     '          NJL_FIEL,1,nr)-1)
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &          NJT,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &          *9999)
              IF(IOTYPE.NE.3) NJ_FIT(nhj,njj)=IDATA(1)+(NJ_LOC(
     '          NJL_FIEL,1,nr)-1)
            ENDDO !nhj
          ELSE IF(TYPE(1:5).EQ.'GAUSS') THEN
c JWF 26/11/02      Set KTYP11=1 for gauss fitting
            KTYP11=1
            DO nhj=1,NUM_FIT(njj)
              WRITE(CHAR2,'(I1)') nhj
              IDEFLT(1)=njj
              WRITE(CHAR3,'(I1)') IDEFLT(1)
              FORMAT='($,'' Specify Gauss pt var to be '
     '          //'fitted for var '//CHAR2(1:1)//' of fit '
     '          //'var '//CHAR1(1:1)//' ['//CHAR3(1:1)//']: '',I2)'
              IF(IOTYPE.EQ.3) IDATA(1)=NG_FIT(nhj,njj)
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &          NIYGM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &          *9999)
              IF(IOTYPE.NE.3) NG_FIT(nhj,njj)=IDATA(1)
            ENDDO !nhj
          ENDIF
        ENDDO !njj
        IF(TYPE(1:5).EQ.'GAUSS') THEN
          l=0
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            l=l+1
            LN(l)=ne
          ENDDO !noelem (ne)
          LN(0)=l
        ENDIF !type

        IF(KTYP8.EQ.8) THEN !patch fitting problem

          IDEFLT(1)=1
          FORMAT='($,'' Specify Gauss point basis function [1]: '',I2)'
          IF(IOTYPE.EQ.3) IDATA(1)=PATCH_BASIS
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NBT,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) PATCH_BASIS=IDATA(1)

          PATCH_INTERNAL=.FALSE.

          IF(NIT(PATCH_BASIS).EQ.2) THEN
            FORMAT='('' Enter polynomial to fit [2]:'''//
     '        '/''   (1) Bilinear P=[1,x,y]'''//
     '        '/''   (2) Bilinear P=[1,x,y,xy]'''//
     '        '/''   (3) Biquadratic'''//
     '        '/$,''    '',I1)'
            IDEFLT(1)=2
            IF(IOTYPE.EQ.3) IDATA(1)=PATCH_POLYNOMIAL
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) PATCH_POLYNOMIAL=IDATA(1)

            FORMAT='($,'' Assemble patch about internal nodes for '
     '        //' boundary nodes [Y]?'',A)'
            IF(IOTYPE.EQ.3) ADATA(1)='Y'
            CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,AYES,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) THEN
              IF(ADATA(1).EQ.'Y') THEN
                PATCH_INTERNAL=.TRUE.
              ENDIF
            ENDIF

          ELSE IF(NIT(PATCH_BASIS).EQ.3) THEN
            FORMAT='('' Enter polynomial to fit [1]:'''//
     '        '/''   (1) trilinear '''//
     '        '/$,''    '',I1)'
            IDEFLT(1)=1
            IF(IOTYPE.EQ.3) IDATA(1)=PATCH_POLYNOMIAL
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,1,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) PATCH_POLYNOMIAL=IDATA(1)
          ELSE
            CALL ASSERT(.FALSE.,'>>1D patches are not implemented',
     '        ERROR,*9999)
          ENDIF

          FORMAT='('' Enter Gauss point pattern for patch [1]:'''//
     '      '/''   (1) All Gauss points'''//
     '      '/''   (2) 4 Gauss points from 3x3 Gauss points'''//
     '      '/''   (3) 8 Gauss points from 3x3x3 Gauss points'''//
     '      '/$,''    '',I1)'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=IDEFLT(1)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) PATCH_PATTERN=IDATA(1)
        ENDIF

C------------------- signal fitting problem ---------------------

      ELSE IF(KTYP8.EQ.4) THEN !signal fitting problem

        NUM_FIT(0)=1
        NUM_FIT(1)=1
        IDEFLT(1)=1
        FORMAT='($,'' Specify the field variable to fit [1]: '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=NLH_FIT(1,2,1)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     '    NJ_LOC(NJL_FIEL,0,nr),LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '    INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          NLH_FIT(1,1,1)=NJ_LOC(NJL_FIEL,IDATA(1),nr)
          NLH_FIT(1,2,1)=IDATA(1)
        ENDIF
        NJ_FIT(1,1)=NJT+1

        FORMAT='(/$,'' Enter the data (electrode) region number '
     '    //'[1]: '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=DATA_REGION
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,9,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) DATA_REGION=IDATA(1)
        FORMAT='(/$,'' Enter the sample skip number [1]: '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=SIGNAL_SKIPSAMPLE
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,99,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) SIGNAL_SKIPSAMPLE=IDATA(1)
C LEO 18-OCT-97
C***    output level for fitting signal

        FORMAT='('' Specify output for fitting signal '
     '    //'[0]: '''//
     '    '/''   (0) No output'''//
     '    '/''   (1) Completed fits and timing'''//
     '    '/''   (2) & Current fit'''//
     '    '/''   (3) & Data and fit information'''//
     '    '/$,''    '',I1)'
        IDEFLT(1)=0
        IF(IOTYPE.EQ.3) IDATA(1)=SIGNAL_OUTPUT
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,3,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) SIGNAL_OUTPUT=IDATA(1)

        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CDEFLT(1)=FILE00(IBEG1:IEND1)
        FORMAT='(/$,'' Enter the input signal filename [current]: '',A)'
        IF(IOTYPE.EQ.3) CDATA(1)=INFILENAME(1:100)
        CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,100,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) INFILENAME=CDATA(1)(1:100)
        FORMAT='($,'' Enter the output history filename [current]: '','
     '    //'A)'
        IF(IOTYPE.EQ.3) CDATA(1)=OUTFILENAME(1:100)
        CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,100,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) OUTFILENAME=CDATA(1)(1:100)

        l=0
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          l=l+1
          LN(l)=ne
        ENDDO !noelem (ne)
        LN(0)=l

C------------------- motion fitting with Fourier basis ----------

      ELSE IF(KTYP8.EQ.5) THEN !motion fitting with Fourier basis

        ERROR='>>Not implemented'
        GOTO 9999

C------------------- data fitting with optimisation ----------

      ELSE IF(KTYP8.EQ.6) THEN !data fitting with optimisation

        CALL ASSERT(NHM.GE.NJM,'>>Increase NHM (must be >= NJM)',
     '    ERROR,*9999)
        NUM_FIT(0)=1
        NUM_FIT(1)=NJT
        DO i=1,2
          DO nhj=1,NJT
            NLH_FIT(nhj,i,1)=nhj
          ENDDO !nhj
        ENDDO !i
        DO nhj=1,NJT
          NJ_FIT(nhj,1)=nhj
        ENDDO !nhj

C------------------- material fitting problem -------------------

      ELSE IF(KTYP8.EQ.7) THEN !material fitting problem
        NUM_FIT(0)=1 !1 fitting problem only
        CALL ASSERT(NUM_FIT(0).LE.NFITVARMX,
     '    '>>Increase NFITVARMX in fit000.cmn',ERROR,*9999)

        FORMAT='($,'' Specify number of sets of residuals [1]: '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=KTYP28
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NIYM-10,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) KTYP28=IDATA(1)

        IDEFLT(1)=NJ_LOC(NJL_FIEL,0,nr)
        WRITE(CHAR1,'(I1)') IDEFLT(1)
        FORMAT='($,'' Specify #material parameters in the fit'
     '    //' ['//CHAR1(1:1)//']: '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=NUM_FIT(1)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,9,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) NUM_FIT(1)=IDATA(1)
        CALL ASSERT(NUM_FIT(1).LE.NHFITMX,
     '    '>>Increase NHFITMX in fit000.cmn',ERROR,*9999)
        NMNO(1,0)=NUM_FIT(1)

        CHAR1=' '
        DO nhj=1,NMNO(1,0)
          IDEFLT(nhj)=nhj
          IBEG1=2*nhj-1
          IEND1=2*nhj
          WRITE(CHAR1(IBEG1:IEND1),'(I2)') IDEFLT(nhj)
        ENDDO !nhj
        CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
        FORMAT='($,'' Enter list of material params in fit ['
     '    //CHAR1(IBEG1:IEND1)//']: '',20I3)'
        IF(IOTYPE.EQ.3) THEN
          DO nhj=1,NMNO(1,0)
            IDATA(nhj)=NMNO(1,nhj)
          ENDDO !nhj
        ENDIF
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '    NMNO(1,0),ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,20,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          DO nhj=1,NMNO(1,0)
            NMNO(1,nhj)=IDATA(nhj)
          ENDDO !nhj
        ENDIF

        DO nhj=1,NUM_FIT(1) !loop over material parameters
          nm=NMNO(1,nhj)
          IDEFLT(1)=nm
          WRITE(CHAR3,'(I1)') IDEFLT(1)
C Specify where result of fit is to be stored
          FORMAT='($,'' Specify field var# to store material param#'
     '        //CHAR3(1:1)//' ['//CHAR3(1:1)//']: '',I1)'
          IF(IOTYPE.EQ.3) IDATA(1)=NLH_FIT(nhj,2,1)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     '      NJ_LOC(NJL_FIEL,0,nr),LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     '      RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            NLH_FIT(nhj,1,1)=NJ_LOC(NJL_FIEL,IDATA(1),nr)
            NLH_FIT(nhj,2,1)=IDATA(1)
          ENDIF

C KAT 23Dec99: For appropriate ny maps need this.
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            nb=NBJ(NLH_FIT(nhj,1,1),ne)
            CALL ASSERT(NKT(0,nb).EQ.1,'Field must be piecewise linear',
     '        ERROR,*9999)
          ENDDO !noelem (ne)
        ENDDO !nhj


C------------------- SPLINE signal fitting problem ---------------------

      ELSE IF(KTYP8.EQ.9) THEN !signal fitting problem

        NUM_FIT(0)=1
        NUM_FIT(1)=1

        IDEFLT(1)=1
        FORMAT='($,'' Specify the field variable to fit [1]: '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=NLH_FIT(1,2,1)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     '    NJ_LOC(NJL_FIEL,0,nr),LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '    INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          NLH_FIT(1,1,1)=NJ_LOC(NJL_FIEL,IDATA(1),nr)
          NLH_FIT(1,2,1)=IDATA(1)
        ENDIF
        NJ_FIT(1,1)=NJT+1

        FORMAT='(/$,'' Enter the data (electrode) region number '
     '    //'[1]: '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=DATA_REGION
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,9,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) DATA_REGION=IDATA(1)
        FORMAT='(/$,'' Enter the sample skip number [1]: '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=SIGNAL_SKIPSAMPLE
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,99,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) SIGNAL_SKIPSAMPLE=IDATA(1)

        FORMAT='('' Specify output for fitting signal '
     '    //'[0]: '''//
     '    '/''   (0) No output'''//
     '    '/''   (1) Horizontal spline information '''//
     '    '/''   (2) & Complete fit timing'''//
     '    '/''   (3) & Current fit'''//
     '    '/''   (4) & Vertical spline information '''//
     '    '/$,''    '',I1)'
        IDEFLT(1)=0
        IF(IOTYPE.EQ.3) IDATA(1)=SIGNAL_OUTPUT
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,4,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) SIGNAL_OUTPUT=IDATA(1)

        CALL STRING_TRIM(FILE00,IBEG1,IEND1)
        CDEFLT(1)=FILE00(IBEG1:IEND1)
        FORMAT='(/$,'' Enter the input signal filename [current]: '',A)'
        IF(IOTYPE.EQ.3) CDATA(1)=INFILENAME(1:100)
        CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,100,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) INFILENAME=CDATA(1)(1:100)
        FORMAT='($,'' Enter the output history filename [current]: '','
     '    //'A)'
        IF(IOTYPE.EQ.3) CDATA(1)=OUTFILENAME(1:100)
        CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,100,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) OUTFILENAME=CDATA(1)(1:100)

C*** add additional information here
      ENDIF !ktyp8 setup

C----------------------------------------------------------------

C***  Set up nh numbers
      nhx_max=0
      DO njj=1,NUM_FIT(0)
        DO nhj=1,NUM_FIT(njj)
          nhx_max=nhx_max+1
          NLH_FIT(nhj,3,njj)=nhx_max

          CALL ASSERT(nhx_max.LE.40,'>>Maximum number of nhx is 40',
     '      ERROR,*9999)

C          CALL ASSERT(nhx_max.LE.9,
C     '      '>>Increase second dimension of NH_FIT in fit000.cmn',
C     '      ERROR,*9999)
C          NH_FIT(1,nhx_max)=njj
C          NH_FIT(2,nhx_max)=nhj
        ENDDO !nhj
      ENDDO !njj

C***  Set up NH_LOC(nhx,nx)
      CALL CALC_NH_LOC(nhx_MAX,nx,ERROR,*9999)

C***Set up dependent variable interpolation information
      DO njj=1,NUM_FIT(0)
        DO nhj=1,NUM_FIT(njj)
          nj=NLH_FIT(nhj,1,njj)
          nhx=NLH_FIT(nhj,3,njj)
          nh=NH_LOC(nhx,nx)
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            NHP(np,nr)=nhx_max
            NVHP(nh,np,1,nr)=NVJP(nj,np)
          ENDDO !nonode (np)
C PM 14Aug02 : No need to treat faces separately at this stage
C          IF(KTYP11.EQ.1) THEN !element fitting
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            nb=NBJ(nj,ne)
            NBH(nh,1,ne)=nb
            NHE(ne)=nhx_max
            DO nn=1,NNT(nb)
              NVHE(nn,nb,nh,ne)=NVJE(nn,nb,nj,ne)
              DO nk=1,NKT(nn,nb)
                NKHE(nk,nn,nh,ne)=nk
              ENDDO !nk
            ENDDO !nn
          ENDDO !noelem (ne)
C new 13/1/2002 implementing face fitting
C From here on the dependent variable elements are actually faces
C This might cause some memory alloc problems because there are typically
C more faces than elements, but we will try to catch them on the way.

C          ELSE IF(KTYP11.EQ.2) THEN !face fitting
C            DO noface=1,LN(0)
C              nf=LN(noface)
C              nb=NBJF(nj,nf)
C              NBH(nh,1,nf)=nb
C              NHE(nf)=NUM_FIT(njj)
C              ne=NPF(6,nf)
C              nef=NPF(8,nf)
C              CALL CALC_FACE_INFORMATION_IND(NBJ(1,ne),NBJF(1,nf),nef,
C     '          NKJE(1,1,1,ne),NKEF,NKJF,NNF,NPNE(1,1,ne),NPNF,nr,
C     '          NVJE(1,1,1,ne),NVJF,SE(1,1,ne),SF,ERROR,*9999)
C              DO nn=1,NNT(nb)
C                NVHE(nn,nb,nh,nf)=NVJF(nn,nb,nj)
C                DO nk=1,NKT(nn,nb)
C                  NKHE(nk,nn,nh,nf)=nk
C                ENDDO !nk
C              ENDDO !nn
C            ENDDO !noface
C          ENDIF

        ENDDO !nhj
      ENDDO !njj

C LKC 30-AUG-98 NHP not dropped down correctly
C cpb 5/2/97 Adding calc_nkh_fit
C      CALL CALC_NKH_FIT(NBH,NEELEM(0,nr),NHP,NKH,NPNE,NPNODE(0,nr),
      CALL CALC_NKH_FIT(NBH,NEELEM(0,nr),NHP(1,nr),
     '  NKH,NPNE,NPNODE(0,nr),nr,nx,ERROR,*9999)

      IF(KTYP8.NE.6) THEN
        DO njj=1,NUM_FIT(0)
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            DO nhj=1,NUM_FIT(njj)
              nj=NJ_FIT(nhj,njj)
              nhx=NLH_FIT(nhj,3,njj)
              nh=NH_LOC(nhx,nx)
              DO nv=1,NVHP(nh,np,1,nr)
                DO nk=1,NKH(nh,np,1,nr)
                  IF(TYPE(1:4).EQ.'DATA') THEN
                    ZP(nk,nv,nh,np,1)=XP(nk,nv,nj,np)
                  ELSE IF(TYPE(1:5).EQ.'GAUSS') THEN
C CS 16/4/2003 Initialising Gauss point fitting the same as
C              data fitting. Can't see why they are treated
C              differently.
C                    ZP(nk,nv,nh,np,1)=0.0d0
                    ZP(nk,nv,nh,np,1)=XP(nk,nv,NLH_FIT(nhj,1,njj),np)
                  ELSE IF(TYPE(1:4).EQ.'GRID') THEN
                    ZP(nk,nv,nh,np,1)=0.0d0
                  ENDIF
                ENDDO !nk
              ENDDO !nv
            ENDDO !nhj
          ENDDO !nonode (np)
        ENDDO !njj

      ENDIF

C LKC 22-MAR-1998
C      IF(KTYP8.NE.8) THEN

      IF(KTYP8.LT.8) THEN !now exludes  KTYP8=9 for splines
        IDEFLT(1)=0
        IF(KTYP8.EQ.6) THEN !Fitting by optimisation
          FORMAT='('' Enter smoothing type [0]:'''
     '      //'/''   (0) None'''
     '      //'/''   (1) Sobolev'''
     '      //'/$,''    '',I1)'
          ILIM=1
        ELSE
          FORMAT='('' Enter smoothing type [0]:'''
     '      //'/''   (0) None'''
     '      //'/''   (1) Sobolev on field'''
     '      //'/''   (2) Sobolev on deviation from initial field'''
     '      //'/''   (3) Strain energy'''
     '      //'/$,''    '',I1)'
          ILIM=3
        ENDIF
        IF(IOTYPE.EQ.3) THEN
          IDATA(1)=KTYP12
        ENDIF
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,ILIM,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) KTYP12=IDATA(1)
      ELSE
        KTYP12=0
      ENDIF

      IF(KTYP12.EQ.0.AND.LN(0).EQ.0) THEN
        l=0
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          CONTINUE=.TRUE.
          nd=1
          DO WHILE(CONTINUE.AND.nd.LE.NDT)
            IF(LD(nd).EQ.ne) THEN
              l=l+1
              LN(l)=ne
              CONTINUE=.FALSE.
            ENDIF
            nd=nd+1
          ENDDO
        ENDDO
        LN(0)=l
      ELSE IF(LN(0).EQ.0)THEN
        l=0
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
C??? cpb 13/12/94 Don't quite understand what this line is trying to do
C          IF(NIT(NBJ(1,ne)).EQ.NIOT)THEN
          l=l+1
          LN(l)=ne
C          ENDIF
        ENDDO
        LN(0)=l
      ENDIF
      IF(DOP) THEN
        WRITE(CHAR1,'(I10)') LN(0)
        CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
        FORMAT='(/'' The '//CHAR1(IBEG1:IEND1)//
     '    ' elements involved in the optimisation are:'')'
        WRITE(OP_STRING,FORMAT)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(1X,10(1X,I5),/:(1X,10(1X,I5)))')
     '    (LN(L),L=1,LN(0))
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

C *** Specify smoothing constraints on each element

      IF(KTYP12.EQ.1.OR.KTYP12.EQ.2) THEN !Sobolev smoothing
C GM 14/7/96 Input factors with generalised element number prompt
C        DO nu=1,NUM
C          RDEFLT(nu)=0.0d0
C        ENDDO
C        DO L=1,LN(0)
C          ne=LN(L)
C          WU(0,ne)=1.0d0 !the scaling factor for the Sobolev weights
C          WRITE(CHAR4,'(I5)') ne
C          CALL STRING_TRIM(CHAR4,IBEG1,IEND1)
C          FORMAT='('' For element '//CHAR4(IBEG1:IEND1)//':'')'
C          CALL GINOUT(IOTYPE,0,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
C     '      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,
C     '      IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C          IF(NIT(NBH(NH_LOC(1,nx),1,NEELEM(1,nr))).EQ.1) THEN
C            NUTOT=2
C            FORMAT='('' The 2 weights on derivs wrt Xi_1/_11/'
C     '        //'are [prev]:'''
C     '        //'/$,''    '',2D9.2)'
C          ELSE IF(NIT(NBH(NH_LOC(1,nx),1,NEELEM(1,nr))).EQ.2) THEN
C            NUTOT=5
C            FORMAT='('' The 5 weights on derivs wrt Xi_1/_11/_2/_22/'
C     '        //'_12 are [prev]:'''
C     '        //'/$,''    '',5D9.2)'
C          ELSE IF(NIT(NBH(NH_LOC(1,nx),1,NEELEM(1,nr))).EQ.3) THEN
C            NUTOT=10
C            FORMAT='('' The 10 weights on derivs wrt Xi_1/_11/_2/_22/'
C     '        //'_12/_3/_33/_23/_31/_123 are [prev]:'''
C     '        //'/$,''    '',10D9.2)'
C          ENDIF
C          IF(IOTYPE.EQ.3) THEN
C            DO nu=1,NUTOT
C              RDATA(1)=WU(nu,ne)
C            ENDDO
C          ENDIF
C          CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
C     '      NUTOT,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,
C     '      IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
C          IF(IOTYPE.NE.3) THEN
C            DO nu=1,NUTOT
C              WU(nu,ne)=RDATA(nu)
C              RDEFLT(nu)=RDATA(nu)
C            ENDDO
C          ENDIF
C          IF(DOP) THEN
C            WRITE(OP_STRING,'('' WU(0..,ne)=''/,10D9.2)')
C     '        (WU(nu,ne),NU=0,NUTOT)
C            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C          ENDIF
C        ENDDO
C GMH 14/7/96 We will use WU(0,ne) as the trigger to see if the user
C has named all the elements.  Set to 0 first, then to 1.0 per element.
C Check at the then to make sure all are 1.0
        DO L=1,LN(0)
          ne=LN(L)
          WU(0,ne)=0.0d0 !the scaling factor for the Sobolev weights
        ENDDO !L
        DO nu=1,10 !maximum number
          RDEFLT(nu)=0.0d0
        ENDDO

        IF(KTYP11.EQ.1) THEN ! Element fitting
C         Get the elements
          noelem=0
 6720     FORMAT='($,'' Enter element #s/name [EXIT]: '',I5)'
          IF(IOTYPE.EQ.3) THEN
            noelem=noelem+1
            IF(noelem.LE.LN(0)) THEN
              ne=LN(noelem)
              IDATA(1)=ne
            ELSE
              IDATA(1)=0
            ENDIF
          ENDIF

 6760     CDATA(1)='ELEMENTS' !for use with group input
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,
     '      0,NET(nr),
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IDATA(1).NE.0) THEN !not default exit
            NELIST(0)=IDATA(0)
            DO n=1,IDATA(0)
              NELIST(n)=IDATA(n)
              ne=IDATA(n)
C             Make sure we are in the region
              IF(.NOT.INLIST(ne,NEELEM(1,nr),
     '          NEELEM(0,nr),N1)) THEN
                WRITE(OP_STRING,'('' >>WARNING: Element '',I5,'' is '
     '            //'not in the current region'')') ne
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                GOTO 6760
              ENDIF
            ENDDO !n

C           Get type for first element in group
            ne=NELIST(1) !rest of group filled at end of loop
            IF(NIT(NBH(NH_LOC(1,nx),1,NEELEM(1,nr))).EQ.1) THEN
              NUTOT=2
              FORMAT='('' The 2 weights on derivs wrt Xi_1/_11/'
     '          //'are [prev]:'''
     '          //'/$,''    '',2D9.2)'
            ELSE IF(NIT(NBH(NH_LOC(1,nx),1,NEELEM(1,nr))).EQ.2) THEN
              NUTOT=5
              FORMAT='('' The 5 weights on derivs wrt Xi_1/_11/_2/_22/'
     '          //'_12 are [prev]:'''
     '          //'/$,''    '',5D9.2)'
            ELSE IF(NIT(NBH(NH_LOC(1,nx),1,NEELEM(1,nr))).EQ.3) THEN
              NUTOT=10
              FORMAT='('' The 10 weights on derivs wrt Xi_1/_11/_2/_22/'
     '          //'_12/_3/_33/_23/_31/_123 are [prev]:'''
     '          //'/$,''    '',10D9.2)'
            ENDIF
C           Set the default
            IF(IOTYPE.EQ.3) THEN
              DO nu=1,NUTOT
                RDATA(nu)=WU(nu+1,ne)
              ENDDO
            ENDIF
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,NUTOT,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &        IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,
     '        INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              DO nu=1,NUTOT
                RDEFLT(nu)=RDATA(nu)
              ENDDO
            ENDIF

C           Apply type to all of elements group
            DO n=1,NELIST(0)
              ne=NELIST(n)
C             Are we in the list of elements?
              IF(INLIST(ne,LN(1),LN(0),N1)) THEN
                WU(0,ne)=1.0d0
                DO nu=1,NUTOT
                  WU(nu+1,ne)=RDATA(nu)
                ENDDO
              ELSE
C               Just warn the user, and don't set the value
                WRITE(OP_STRING,'('' Element '',I5,'' is not '
     '            //'in the fit'')') ne
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO !n

            GO TO 6720 !for more elements
          ENDIF !idata(1).ne.0
          ALLSET=.TRUE.
          DO L=1,LN(0)
            ne=LN(L)
            IF(WU(0,ne).EQ.0.0d0) THEN
              ALLSET=.FALSE.
            ENDIF
            IF(DOP) THEN
              WRITE(OP_STRING,'('' WU(0..,'',I5,'')=''/,10D9.2)')
     '          ne,(WU(nu,ne),NU=0,NUTOT+1)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO !L
          CALL ASSERT(ALLSET,'>>Must set Sobolev factors for all'
     '      //' elements',ERROR,*9999)

C new CS 13/1/2002 implementing face fitting
        ELSE IF(KTYP11.EQ.2) THEN ! Face fitting
C         Get the faces
          noface=0
 7720     FORMAT='($,'' Enter face #s/name [EXIT]: '',I5)'
          IF(IOTYPE.EQ.3) THEN
            noface=noface+1
            IF(noface.LE.LN(0)) THEN
              nf=LN(noface)
              IDATA(1)=nf
            ELSE
              IDATA(1)=0
            ENDIF
          ENDIF
 7760     CDATA(1)='FACES' !for use with group input
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,
     '      0,NFT,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IDATA(1).NE.0) THEN !not default exit
            NFLIST(0)=IDATA(0)
            DO n=1,IDATA(0)
              NFLIST(n)=IDATA(n)
              nf=IDATA(n)
C             Make sure face exists
              IF(.NOT.INLIST(nf,NFFACE(1,nr),
     '          NFFACE(0,nr),N1)) THEN
                WRITE(OP_STRING,'('' >>WARNING: Face '',I5,'' does '
     '            //'not exist in the current region'')') nf
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                GOTO 7760
              ENDIF
            ENDDO !n

C           Get type for first face in group
            nf=NFLIST(1) !rest of group filled at end of loop
            NUTOT=5
            FORMAT='('' The 5 weights on derivs wrt direction'
     '        //'_1/_11/_2/_22/_12 are [prev]:'''
     '        //'/$,''    '',5D9.2)'
C           Set the default
            IF(IOTYPE.EQ.3) THEN
              DO nu=1,NUTOT
                RDATA(nu)=WU(nu+1,nf)
              ENDDO
            ENDIF
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,NUTOT,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &        IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,
     '        INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              DO nu=1,NUTOT
                RDEFLT(nu)=RDATA(nu)
              ENDDO
            ENDIF

C           Apply type to all of faces group
            DO n=1,NFLIST(0)
              nf=NFLIST(n)
              IF(INLIST(nf,NFFACE(1,nr),
     '          NFFACE(0,nr),N1)) THEN
                WU(0,nf)=1.0d0
                DO nu=1,NUTOT
                  WU(nu+1,nf)=RDATA(nu)
                ENDDO
              ELSE
C               Warn the user
                WRITE(OP_STRING,'('' Face '',I5,'' does not '
     '            //'exist in the current region'')') nf
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO !n

            GO TO 7720 !for more faces
          ENDIF !idata(1).nf.0
          ALLSET=.TRUE.
          DO L=1,LN(0)
            nf=LN(L)
            IF(WU(0,nf).EQ.0.0d0) THEN
              ALLSET=.FALSE.
            ENDIF
            IF(DOP) THEN
              WRITE(OP_STRING,'('' WU(0..,'',I5,'')=''/,10D9.2)')
     '          nf,(WU(nu,nf),NU=0,NUTOT+1)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO !L
          CALL ASSERT(ALLSET,'>>Must set Sobolev factors for all'
     '      //' faces',ERROR,*9999)

        ENDIF !face fitting

      ELSE IF(KTYP12.EQ.3) THEN !Strain energy smoothing
        DO L=1,LN(0)
          ne=LN(L)
          WRITE(CHAR4,'(I5)') ne
          CALL STRING_TRIM(CHAR4,IBEG1,IEND1)
          FORMAT='($,'' Enter weight for element '//CHAR4(IBEG1:IEND1)
     '      //' [0]: '',D12.3)'
          IF(IOTYPE.EQ.3) RDATA(1)=WU(2,ne)
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            WU(0,ne)=1.0d0 !total weight
            WU(2,ne)=RDATA(1) !Xi(1) weight
            WU(4,ne)=RDATA(1) !Xi(2) weight
            WU(7,ne)=RDATA(1) !Xi(3) weight
            WU(1,ne)=0.0d0
            WU(3,ne)=0.0d0
            WU(5,ne)=0.0d0
            WU(6,ne)=0.0d0
            WU(8,ne)=0.0d0
            WU(9,ne)=0.0d0
            WU(10,ne)=0.0d0
          ENDIF
        ENDDO !l (ne)
        CALL ASSERT(NHM.GE.NJM,'>>NHM must be >= NJM',ERROR,*9999)
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          DO nhx=1,NHP(np,nr)
            nh=NH_LOC(nhx,nx)
            DO nv=1,NVHP(nh,np,1,nr)
              DO nk=1,NKH(nh,np,1,nr)
                ZP(nk,nv,nh,np,2)=ZP(nk,nv,nh,np,1)
              ENDDO !nk
            ENDDO !nv
          ENDDO !nh
        ENDDO !nonode (np)
        WRITE(OP_STRING,'('' >>WARNING: Deformed coordinates copied '
     '    //'to nc=2 of ZP'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ELSE
C CPB 25/5/97 set WU(0,...) to 1.0 as it is used for fitting without
C smoothing
        DO l=1,LN(0)
          ne=LN(l)
          WU(0,ne)=1.0d0
        ENDDO
      ENDIF

C JWF 10/08/04 Question for volume conservation on host-mesh fitting.
C Only asked if ktyp8=1 (geometric fitting) and 
C NRT>=2 (2 or more regions), which is required for host-mesh.
C KTYP1V is defined here. Is 1/0 for volume/no constraint respectively.
C VOL_WT is setup here and contains the volume constraint penalty. (fit000.cmn)

      KTYP1V=0 !initialise
      VOL_WT=0.0d0 !initialise    
      IF((KTYP8.EQ.1).AND.(NRT.GE.2)) THEN !volume constraint
        ADEFLT(1)='N'
        FORMAT='(/$,'' Do you want volume conservation '
     &    //' [N]? '',A)'
        IF(IOTYPE.EQ.3) THEN
          ADATA(1)='N'
        ENDIF
        CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,
     &    FILEIP,FORMAT,1,ADATA,ANO,CDATA,CDEFLT,
     &    ICHAR,IDATA,IDEFLT,0,IMAX,LDATA,LDEFLT,
     &    RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          IF(ADATA(1).EQ.'Y'.OR.ADATA(1).EQ.'y') THEN
            KTYP1V=1           
          ELSE
            KTYP1V=0
          ENDIF
        ENDIF
        IF(KTYP1V.EQ.1) THEN ! get penalty value
          RDEFLT(1)=1.0d0
          FORMAT='($,'' Enter the penalty weight [1.0]: '',D5.2)'
          IF(IOTYPE.EQ.3) THEN
            RDATA(1)=1.0d0
            VOL_WT=RDATA(1)
          ENDIF
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     &      LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            VOL_WT=RDATA(1)
          ENDIF
        ENDIF ! set penalty
      ENDIF !volume constraint

C KAT 15/11/00: If cross derivatives are to be set to zero then this
C     should be specified in the basis function file.  See documentation
C     on KTYP93.
CC AJPs 17/8/99 Adding the ability fo fit bicubic hermite fields for
CC      later use with BEM transfer matrices.
C      HERMITE_3D = .FALSE.
C      IF(KTYP8.EQ.4) THEN !signal fitting
C! Check if there are any bicubic hermite elements.  If so determine
C! whether cross derivative coefficients are to be included in the
C! interpolation of the dependent variable (if so BE matrix system
C! likely to be poorly conditioned).
C
C        HERMITE_3D = .FALSE.
C        DO noelem=1,NEELEM(0,nr)
C          ne=NEELEM(noelem,nr) !This identifies an elem.
C          nb=NBASEF(NBH(NH_LOC(NLH_FIT(1,3,1),nx),1,ne),1)
C          !First dependent variable
CC AJP 19/8/99
C!          IF(NJ_LOC(NJL_FIEL,0,nr).EQ.3) THEN
CC AJPe 19/8/99
C          IF(NIM.EQ.1) THEN
C            IF(IBT(1,1,nb).EQ.3.AND.IBT(2,1,nb).EQ.4) THEN
C              !Hermite in each direction or hermite simplex or hermite sec
C              HERMITE_3D=.TRUE.
C              GOTO 100
C            ENDIF
C          ELSE
C            IF((IBT(1,1,nb).EQ.2.AND.IBT(1,2,nb).EQ.2).OR.
C     '        (IBT(1,1,nb).EQ.3.AND.IBT(2,1,nb).EQ.4).OR.
C     '        (((IBT(1,1,nb).EQ.5.OR.IBT(1,1,nb).EQ.6).AND.
C     '        IBT(2,1,nb).EQ.4).AND.IBT(1,2,nb).EQ.2).OR.
C     '        (((IBT(1,2,nb).EQ.5.OR.IBT(1,2,nb).EQ.6).AND.
C     '        IBT(2,2,nb).EQ.4).AND.IBT(1,1,nb).EQ.2)) THEN
C              !Hermite in each direction or hermite simplex or hermite sec
C              HERMITE_3D=.TRUE.
C              GOTO 100
C            ENDIF
C          ENDIF
C!        ENDIF
C        ENDDO !noelem
C 100    CONTINUE
C        IF(HERMITE_3D) THEN
C          FORMAT='(/$,'' Do you want to set cross derivatives of'//
C     '      ' the dependent variables to zero [Y]? '',A)'
C          IF(IOTYPE.EQ.3) THEN
C            IF(KTYP93(1,nr).EQ.0) THEN
C              ADATA(1)='N'
C            ELSE
C              ADATA(1)='Y'
C            ENDIF
C          ENDIF
C          CALL GINOUT(IOTYPE,1,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
C     '      1,ADATA,AYES,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
C     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C          IF(ADATA(1).EQ.'Y') THEN
C            KTYP93(1,nr)=1
C          ELSE
C            KTYP93(1,nr)=0
C          ENDIF
C        ELSE
C          KTYP93(1,nr)=0
C        ENDIF
C      ELSE
C        KTYP93(1,nr)=0
C      ENDIF
CC AJPe 17/8/99

C*** Calculate ny maps
      IF(KTYP8.NE.8) THEN
        CALL CALC_NY_MAPS_IND(NKH,NP_INTERFACE,NPNODE,NPNY,nr,NVHP,nx,
     '    NYNP,NYNR,ERROR,*9999)
      ENDIF

      IF((KTYP8.NE.6).AND.(KTYP8.NE.8).AND.(KTYP8.NE.9)) THEN
        IF(nr.EQ.NRLIST(NRLIST(0))) THEN
C         only calc sparsity after ny's for last region have been set up
          IF(KTYP8.EQ.7) THEN !material fit
            ADEFLT(1)='N'
          ELSE
            ADEFLT(1)='Y'
          ENDIF
          FORMAT='(/$,'' Do you want the global matrices stored as '
     '      //'sparse matrices ['//ADEFLT(1)(1:1)//']? '',A)'
          IF(IOTYPE.EQ.3) THEN
            IF((KTYP24.EQ.0).OR.(KTYP8.EQ.7)) THEN
              ADATA(1)='N'
            ELSE
              ADATA(1)='Y'
            ENDIF
          ENDIF
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(ADATA(1).EQ.'Y') THEN
            KTYP24=1
            IF(IOTYPE.NE.3) CALL ASSERT(USE_SPARSE.NE.0,
     '        '>>Set USE_SPARSE to 1 to use sparse matrices',
     '        ERROR,*9999)
            IF(KTYP8.EQ.7) THEN
              WRITE(OP_STRING,'('' >>Warning: material fit matricies '
     '          //'are not sparse'')')
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ELSE
            KTYP24=0
          ENDIF

          CALC_SPARSITY=.FALSE.
          IF(KTYP24.EQ.1) THEN !using sparse matrix storage
            FORMAT='($,'' Do you want to calculate the sparsity '
     '        //'pattern for the global matrices [Y]? '',A)'
            ADEFLT(1)='Y'
            IF(IOTYPE.EQ.3) THEN
              ADATA(1)='Y'
              CALC_SPARSITY=.TRUE.
            ENDIF
            CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &        IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &        *9999)
            IF(ADATA(1).EQ.'Y') THEN
              IF(IOTYPE.NE.3) THEN
                CALC_SPARSITY=.TRUE.
              ENDIF
            ELSE
              WRITE(OP_STRING,'('' >>Remember to "read matrix;<FILE> '
     '          //'sparsity <ascii/binary>"'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            ENDIF

C CPB 6/11/95 Temporary work array allocation

C KAT 14Jan00: Calling CALC_SPARSE_FIT even for KTYP24=0 to check NZ_GK_M
C            IF(KTYP24.NE.0) THEN
            WORK_PTR=0
            CALL ALLOCATE_MEMORY(NYT(1,1,nx)*NYT(2,1,nx),1,CHARTYPE,
     '        WORK_PTR,MEM_INIT,ERROR,*9999)
          ENDIF !KTYP24.EQ.1
          CALL CALC_SPARSE_FIT(NISC_GKM,NISR_GKM,ISC_GK,ISR_GK,LGE,
     '      NYT(1,1,nx),NYT(2,1,nx),NBH,NEELEM,NPNE,NRLIST,NVHE,nx,NYNE,
     '      NYNP,KTYP24,%VAL(WORK_PTR),CALC_SPARSITY,ERROR,*9999)
          IF(KTYP24.NE.0) CALL FREE_MEMORY(WORK_PTR,ERROR,*9999)
C KAT 14Jan00
C          ENDIF
        ENDIF !nr
      ENDIF !KTYP8.NE.6 and .NOT. KTYP8.NE.8

C LKC 22-MAR-1998      IF(KTYP8.NE.8) THEN
      IF(KTYP8.LT.8) THEN
        IF(ITYP6(nr,nx).EQ.1) THEN !Linear
          ADEFLT(1)='N'
          FORMAT='(/$,'' Do you want to enter the coupling coefficients'
     '      //' [N]? '',A)'
          IF(IOTYPE.EQ.3) THEN
            IF(ENTERCOUPLING) THEN
              ADATA(1)='Y'
            ELSE
              ADATA(1)='N'
            ENDIF
          ENDIF
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            IF(ADATA(1).EQ.'Y'.OR.ADATA(1).EQ.'y') THEN
              ENTERCOUPLING=.TRUE.
              ERROR='>>Not implemented'
              GOTO 9999
            ELSE
              ENTERCOUPLING=.FALSE.
            ENDIF
          ENDIF
        ENDIF

C***    Initialise arrays

        IF(IOTYPE.NE.3) THEN
C***      Set all ny's in the region to being fixed (ie out of the fit)
          DO no_nynr=1,NYNR(0,0,1,nr) !loop over global variables
            ny=NYNR(no_nynr,0,1,nr) !is global variable number
            DO iy=1,NIYFIXM
              FIX(ny,iy)=.TRUE.
            ENDDO !iy
          ENDDO !no_nynr

C***     Set all ny's in the region and in the nj_fit to being free (ie
C***     in the fit)
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            DO njj=1,NUM_FIT(0)
              DO nhj=1,NUM_FIT(njj)
                nhx=NLH_FIT(nhj,3,njj)
                nh=NH_LOC(nhx,nx)
                DO nv=1,NVHP(nh,np,1,nr)
                  DO nk=1,NKH(nh,np,1,nr)
                    ny=NYNP(nk,nv,nh,np,0,1,nr)
                    DO iy=1,NIYFIXM
                      FIX(ny,iy)=.FALSE.
                    ENDDO !iy
                  ENDDO !nk
                ENDDO !nv
              ENDDO !nhj
            ENDDO !njj
          ENDDO !nonode (np)

C PM 14Aug02 : See IF-ENDIF below
C        IF(KTYP11.EQ.2) THEN !faces, fix unused versions
C          DO nonode=1,NPNODE(0,nr)
C            np=NPNODE(nonode,nr)
C            DO njj=1,NUM_FIT(0)
C              DO nhj=1,NUM_FIT(njj)
C                nhx=NLH_FIT(nhj,3,njj)
C                nh=NH_LOC(nhx,nx)
C                DO nv=1,NVHP(nh,np,1,nr)
C                  USED=.FALSE.
C                  DO noface=1,LN(0)
C                    nf=LN(noface)
C                    ne=NPF(6,nf)
C                    nef=NPF(8,nf)
C                    CALL CALC_FACE_INFORMATION_IND(NBJ(1,ne),
C     '                NBJF(1,nf),nef,NKJE(1,1,1,ne),NKEF,
C     '                NKJF,NNF,NPNE(1,1,ne),NPNF,nr,
C     '                NVJE(1,1,1,ne),NVJF,SE(1,1,ne),SF,ERROR,*9999)
C                    nbf=NBJF(nh,nf)
C                    DO nn=1,NNT(nbf)
C                      np2=NPNF(nn,nbf)
C                      nv2=NVJF(nn,nbf,nh)
C                      IF((nv2.EQ.nv).AND.(np2.EQ.np)) USED=.TRUE.
C                    ENDDO !nn
C                  ENDDO !nf
C                  IF(.NOT.USED) THEN
C                    DO nk=1,NKH(nh,np,1,nr)
C                      ny=NYNP(nk,nv,nh,np,0,1,nr)
C                      DO iy=1,NIYFIXM
C                        FIX(ny,iy)=.TRUE.
C                      ENDDO !iy
C                    ENDDO !nk
C                  ENDIF
C                ENDDO !nv
C              ENDDO !nhj
C            ENDDO !njj
C          ENDDO !nonode (np)
C        ENDIF

        IF(KTYP11.EQ.2) THEN !face fitting - initialise FIX
          DO njj=1,NUM_FIT(0)
            DO nhj=1,NUM_FIT(njj)
              nhx=NLH_FIT(nhj,3,njj)
              nh=NH_LOC(nhx,nx)
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                DO nv=1,NVHP(nh,np,1,nr)
                  DO nk=1,NKH(nh,np,1,nr)
                    ny=NYNP(nk,nv,nh,np,0,1,nr)
                    DO iy=1,NIYFIXM
                      FIX(ny,iy)=.TRUE.
                    ENDDO !iy
                  ENDDO !nk
                ENDDO !nv
              ENDDO !nonode
            ENDDO !nhj
          ENDDO !njj

C  Unfixing derivatives which are involved in the face fitting

          DO njj=1,NUM_FIT(0)
            DO nhj=1,NUM_FIT(njj)
              nhx=NLH_FIT(nhj,3,njj)
              nh=NH_LOC(nhx,nx)
              DO noface=1,LN(0)
                nf=LN(noface)
                ne=NPF(6,nf)
                nef=NPF(8,nf)
                CALL CALC_FACE_INFORMATION_IND(NBJ(1,ne),
     '            NBJF(1,nf),nef,NKJE(1,1,1,ne),NKEF,NKJF,NNF,
     '            NPNE(1,1,ne),NPNF,nr,NVJE(1,1,1,ne),NVJF,
     '            SE(1,1,ne),SF,ERROR,*9999)
                nbf=NBJF(nh,nf)
                DO nn=1,NNT(nbf)
                  np=NPNF(nn,nbf)
                  nv=NVJF(nn,nbf,nh)
                  DO nk=1,NKT(nn,nbf)
                    ny=NYNP(NKJF(nk,nn,nh),nv,nh,np,0,1,nr)
                    IF (FIX(ny,1)) THEN
                      DO iy=1,NIYFIXM
                        FIX(ny,iy)=.FALSE.
                      ENDDO !iy
                    ENDIF
                  ENDDO !nk
                ENDDO !nn
              ENDDO !noface
            ENDDO  !nhj
          ENDDO  !njj

          IF(DOP) THEN
            CALL WRITE_LINE(IODI,'FIX:',ERROR,*9999)
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
              DO njj=1,NUM_FIT(0)
                DO nhj=1,NUM_FIT(njj)
                  nj=NLH_FIT(nhj,1,njj)
                  nhx=NLH_FIT(nhj,3,njj)
                  nh=NH_LOC(nhx,nx)
                  DO nv=1,NVHP(nh,np,1,nr)
                    WRITE(OP_STRING(1),'(''np='',I5,'
     '                //''', nj='',I2,'', nv='',I2,'
     '                //''': '',99L1)') np,nj,nv,
     '                (FIX(NYNP(nk,nv,nh,np,0,1,nr),1),
     '                nk=1,NKH(nh,np,1,nr))
                    CALL WRITE_LINE(IODI,
     '                OP_STRING(1)(:LEN_TRIM(OP_STRING(1))),ERROR,*9999)
                  ENDDO !nv
                ENDDO !nhj
              ENDDO !njj
            ENDDO !nonode (np)
          ENDIF !DOP

        ENDIF ! KTYP11 == 2
        ENDIF ! iotype != 3

C ***   Enter constraints on fitting
        IF(KTYP8.NE.5) THEN
C         geom/field/fibre/sheet/patch/optimisation/signal

C 22/7/97 LC archived section : GMH 13/3/96 Add #s/name to fit stuff

          CONTINUE=.TRUE.
          nonode=0
 6100     IF(ENTERCOUPLING) THEN
            FORMAT=
     '        '(/$,'' Enter node #s/name to fix or to specify coupling'
     '        //' [EXIT]: '',I4)'
          ELSE
            FORMAT='(/$,'' Enter node #s/name to fix [EXIT]: '',I4)'
          ENDIF
          IF(IOTYPE.EQ.3) THEN
            nonode=nonode+1
            IF(nonode.LE.NPNODE(0,nr)) THEN
              np=NPNODE(nonode,nr)
              IDATA(1)=np
            ELSE
              IDATA(1)=0
            ENDIF
          ENDIF
 6500     CDATA(1)='NODES' !for use with group input
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,
     '      NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,
     '      ICHAR,IDATA,IZERO,0,NPT(nr),
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IDATA(1).NE.0) THEN !not default exit
            NPLIST(0)=IDATA(0)
            DO n=1,IDATA(0)
              NPLIST(n)=IDATA(n)
              np=IDATA(n)
              IF(.NOT.INLIST(np,NPNODE(1,nr),
     '          NPNODE(0,nr),N1)) THEN
                WRITE(OP_STRING,'('' Node '',I5,'' does not '
     '            //'belong to the current region'')') np
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                GOTO 6500
              ENDIF
            ENDDO !n

C           Define ft condition for first node in group
            np=NPLIST(1) !rest of group filled at end of njj loop
            DO njj=1,NUM_FIT(0)
              WRITE(CHAR1,'(I10)') njj
              CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
              DO nhj=1,NUM_FIT(njj)
                WRITE(CHAR2,'(I1)') nhj
                nhx=NLH_FIT(nhj,3,njj)
                nh=NH_LOC(nhx,nx)
                ADEFLT(1)='N'
                IF(NVHP(nh,np,1,nr).EQ.1.AND.NKH(nh,np,1,nr).EQ.1) THEN
                  FORMAT='($,'' Is variable '//CHAR2(1:1)//' of fit '
     '              //'variable '//CHAR1(IBEG1:IEND1)//' fixed '
     '              //'[N]?: '',A)'
                ELSE
                  FORMAT='($,'' Are any variables for variable '
     '              //CHAR2(1:1)//' of fit variable '
     '              //CHAR1(IBEG1:IEND1)//' fixed [N]?: '',A)'
                ENDIF
                IF(IOTYPE.EQ.3) THEN
                  ISFIXED=.FALSE.
                  DO nv=1,NVHP(nh,np,1,nr)
                    DO nk=1,NKH(nh,np,1,nr)
                      ny=NYNP(nk,nv,nh,np,0,1,nr)
                      IF(FIX(ny,1)) ISFIXED=.TRUE.
                    ENDDO !nk
                  ENDDO !nv
                  IF(ISFIXED) THEN
                    ADATA(1)='Y'
                  ELSE
                    ADATA(1)='N'
                  ENDIF
                ENDIF
                CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '            0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '            ERROR,*9999)
                IF(IOTYPE.NE.3) THEN
                  IF(ADATA(1).EQ.'Y'.OR.ADATA(1).EQ.'y') THEN
                    ISFIXED=.TRUE.
                  ELSE
                    ISFIXED=.FALSE.
                    DO nv=1,NVHP(nh,np,1,nr)
                      DO nk=1,NKH(nh,np,1,nr)
                        ny=NYNP(nk,nv,nh,np,0,1,nr)
                        FIX(ny,1)=.FALSE.
                      ENDDO !nk
                    ENDDO !nv
                  ENDIF
                ENDIF
                IF(ISFIXED) THEN
                  IF(NVHP(nh,np,1,nr).EQ.1.AND.NKH(nh,np,1,nr).EQ.1)
     '              THEN
C*** Don't need to ask a question as it has already been answered
                    ny=NYNP(1,1,nh,np,0,1,nr)
                    FIX(ny,1)=.TRUE.
                  ELSE
                    DO nv=1,NVHP(nh,np,1,nr) !loop over versions
                      IF(NVHP(nh,np,1,nr).GT.1) THEN
                        WRITE(CHAR3,'(I2)') nv
                        FORMAT='('' For version number '//CHAR3(1:2)
     '                    //':'')'
                        CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,1,0,NOQUES,
     '                    FILEIP,FORMAT,1,ADATA,ANO,CDATA,CDEFLT,
     '                    ICHAR,IDATA,IDEFLT,0,IMAX,LDATA,LDEFLT,
     '                    RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
                      ENDIF
                      DO nk=1,NKH(nh,np,1,nr) !loop over spatial derivs
                        ny=NYNP(nk,nv,nh,np,0,1,nr)
                        ADEFLT(1)='N'
                        IF(nk.EQ.1) THEN
                          FORMAT='($,'' Is variable '//CHAR2(1:1)
     '                      //' of fit variable '//CHAR1(IBEG1:IEND1)
     '                      //' fixed [N]?: '',A)'
                        ELSE IF(nk.GT.1) THEN
                          WRITE(CHAR3,'(I1)') nk
                          FORMAT='($,'' Is variable '//CHAR2(1:1)
     '                      //' of fit variable '
     '                      //CHAR1(IBEG1:IEND1)//' derivative '
     '                      //CHAR3(1:1)//' fixed [N]?: '',A)'
                        ENDIF
                        IF(IOTYPE.EQ.3) THEN
                          IF(FIX(ny,1)) THEN
                            ADATA(1)='Y'
                          ELSE
                            ADATA(1)='N'
                          ENDIF
                        ENDIF
                        CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,
     '                    FILEIP,FORMAT,1,ADATA,ANO,CDATA,CDEFLT,
     '                    ICHAR,IDATA,IDEFLT,0,IMAX,LDATA,LDEFLT,
     '                    RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
                        IF(IOTYPE.NE.3) THEN
                          IF(ADATA(1).EQ.'Y'.OR.ADATA(1).EQ.'y') THEN
                            FIX(ny,1)=.TRUE.
                          ELSE

C PM 14Aug02: To make user aware of incorrect derivative fixing in face-
C             fitting.
C                            FIX(ny,1)=.FALSE.

                            IF(KTYP11.EQ.1) THEN ! element fitting
                              FIX(ny,1)=.FALSE.
                            ELSE IF(KTYP11.EQ.2) THEN ! face fitting
                              IF(FIX(ny,1)) THEN
                                WRITE(OP_STRING,'(/'' Version '',i2,
     '                            '' of Derivative '',i2,'' of Fit '
     '                            //'variable '',i1,'' in Node '',i4,
     '                            '' cannot be unfixed !'')')
     '                            NPNY(2,ny,1),NPNY(1,ny,1),
     '                            NPNY(3,ny,1),NPNY(4,ny,1)
                                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                                WRITE(OP_STRING,'(/'' The variable is '
     '                            //'not in the Fit - Correct ipfit '
     '                            //'file !!!'')')
                                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                                FACESIG=.FALSE.
                              ENDIF
                            ENDIF

                          ENDIF
                        ENDIF
                        IF(FIX(ny,1).AND.ENTERCOUPLING) THEN
C*** variable is fixed - find out if it is coupled to another var.
                          ERROR='>> Not implemented'
                          GOTO 9999
                        ENDIF
                      ENDDO !nk
                    ENDDO !nv
                  ENDIF
                ENDIF !isfixed
              ENDDO !nhj
            ENDDO !njj
C           Apply fit conditions to rest of nodes group
C KAT 10Jan99: Only apply fix derivatives/versions fixed on node 1.
            DO n=2,NPLIST(0)
              np2=NPLIST(n)
              DO njj=1,NUM_FIT(0)
                DO nhj=1,NUM_FIT(njj)
                  nhx=NLH_FIT(nhj,3,njj)
                  nh=NH_LOC(nhx,nx)
                  DO nv=1,NVHP(nh,np,1,nr)
                    DO nk=1,NKH(nh,np,1,nr) !loop over spatial derivs
                      ny=NYNP(nk,nv,nh,np2,0,1,nr)
                      ny_first=NYNP(nk,nv,nh,np,0,1,nr)

C PM 18-Sep-02: To avoid incorrect derivative fixing in face fitting
C                      IF(ny.NE.0) FIX(ny,1)=FIX(ny_first,1)
                      IF (ny.NE.0) THEN
                        IF (KTYP11.EQ.1) THEN !element fitting
                          FIX(ny,1)=FIX(ny_first,1)
                        ELSEIF (KTYP11.EQ.2) THEN !face fitting
                          IF (FIX(ny,1).AND.(.NOT.FIX(ny_first,1))) THEN
                            WRITE(OP_STRING,'(/'' Version '',i2,
     '                        '' of Derivative '',i2,'' of Fit '
     '                        //'variable '',i1,'' in Node '',i4,
     '                        '' cannot be unfixed !'')') nv,nk,nh,np2
                            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                            WRITE(OP_STRING,'(/'' The variable is '
     '                        //'not in the Fit - Correct ipfit '
     '                        //'file !!!'')')
                            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                            FACESIG=.FALSE.
                          ELSE
                            FIX(ny,1)=FIX(ny_first,1)
                          ENDIF
                        ENDIF !element/face fitting
                      ENDIF !ny
                    ENDDO !nk
                  ENDDO !nv
                ENDDO !nhj
              ENDDO !njj
            ENDDO !n
            GO TO 6100 !for more nodes
          ENDIF !default exit

        ELSE IF(KTYP8.EQ.5) THEN !motion fitting with Fourier basis
        ENDIF
      ENDIF ! KTYP8.LT.8

C     Fix unused versions,
      IF((KTYP11.EQ.1).AND.(KTYP8.NE.3).AND.(KTYP8.NE.4).AND.
     '  (KTYP8.NE.9)) THEN

C     Not sure why this doesn't work for signal fitting
C     elements, inefficient but does the trick for now
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          DO njj=1,NUM_FIT(0)
            DO nhj=1,NUM_FIT(njj)
              nhx=NLH_FIT(nhj,3,njj)
              nh=NH_LOC(nhx,nx)
C VYW&MPN  10/03/2010/ Need to correct geometric basis functions
	      njj2=NLH_FIT(nhj,1,njj)
              DO nv=1,NVHP(nh,np,1,nr)
                USED=.FALSE.
                DO noface=1,NEELEM(0,nr)
                  ne=NEELEM(noface,nr)
C VYW&MPN  10/03/2010/ Need to correct geometric basis functions
                  nb=NBJ(njj2,ne)
C                  nb=NBJ(nh,ne)
                  DO nn=1,NNT(nb)
                    np2=NPNE(nn,nb,ne)
                    nv2=NVHE(nn,nb,nh,ne)
                    IF((nv2.EQ.nv).AND.(np2.EQ.np)) USED=.TRUE.
                  ENDDO !nn
                ENDDO !nf
                IF(.NOT.USED) THEN
                  DO nk=1,NKH(nh,np,1,nr)
                    ny=NYNP(nk,nv,nh,np,0,1,nr)
                    DO iy=1,NIYFIXM
                      FIX(ny,iy)=.TRUE.
                    ENDDO !iy
                  ENDDO !nk
                ENDIF
              ENDDO !nv
            ENDDO !nhj
          ENDDO !njj
        ENDDO !nonode (np)
      ELSE IF(KTYP11.EQ.2) THEN !faces
C PM 14Aug02 : I think this section is not necessary.
C        DO nonode=1,NPNODE(0,nr)
C          np=NPNODE(nonode,nr)
C          DO njj=1,NUM_FIT(0)
C            DO nhj=1,NUM_FIT(njj)
C              nhx=NLH_FIT(nhj,3,njj)
C              nh=NH_LOC(nhx,nx)
C              DO nv=1,NVHP(nh,np,1,nr)
C                USED=.FALSE.
C                DO noface=1,LN(0)
C                  nf=LN(noface)
C                  ne=NPF(6,nf)
C                  nef=NPF(8,nf)
C                  CALL CALC_FACE_INFORMATION_IND(NBJ(1,ne),
C     '              NBJF(1,nf),nef,NKJE(1,1,1,ne),NKEF,
C     '              NKJF,NNF,NPNE(1,1,ne),NPNF,nr,
C     '              NVJE(1,1,1,ne),NVJF,SE(1,1,ne),SF,ERROR,*9999)
C                  nbf=NBJF(nh,nf)
C                  DO nn=1,NNT(nbf)
C                    np2=NPNF(nn,nbf)
C                    nv2=NVJF(nn,nbf,nh)
C                    IF((nv2.EQ.nv).AND.(np2.EQ.np)) USED=.TRUE.
C                  ENDDO !nn
C                ENDDO !nf
C                IF(.NOT.USED) THEN
C                  DO nk=1,NKH(nh,np,1,nr)
C                    ny=NYNP(nk,nv,nh,np,0,1,nr)
C                    DO iy=1,NIYFIXM
C                      FIX(ny,iy)=.TRUE.
C                    ENDDO !iy
C                  ENDDO !nk
C                ENDIF
C              ENDDO !nv
C            ENDDO !nhj
C          ENDDO !njj
C        ENDDO !nonode (np)
      ENDIF

      IF(KTYP11.EQ.2) THEN
        CALL ASSERT(FACESIG,'Correct IPFIT file',ERROR,*9999)
      ENDIF

C *** Define solution parameters
C     For all cases other than optimisation and patch fitting
      IF((KTYP8.NE.6).AND.(KTYP8.NE.8).AND.(KTYP8.NE.9)) THEN
        CALL IPSOLU(nr,nx,ERROR,*9999)
      ENDIF

      CALL EXITS('IPFIT')
      RETURN
 9999 CALL ERRORS('IPFIT',ERROR)
      CALL EXITS('IPFIT')
      RETURN 1
      END


