      SUBROUTINE IPMAT3(CPBASIS,IBT,IDO,INP,NBJ,NEELEM,NELIST,
     '  NELIST2,NENP,NENQ,NNB,NPLIST,NPNE,NPNODE,NQET,NQNE,NQS,NQXI,nr,
     '  nx,NXI,POINTS,CE,CP,CQ,XE,YP,ALL_REGIONS,FIX,ERROR,*)

C#### Subroutine: IPMAT3
C###  Description:
C###    IPMAT3 inputs material parameters.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'b32.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'coro00.cmn'
      INCLUDE 'docu00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'iltot00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'titl30.cmn'
!     Parameter List
      INTEGER CPBASIS,IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NELIST(0:NEM),NELIST2(0:NEM),NENP(NPM,0:NEPM,0:NRM),
     '  NENQ(0:8,NQM),NNB(4,4,4,NBFM),NPLIST(0:NPM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NQET(NQSCM),NQNE(NEQM,NQEM),NQS(NEQM),
     '  NQXI(0:NIM,NQSCM),nr,nx,NXI(-NIM:NIM,0:NEIM,0:NEM),POINTS
      REAL*8 CE(NMM,NEM),CP(NMM,NPM),CQ(NMM,NQM),XE(NSM,NJM),
     '  YP(NYM,NIYM)
      CHARACTER ERROR*(*)
      LOGICAL ALL_REGIONS,FIX(NYM,NIYFIXM)
!     Local Variables
      INTEGER IB,IB2,IB5,IBEG,IBEG1,ICHAR,IE,IE2,IE5,IEND,IEND1,i,ii,ij,
     '  ik,il,INFO,K,nb,NBMIN,ne,neq,nh,NH_GROUP_TOT,nii,nij,nik,n,nk,
     '  nn,noelem,nonode,NOQUES,No_elem_vertex,No_transmural_elems,nm,
     '  NOTUSED,np,np1,nq,nqq,nqsc,ns
!     SMAR009 23/12/98 ,IDEFTYPE
      REAL*8 C3M,C3S,PXI,RSAMP,XI(3)
      CHARACTER CHAR*100,CHAR1*1,CHAR11*11,CHAR2*2,CHAR3*5,
     '  CHAR4*12,CHAR5*5,CHAR8*8,CHARNH*5,Elem_position*4
      LOGICAL ALLSET,FILEIP,NEEDFACES
!     Functions
      LOGICAL INLIST

      INTEGER TEMP_ITYP19,TEMP_ITYP3
      LOGICAL COUPLED_MODEL_HACK

C *** DPN 05 October 1999
C ***   Need another method of input distributed material parameters
C ***   for cellular based modelling - i.e. more than just ITYP2(nr,nx)
C ***   and ITYP3(nr,nx) dependent. Need to consider the type of model
C ***   (electrical, mechanical, coupled, etc..) as well as the
C ***   individual models (LR, JRW, HMT, etc...) e.g. want something
C ***   dependent on ITYP19(nr,nx) as well as ITYP2 and ITYP3.
C ***
C ***   For now, single cell mechanics doesn't use "fem def mate",
C ***   there is a hack for coupled models.
C ***

      CALL ENTERS('IPMAT3',*9999)

      CALL ASSERT(ITYP2(nr,nx).NE.0,
     '  '>>Solution type has not been defined',ERROR,*9999)
      CALL ASSERT(ITYP2(nr,nx).GT.2,'>>Solution type is incorrect',
     '  ERROR,*9999)

C *** DPN 05 October 1999 - hack for coupled models
      IF(ITYP19(nr,nx).EQ.7) THEN !coupled_model
        COUPLED_MODEL_HACK=.TRUE.
        TEMP_ITYP19=ITYP19(nr,nx)
        TEMP_ITYP3=ITYP3(nr,nx)
        !now pretend a straight electrical model
        ITYP19(nr,nx)=1
        !and a real electrical model at that!
        IF(ITYP3(nr,nx).EQ.1) THEN !LR-HMT
          ITYP3(nr,nx)=6
        ELSE !use a user defined model
          ITYP3(nr,nx)=10
        ENDIF
      ELSE
        COUPLED_MODEL_HACK=.FALSE.
      ENDIF

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      WRITE(CHAR1,'(I1)') ITYP5(nr,nx)
      WRITE(CHAR2,'(I2)') ITYP2(nr,nx)
      IF(CHAR2(1:1).EQ.' ') CHAR2(1:1)='0'
      WRITE(CHAR5,'(I1)') ITYP3(nr,nx)
      DOCFILE='FE'//CHAR1(1:1)//CHAR2(1:2)//CHAR5(1:1) !is docum.n file
      DOCLABEL='Parameter number'                      !is label string

      IF(ITYP5(nr,nx).EQ.1) THEN        !static analysis
        ILT(1,nr,nx)=ILTOT1(NJT,ITYP2(nr,nx),ITYP3(nr,nx))
      ELSE IF(ITYP5(nr,nx).EQ.2) THEN   !time integration
        ILT(1,nr,nx)=ILTOT2(NJT,ITYP2(nr,nx),ITYP3(nr,nx))
      ELSE IF(ITYP5(nr,nx).EQ.3) THEN   !modal analysis
        ILT(1,nr,nx)=ILTOT3(NJT,ITYP2(nr,nx),ITYP3(nr,nx))
      ELSE IF(ITYP5(nr,nx).EQ.4) THEN   !Quasi-static analysis
        ILT(1,nr,nx)=ILTOT4(NJT,ITYP2(nr,nx),ITYP3(nr,nx))
      ELSE IF(ITYP5(nr,nx).EQ.5) THEN   !Wavefront path analysis
        ILT(1,nr,nx)=ILTOT5(NJT,ITYP2(nr,nx),ITYP3(nr,nx))
      ELSE IF(ITYP5(nr,nx).EQ.6) THEN   !buckling analysis
        ILT(1,nr,nx)=ILTOT6(NJT,ITYP2(nr,nx),ITYP3(nr,nx))
      ENDIF
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' ILT(1,nr,nx)='',I2)') ILT(1,nr,nx)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF
c SMAR009 23/12/98      IDEFTYPE=1

c cpb 6/9/95 Isotropic and orthotropic materials cannot be handled
C until the title and total arrays have an extra index to indicate
C what the "trophy" of the material is

Cc cpb 4/9/95 Adding anisotropic materials for generalised Laplace
C      IMT(1)=0
C      IMT_OFFSET=0
C      IMT_EXTRA=0
C      IF((ITYP5(nr,nx).EQ.1.OR.ITYP5(nr,nx).EQ.4).AND.
C     '  ITYP2(nr,nx).EQ.3.AND.ITYP3(nr,nx).EQ.2) THEN
CC       Generalised Laplace
C        FORMAT='('' Specify whether the material is [1]: '''//
C     '    '/''  (1) Isotropic'''//
C     '    '/''  (2) Orthotropic '''//
C     '    '/$,''   '',I1)'
C        IF(IOTYPE.EQ.3) IDATA(1)=IMT(1)
C        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,2,
C     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C        IF(IOTYPE.NE.3) IMT(1)=IDATA(1)
C        IF(IMT(1).EQ.2) THEN
C          CALL ASSERT(NJ_LOC(njl_fibr,0,nr).GT.0,
C     '      '>>Fibre field has not been defined',
C     '      ERROR,*9999)
C          IMT_OFFSET=1
C          IMT_EXTRA=NJT-1
C        ENDIF
C      ENDIF

      CALL ASSERT(ILT(1,nr,nx).LE.NMM,'>>Increase NMM',ERROR,*9999)

C SGM 26Oct2000 grid-based Finite element also
C MLT 29Nov02 Grid FV also
      IF((ITYP2(nr,nx).EQ.9).AND.
     '  (ITYP4(nr,nx).EQ.4.OR.ITYP4(nr,nx).EQ.6.OR.
     '   ITYP4(nr,nx).EQ.7).AND.
     '  (ITYP5(nr,nx).EQ.2)) THEN
        ILT(0,nr,nx)=3
        nm=2
      ELSE
        ILT(0,nr,nx)=1
        nm=0
      ENDIF

C DMAL 18-JULY-2002 Adding loop over dependent variables
      IF(ITYP5(nr,nx).EQ.2) THEN !time-dependent
        IF(ITYP2(nr,nx).EQ.3) THEN !advection-diffusion
          NH_GROUP_TOT=KTYP3A(nx)
        ELSE
          NH_GROUP_TOT=1
        ENDIF
      ELSE
        NH_GROUP_TOT=1
      ENDIF

      DO nh=1,NH_GROUP_TOT !dependent variables
        DO il=ILT(0,nr,nx),ILT(1,nr,nx)
          nm=nm+1
          CALL ASSERT(nm.LE.NMM,'>>Increase NMM',ERROR,*9999)
          IF(ITYP5(nr,nx).EQ.1) THEN        !static analysis
            CHAR=TITL31(il,ITYP2(nr,nx),ITYP3(nr,nx))
          ELSE IF(ITYP5(nr,nx).EQ.2) THEN   !time integration
            CHAR=TITL32(il,ITYP2(nr,nx),ITYP3(nr,nx))
          ELSE IF(ITYP5(nr,nx).EQ.3) THEN   !modal analysis
            CHAR=TITL33(il,ITYP2(nr,nx),ITYP3(nr,nx))
          ELSE IF(ITYP5(nr,nx).EQ.4) THEN   !Quasi-static analysis
            CHAR=TITL31(il,ITYP2(nr,nx),ITYP3(nr,nx))
          ELSE IF(ITYP5(nr,nx).EQ.5) THEN   !Wavefront path analysis
            CHAR=TITL35(il,ITYP2(nr,nx),ITYP3(nr,nx))
          ELSE IF(ITYP5(nr,nx).EQ.6) THEN   !buckling analysis
            CHAR=TITL36(il,ITYP2(nr,nx),ITYP3(nr,nx))
          ENDIF
          CALL STRING_TRIM(CHAR,IBEG,IEND)
C DPN 23 November 1998 - change the printed default choice to match
C          call to ginout()
C          WRITE(CHAR1,'(I1)') IDEFTYPE
          IDEFLT(1)=POINTS
          WRITE(CHAR1,'(I1)') IDEFLT(1)
          WRITE(CHARNH,'(I5)') nh
          CALL STRING_TRIM(CHARNH,IBEG1,IEND1)
          IF(NH_GROUP_TOT.GT.1) THEN
            FORMAT='(/'' '//CHAR(IBEG:IEND)//' for dependent variable'//
     '        ' '//CHARNH(IBEG1:IEND1)//' is ['//CHAR1//']: '''//
     '        '/''   (1) Constant spatially       '''//
     '        '/''   (2) Piecewise constant (defined by elements)'''//
     '        '/''   (3) Piecewise linear   (defined by nodes)   '''//
     '        '/''   (4) Defined by grid points (CQ)       '''//
     '        '/''   (5) Defined elsewhere by Gauss pt array (YG)'''//
     '        '/$,'' '',I1)'
          ELSE
            FORMAT='(/'' '//CHAR(IBEG:IEND)//' is ['//CHAR1//']: '''//
     '        '/''   (1) Constant spatially       '''//
     '        '/''   (2) Piecewise constant (defined by elements)'''//
     '        '/''   (3) Piecewise linear   (defined by nodes)   '''//
     '        '/''   (4) Defined by grid points (CQ)       '''//
     '        '/''   (5) Defined elsewhere by Gauss pt array (YG)'''//
     '        '/$,'' '',I1)'
          ENDIF
          IF(IOTYPE.EQ.3) IDATA(1)=ILP(nm,1,nr,nx)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,5,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,nm,ERROR,*9999)
          IF(iotype.ne.3) ILP(nm,1,nr,nx)=IDATA(1)
!       SMAR009 23/12/98 IDEFTYPE=IDATA(1)

          IF(KTYP3_mate(nx).GT.1) THEN
            FORMAT='($,'' Is this parameter time-varying [Y]? '',A)'
            IF(IOTYPE.EQ.3) THEN
              IF(ILP(nm,1,nr,nx).GT.0) THEN
                ADATA(1)='N'
              ELSE IF(ILP(nm,1,nr,nx).LE.0) THEN
                ADATA(1)='Y'
              ENDIF
            ENDIF
            CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,AYES,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,0,ERROR,*9999)
            IF(iotype.ne.3) THEN
              IF(ADATA(1).EQ.'Y') ILP(nm,1,nr,nx)=-ILP(nm,1,nr,nx)
            ENDIF
          ENDIF

C KSB 12/07/04 - moved RDEFLT out here to ensure default defined for all problems
          RDEFLT(1)=RDMATE(nm,ITYP5(nr,nx),ITYP2(nr,nx),ITYP3(nr,nx))
          WRITE(CHAR11,'(E11.4)') RDEFLT(1)
          CALL STRING_TRIM(CHAR11,IB,IE)
          
          IF(ILP(nm,1,nr,nx).EQ.1) THEN !constant spatially
C            RDEFLT(1)=RDMATE(nm,ITYP5(nr,nx),ITYP2(nr,nx),ITYP3(nr,nx))
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
              WRITE(OP_STRING,'('' RDEFLT(1)='',E12.3)') RDEFLT(1)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
            ENDIF
C            WRITE(CHAR11,'(E11.4)') RDEFLT(1)
C            CALL STRING_TRIM(CHAR11,IB,IE)
            FORMAT='($,'' The value is ['//CHAR11(IB:IE)//']: '',E11.4)'
            IF(IOTYPE.EQ.3) RDATA(1)=CE(nm,NEELEM(1,nr))
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &        IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,il,ERROR,
     &        *9999)
            IF(iotype.ne.3) THEN
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                CE(nm,ne)=RDATA(1)
              ENDDO
C SGM 26Oct2000 grid-based Finite element also
C MLT 29Nov02 Grid FV also
              IF((ITYP4(nr,nx).EQ.4).OR.
     '           (ITYP4(nr,nx).EQ.6).OR.(ITYP4(nr,nx).EQ.7).OR.
     '           (ITYP4(nr,nx).EQ.3)) THEN 
C collocation or grid-based FE or finite differenes
C              DO nq=1,NQT
C                CQ(nm,nq)=RDATA(1)
C              ENDDO
                DO nq=NQR(1,nr),NQR(2,nr)
                  CQ(nm,nq)=RDATA(1)
                ENDDO
              ENDIF
            ENDIF

          ELSE IF(ILP(nm,1,nr,nx).EQ.2) THEN !defined by elements

C BGW 9/11/2001  Now can use element groups and multiple elements
            IF (IOTYPE.NE.3) THEN !Not write
              DO noelem=0,NEELEM(0,nr)
                NELIST2(noelem)=0
              ENDDO
            ENDIF
            DO noelem=0,NEELEM(0,nr)
              NELIST(noelem)=NEELEM(noelem,nr)
            ENDDO

 5000       IF (NELIST(1).NE.0) THEN !Default exit or exit
              IF (IOTYPE.EQ.3) THEN  !Writing out
                WRITE(CHAR11,'(E11.4)') 0
                DO noelem=1,NEELEM(0,nr) !Write nodes individually
                  NELIST(0)=1
                  ne=NEELEM(noelem,nr)
                  NELIST(1)=ne
                  FORMAT='($,'' Enter element #s/name [EXIT]: '',I5)'
                  CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,
     &              FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &              NELIST,IDEFLT,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,
     &              -RMAX,RMAX,il,ERROR,*9999)
                  RDATA(1)=CE(nm,ne)
                  CALL STRING_TRIM(CHAR11,IB,IE)
                  FORMAT=
     '              '($,'' Enter value for element(s) ['
     '              //CHAR11(IB:IE)//']: '',E12.5)'
                  CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     &              FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &              IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,
     &              -RMAX,RMAX,il,ERROR,*9999)
                    WRITE(CHAR11,'(E11.4)') RDATA(1)
                ENDDO
                NELIST(0)=1  !Exit
                NELIST(1)=0
                FORMAT='($,'' Enter element #s/name [EXIT]: '',I5)'
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,NELIST,
     '            IDEFLT,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,
     '            il,ERROR,*9999)
              ELSE !IOTYPE.ne.3

                CDATA(1)='ELEMENTS' !for use with group input
                IDEFLT(1)=0 !Default is exit
                FORMAT='($,'' Enter element #s/name [EXIT]: '',I5)'
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,NELIST,
     '            IDEFLT,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,
     '            il,ERROR,*9999)

                IF(NELIST(1).NE.0) THEN !not default exit
                  WRITE(CHAR11,'(E11.4)') RDEFLT(1)
                  CALL STRING_TRIM(CHAR11,IB,IE)
                  FORMAT='($,'' Enter value for element(s) ['
     '             //CHAR11(IB:IE)//']: '',E12.5)'

                  CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     &              FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &              IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,
     &              -RMAX,RMAX,il,ERROR,*9999)
                  RDEFLT(1)=RDATA(1)

           !Loop over requested elements
                  DO noelem=1,NELIST(0)
                    ne=NELIST(noelem)
                    CE(nm,ne)=RDATA(1)
                    NELIST2(ne)=1  !Element has been set

C SGM 26Oct2000 grid-based Finite element also
C PM 03May2001 finite difference added
C MLT 29Nov02 grid Finite Volumes added

                    IF(ITYP4(nr,nx).EQ.3.OR.ITYP4(nr,nx).EQ.4.OR.
     '                ITYP4(nr,nx).EQ.6.OR.ITYP4(nr,nx).EQ.7) THEN  !collocation/grid-based FE/FD/grid FV. 
                      DO nqq=1,NQET(NQS(ne))
                        nq=NQNE(ne,nqq)
                        CQ(nm,nq)=CE(nm,ne) !RDATA(1)
                      ENDDO
                    ENDIF

                  ENDDO !noelem
                ENDIF !not exit
                GOTO 5000
              ENDIF !nelist(1).NE.0
            ENDIF !IOTYPE=3

C           DO noelem=1,NEELEM(0,nr)
C           ne=NEELEM(noelem,nr)
C              IF(ne.EQ.1) RDEFLT(1)=RDMATE(il,ITYP5(nr,nx),
C     '   ITYP2(nr,nx),ITYP3(nr,nx))
C              IF(ne.GT.1) RDEFLT(1)=CE(il,ne-1)
C              WRITE(CHAR5,'(I5)') ne
C              CALL STRING_TRIM(CHAR5,IB5,IE5)
C              WRITE(CHAR11,'(E11.4)') RDEFLT(1)
C              CALL STRING_TRIM(CHAR11,IB,IE)
C              FORMAT='($,'' The value in element '//CHAR5(IB5:IE5)//
C     '   ' is ['//CHAR11(IB:IE)//']: '',E11.4)'
C              IF(IOTYPE.EQ.3) RDATA(1)=CE(il,ne)
C              CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '   ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
C     '   LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,il,ERROR,*9999)
C              IF(iotype.ne.3) THEN
C         CE(il,ne)=RDATA(1)
C              ENDIF

CC SGM 26Oct2000 grid-based Finite element also
CC PM 03May2001 finite difference added
C            IF(ITYP4(nr,nx).EQ.3.OR.ITYP4(nr,nx).EQ.4.OR.
C     ,        ITYP4(nr,nx).EQ.6) THEN  !collocation/grid-based FE/finite diff.
C              DO nqq=1,NQET(NQS(ne))
C                nq=NQNE(ne,nqq)
C                CQ(il,nq)=RDATA(1)
C              ENDDO

CC MLB 18 August 1997
CC old
CC     nb_extended=1
CC     DO WHILE (nb_extended.LE.NBT.AND.NBC(nb_extended).NE.7)
CC       nb_extended=nb_extended+1
CC     ENDDO
CC     CALL ASSERT((NBC(nb_extended).EQ.7),
CC '     'Extended basis function not defined',ERROR,*9999)
CC     DO ng=1,NGT(nb_extended)
CC       nq=NQGE(ng,ne,nb_extended)
CC       CQ(il,nq)=RDATA(1)
CC     ENDDO

C        ENDIF
C      ENDDO


C BGW 2/11/2001
            !Check to see if all elements are specified
            IF (IOTYPE.NE.3) THEN !Not write
              ALLSET=.TRUE.
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                IF ((NELIST2(ne).EQ.0).AND.ALLSET) THEN
                  ALLSET=.FALSE.
                  WRITE(CHAR5,'(I5)') ne
                ENDIF
              ENDDO
              CALL ASSERT(ALLSET,'>>Element '//CHAR5//' is not set',
     '          ERROR,*9999)
            ENDIF

          ELSE IF(ILP(nm,1,nr,nx).EQ.3) THEN !defined by nodes
            NEEDFACES=ITYP5(nr,nx).EQ.5 !Wavefront path analysis
     '        .AND.ITYP15(nr,nx).NE.0 !and upwind
            IF(NEEDFACES) THEN
              FORMAT='($,'' Enter element basis type number [1]: '',I1)'
              NBMIN=1
            ELSE
              IDEFLT(1)=CPBASIS
              WRITE(CHAR1,'(I1)') IDEFLT(1)
              FORMAT='($,'' Enter basis type number ['//CHAR1//']: '',
     '                 I1)'
              NBMIN=0
            ENDIF
            IF(IOTYPE.EQ.3) IDATA(1)=NMB(nm,1,nx)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &        NBMIN,NBFM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,0,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) NMB(nm,1,nx)=IDATA(1)

            IF(NEEDFACES) THEN
              FORMAT='($,'' Enter face basis type number [1]: '',I1)'
              IF(IOTYPE.EQ.3) IDATA(1)=NMB(nm,2,nx)
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,
     &          NBFM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,0,ERROR,*9999)
              IF(IOTYPE.NE.3) NMB(nm,2,nx)=IDATA(1)

C             Prompt for nodal parameters allowing for node group input
              IF(IOTYPE.EQ.3) THEN
                nonode=0 !init node loop for writing params
              ELSE IF(IOTYPE.NE.3) THEN
C             init CP so can check that each node has been set later
                DO nonode=1,NPNODE(0,nr)
                  np=NPNODE(nonode,nr)
                  CP(nm,np)=-RMAX
                ENDDO !nonode (np)
              ENDIF
 7000         FORMAT='($,'' Enter node #s/name [EXIT]: '',I5)'
              IF(IOTYPE.EQ.3) THEN
              nonode=nonode+1
                IF(nonode.LE.NPNODE(0,nr)) THEN
                  np=NPNODE(nonode,nr)
                  IDATA(1)=np
                ELSE
                  IDATA(1)=0 !to terminate node loop
                ENDIF
                IDATA(0)=1 !write out one node at a time
              ENDIF !iotype=3
 7100         CDATA(1)='NODES' !for use with group input
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,
     '          NPT(nr),LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     '          *9999)
              IF(IDATA(1).NE.0) THEN !not default exit
                NPLIST(0)=IDATA(0)
                DO n=1,IDATA(0)
                  NPLIST(n)=IDATA(n)
                  np=IDATA(n)
                  IF(.NOT.INLIST(np,NPNODE(1,nr),NPNODE(0,nr),NOTUSED))
     '              THEN
                    WRITE(OP_STRING,'('' Node '',I5,'' does not '
     '                //'belong to the current region'')') np
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                    GOTO 7100
                  ENDIF
                ENDDO !n
C               Define parameter for first node in group
                np=NPLIST(1) !rest of group filled at end
                FORMAT='($,'' The value is [0]: '',D11.4)'
                IF(IOTYPE.EQ.3) RDATA(1)=CP(nm,np)
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IDEFLT,IMIN,IMAX,
     '            LDATA,LDEFLT,RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
                IF(IOTYPE.NE.3) THEN
                  CP(nm,np)=RDATA(1)
C                 Set parameters for rest of nodes in group
                  DO n=2,NPLIST(0)
                    np1=NPLIST(n)
                    CP(nm,np1)=CP(nm,np)
                  ENDDO !n
                ENDIF
                GO TO 7000 !for more nodes
              ENDIF !idata(1).NE.0

              IF(IOTYPE.NE.3) THEN
C        check that CP has been set for each node
                ALLSET=.TRUE.
                DO nonode=1,NPNODE(0,nr)
                  np=NPNODE(nonode,nr)
                  IF(CP(nm,np).EQ.-RMAX) THEN
                    ALLSET=.FALSE.
                    CP(nm,np)=0.0d0
                  ENDIF
                ENDDO !nonode (np)
                IF(.NOT.ALLSET) THEN
                  WRITE(OP_STRING,'('' >>WARNING: Parameter has been '
     '              //'assigned zero for nodes not specified'')')
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDIF
              ENDIF

C SGM 30 Oct 2000 ignore grid-based Finite element also
            ELSE IF(.FALSE.)
     '          THEN !transverse activation wavespeed only
C DMAL 22-OCT-2002 Could not work out what this 'elseif' condition was
C to achieve and it wasn't precise enough  to ensure it doesn't effect
C other situation. So Karl suggested I set it to false and if someone
C comes forward they can explain what it is trying to do.
C            ELSE IF((NKT(0,NMB(nm,1,nx)).EQ.2).AND.
C     '          (ITYP4(nr,nx).NE.4.AND.ITYP4(nr,nx).NE.6))
C     '          THEN !transverse activation wavespeed only
              WRITE(OP_STRING,'('' Note: This should only be used'
     '          //' for transverse wavespeed'')')
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              CALL NENXI(NBJ,NEELEM,NENP,NNB,NPNE,NXI,ERROR,*9999) !to return elements surrounding ne in NXI
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
!             Find whether one or two transmural elements at this node
                No_transmural_elems=0
                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  nb=NBJ(1,ne)
                  DO nn=1,NNT(nb)
                    IF(NPNE(nn,nb,ne).eq.np) THEN
                      IF(NXI(3,1,ne).GT.0) THEN
                        Elem_position='endo'
                        No_transmural_elems=2
                        No_elem_vertex=nn
                      ELSE IF(NXI(-3,1,ne).GT.0) THEN
                        Elem_position='epic'
                        No_transmural_elems=2
                        No_elem_vertex=nn
                      ELSE
                        IF(no_transmural_elems.ne.2) THEN
                          No_transmural_elems=1
                          No_elem_vertex=nn
                        ENDIF
                      ENDIF
                      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
                        WRITE(OP_STRING,
     '                    '('' No_transmural_elems='',I1)')
     '                    No_transmural_elems
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                        WRITE(OP_STRING,'('' Elem_position='',A)')
     '                    Elem_position
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                        WRITE(OP_STRING,'('' No_elem_vertex='',I2)')
     '                    No_elem_vertex
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
                      ENDIF
                      GO TO 200
                    ENDIF
                  ENDDO
 200              CONTINUE
                ENDDO
                C3M=DMAX1(CE(1,1)/3.0d0,CE(2,1)/2.0d0) !centre point transv. veloc.
                C3S=DMAX1(CE(1,1)/2.0d0,CE(2,1)/1.5d0) !outside transverse velocity
                IF(No_transmural_elems.EQ.1) THEN
                  CP(nm,2*np-1)=C3S !velocity
                  IF(No_elem_vertex.LE.4) THEN      !endocardial node
                    CP(nm,2*np)= 4.0d0*(C3M-C3S) !velocity derivative
                  ELSE IF(No_elem_vertex.GT.4) THEN !epicardial node
                    CP(nm,2*np)=-4.0d0*(C3M-C3S) !velocity derivative
                  ENDIF
                ELSE IF(No_transmural_elems.EQ.2) THEN
                  IF(Elem_position.EQ.'endo') THEN
                    IF(No_elem_vertex.LE.4) THEN      !endocardial node
                      CP(nm,2*np-1)=C3S   !velocity
                      CP(nm,2*np)= 2.0d0*(C3M-C3S) !velocity derivative
                    ELSE IF(No_elem_vertex.GT.4) THEN !midwall node
                      CP(nm,2*np-1)=C3M    !velocity
                      CP(nm,2*np)= 0.0d0    !velocity derivative
                    ENDIF
                  ELSE IF(Elem_position.EQ.'epic') THEN
                    IF(No_elem_vertex.LE.4) THEN      !midwall node
                      CP(nm,2*np-1)=C3M   !velocity
                      CP(nm,2*np)= 0.0d0    !velocity derivative
                    ELSE IF(No_elem_vertex.GT.4) THEN !epicardial node
                      CP(nm,2*np-1)=C3S   !velocity
                      CP(nm,2*np)=-2.0d0*(C3M-C3S) !velocity derivative
                    ENDIF
                  ENDIF
                ENDIF
              ENDDO
C DMAL 14-NOV-2002 This code doen't seem to do anything so I have commented
C it out with CAPITAL "C", the lowercase "c" were already there before
C I commented this section of code.
C            ELSE IF(NMB(nm,1,nx).EQ.0) THEN  !defined at nodes & not interpolated
C              nonode=0
C 6100         FORMAT='($,'' Enter node number [EXIT]: '',I5)'
C              IF(IOTYPE.EQ.3) THEN
C                nonode=nonode+1
c          IDATA(1)=NODE_SOURCE(nonode)
C              ENDIF
C              CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
C     '          1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,NPM,
C     '          LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C              IF(IDATA(1).ne.0) THEN
C                IF(iotype.ne.3) THEN
C                  nonode=nonode+1
c            NODE_SOURCE(nonode)=IDATA(1)
C                ENDIF
C                GO TO 6100
C              ENDIF
            ELSE IF(NMB(nm,1,nx).GE.0) THEN !defined at nodes & interpolated
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                IF(nonode.EQ.1) RDEFLT(1)=RDMATE(nm,ITYP5(nr,nx),
     '          ITYP2(nr,nx),ITYP3(nr,nx))
                IF(nonode.GT.1) RDEFLT(1)=CP(nm,
     '          NPNODE(nonode-1,nr))
                WRITE(CHAR5,'(I5)') np
                CALL STRING_TRIM(CHAR5,IB5,IE5)
                WRITE(CHAR11,'(E11.4)') RDEFLT(1)
                CALL STRING_TRIM(CHAR11,IB,IE)
                FORMAT='($,'' The value at node '//CHAR5(IB5:IE5)//
     '            ' is ['//CHAR11(IB:IE)//']: '',E11.4)'
                IF(IOTYPE.EQ.3) RDATA(1)=CP(nm,np)
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '            IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,il,
     '            ERROR,*9999)
                IF(iotype.ne.3) CP(nm,np)=RDATA(1)
              ENDDO
            ENDIF

C SGM 26Oct2000 grid-based Finite element also
C MLT 29Nov02 Grid Finite Volume also
            IF(ITYP4(nr,nx).EQ.4.OR.ITYP4(nr,nx).EQ.6.OR.
     '         ITYP4(nr,nx).EQ.7) THEN !Grid based
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                nb=NMB(nm,1,nx)
                ns=0
                DO nn=1,NNT(nb)
                  DO nk=1,NKT(nn,nb)
                    ns=ns+1
C!!! Assuming linear basis!!!! (as above i.e. doesn't ask for derivatives)
                    IF(nk.EQ.1) THEN
                      XE(ns,1)=CP(nm,NPNE(nn,nb,ne))
                    ELSE
                      XE(ns,1)=0.0d0
                    ENDIF
                  ENDDO !nk
                ENDDO !nn
                nqsc=NQS(ne)
                II=MAX(1,NQXI(1,nqsc))
                IJ=1
                IK=1
                IF(NQXI(0,nqsc).GT.1) IJ=MAX(1,NQXI(2,nqsc))
                IF(NQXI(0,nqsc).GT.2) IK=MAX(1,NQXI(3,nqsc))
                DO i=1,3
                  XI(i)=0.0d0
                ENDDO !i
! Loop over the grid points in each element
                DO nik=1,IK
                  DO nij=1,IJ
                    DO nii=1,II
                      neq=nii+((nij-1)*NQXI(1,nqsc))
                      IF(NQXI(0,nqsc).GT.1) neq=neq+((nik-1)*
     '                  NQXI(1,nqsc)*NQXI(2,nqsc))
                      nq=NQNE(ne,neq)
! Evaluate each grid point only once
                      IF(NENQ(1,nq).EQ.ne) THEN
! Local xi coordinates of grid point nq in element ne
                        IF(II.NE.1) XI(1)=DBLE(nii-1)/DBLE(II-1)
                        IF(IJ.NE.1) XI(2)=DBLE(nij-1)/DBLE(IJ-1)
                        IF(IK.NE.1) XI(3)=DBLE(nik-1)/DBLE(IK-1)

                        CQ(nm,nq)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                    INP(1,1,nb),nb,1,XI,XE(1,1))

                      ENDIF

                    ENDDO !nii
                  ENDDO !nij
                ENDDO !nik
              ENDDO !noelem
            ENDIF

          ELSE IF(ILP(nm,1,nr,nx).EQ.4) THEN !defined by grid points
            CALL ASSERT(NQT.GT.0,'>>No grid points defined',ERROR,*9999)

            IF(IOTYPE.NE.3) THEN
              FORMAT='($,'' Enter grid values by element [Y]? '',A)'
              ADEFLT(1)='Y'
              ADATA(1)='Y'
              CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,1,
     &          IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &          *9999)

              IF(ADATA(1).EQ.'Y') THEN
 802            FORMAT='($,'' Enter element #/s [EXIT]: '',I5)'
                CDATA(1)='ELEMENTS' !for use with group input
                CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IZERO,0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,
     '            INFO,ERROR,*9999)

                IF(IDATA(1).NE.0) THEN !not default exit
                  NELIST(0)=IDATA(0)
                  DO ne=1,NELIST(0)
                    NELIST(ne)=IDATA(ne)
                  ENDDO

                  FORMAT='($,'' Is the value constant within '
     '              //'element(s) [Y]? '',A)'
                  ADEFLT(1)='Y'
                  ADATA(1)='Y'
                  CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,
     &              FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &              IDATA,IZERO,1,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,
     &              RMAX,INFO,ERROR,*9999)

                  IF(ADATA(1).EQ.'Y') THEN
                    RDEFLT(1)=RDMATE(nm,ITYP5(nr,nx),ITYP2(nr,nx),
     '              ITYP3(nr,nx))
                    WRITE(CHAR11,'(E11.4)') RDEFLT(1)
                    CALL STRING_TRIM(CHAR11,IB,IE)
                    FORMAT='($,'' The value in the selected elements'
     '                //' is ['//CHAR11(IB:IE)//']: '',E11.4)'
                    CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     &                FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &                IDATA,IZERO,1,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,
     &                -RMAX,RMAX,INFO,ERROR,*9999)

                    DO noelem=1,NELIST(0)
                      ne=NELIST(noelem)
                      DO nqq=1,NQET(NQS(ne))
                        nq=NQNE(ne,nqq)
                        CQ(nm,nq)=RDATA(1)
                      ENDDO
                    ENDDO
                    GO TO 802
                  ELSE
                    DO noelem=1,NELIST(0)
                      ne=NELIST(noelem)
                      DO nqq=1,NQET(NQS(ne))
                        nq=NQNE(ne,nqq)
                        IF(nqq.EQ.1) THEN
                          RDEFLT(1)=RDMATE(nm,ITYP5(nr,nx),ITYP2(nr,nx),
     '                      ITYP3(nr,nx))
                        ELSE
                          RDEFLT(1)=CQ(nm,NQNE(ne,nqq-1))
                        ENDIF
                        WRITE(CHAR5,'(I5)') nq
                        CALL STRING_TRIM(CHAR5,IB5,IE5)
                        WRITE(CHAR11,'(E11.4)') RDEFLT(1)
                        CALL STRING_TRIM(CHAR11,IB,IE)
                        FORMAT='($,'' The value for grid '
     '                    //CHAR5(IB5:IE5)//
     '                    ' is ['//CHAR11(IB:IE)//']: '',E11.4)'
                        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     '                    FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,
     '                    ICHAR,IDATA,IZERO,1,IMAX,LDATA,LDEFLT,RDATA,
     '                    RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
                        CQ(nm,nq)=RDATA(1)
                      ENDDO
                    ENDDO
                    GOTO 802
                  ENDIF
                ENDIF
              ELSE
                DO nq=1,NQT
                  IF(nq.EQ.1) THEN
                    RDEFLT(1)=RDMATE(nm,ITYP5(nr,nx),ITYP2(nr,nx),
     '                ITYP3(nr,nx))
                  ELSE
                    RDEFLT(1)=CQ(nm,nq-1)
                  ENDIF
                  WRITE(CHAR8,'(I8)') nq
                  CALL STRING_TRIM(CHAR8,IB5,IE5)
                  WRITE(CHAR11,'(E11.4)') RDEFLT(1)
                  CALL STRING_TRIM(CHAR11,IB,IE)
                  FORMAT='($,'' The value for grid '//CHAR8(IB5:IE5)//
     '              ' is ['//CHAR11(IB:IE)//']: '',E11.4)'
                  CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     &              FILEIP,FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,
     &              IDATA,IZERO,1,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '              RMAX,INFO,ERROR,*9999)
                  CQ(nm,nq)=RDATA(1)
                ENDDO
              ENDIF
            ELSE IF(IOTYPE.EQ.3) THEN
              FORMAT='($,'' Enter grid values by element [Y]? '',A)'
              ADEFLT(1)='N'
              ADATA(1)='N'
              CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,1,
     &          IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &          *9999)
              DO nq=1,NQT
                RDEFLT(1)=CQ(nm,nq)
                WRITE(CHAR8,'(I8)') nq
                CALL STRING_TRIM(CHAR8,IB5,IE5)
                WRITE(CHAR11,'(E11.4)') RDEFLT(1)
                CALL STRING_TRIM(CHAR11,IB,IE)
                FORMAT='($,'' The value for grid '//CHAR8(IB5:IE5)//
     '            ' is ['//CHAR11(IB:IE)//']: '',E11.4)'
                RDATA(1)=CQ(nm,nq)
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '            FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     '            IZERO,1,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,
     '            RMAX,INFO,ERROR,*9999)
              ENDDO
            ENDIF

          ELSE IF(ILP(nm,1,nr,nx).EQ.5) THEN !defined elsewhere by elements

          ENDIF
       ENDDO ! il
      ENDDO !nh

      IF(KTYP14.GT.0) THEN
        FORMAT='($,'' Specify material parameter to be incremented'
     '    //' [1]: '',I2)'
        IF(IOTYPE.EQ.3) IDATA(1)=KTYP14
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NMM,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,0,ERROR,*9999)
        IF(iotype.ne.3) KTYP14=IDATA(1)
      ENDIF

c      IF(ITYP5(nr,nx).EQ.2.AND.ITYP2(nr,nx).EQ.9.AND.
c     '  ITYP3(nr,nx).EQ.6) THEN
cC *** Luo-Rudy ionic current model
c        Nao=140.d0
c        Nai=18.d0
c        Ki=145.d0
c        DO nq=1,NQT
c          Ko=CQ(13,nq) !Extracellular potassium concentration
c          CQ(17,nq) = 0.282d0*DSQRT(Ko/5.4d0)
c          CQ(18,nq) = 0.6047d0*DSQRT(Ko/5.4d0)
c          CQ(19,nq) = -26.7d0*DLOG((Ki+0.01833d0*Nai)/
c     '      (Ko+0.01833d0*Nao))
c          CQ(20,nq) = -26.7d0*DLOG(Ki/Ko)
c        ENDDO
c      ENDIF !Luo-Rudy ionic current model


      IF((ITYP4(nr,nx).EQ.3).AND.(ITYP3(nr,nx).EQ.1)) THEN !coronaries

        FORMAT='('' Enter the 6 numerator values for Ra'')'
        CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,0,0,
     '    NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,
     '    ICHAR,IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
     '    RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)


        NUM_BOX(1,1)=0.00051137d0
        NUM_BOX(2,1)=-7.3677d-05
        NUM_BOX(3,1)=0.0035752d0
        NUM_BOX(4,1)=-0.0021702d0
        NUM_BOX(5,1)=4.2085d-06
        NUM_BOX(6,1)=1.0d0


        DO nm=1,6
          RDEFLT(1)=NUM_BOX(nm,1)
          WRITE(CHAR3,'(I1)') nm
          WRITE(CHAR4,'(E12.5)') NUM_BOX(nm,1)
C LKC 6-DEC-2000 Zero length string
C          FORMAT='($,'' numerator value '//CHAR3(1:1)//''//
C     '      ' is ['//CHAR4(1:12)//']: '',G25.17)'
          FORMAT='($,'' numerator value '//CHAR3(1:1)//
     '      ' is ['//CHAR4(1:12)//']: '',G25.17)'
          IF(IOTYPE.EQ.3) THEN
            RDATA(1)=NUM_BOX(nm,1)
          ENDIF
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,
     '      NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,
     '      ICHAR,IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
     '      RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            NUM_BOX(nm,1)=RDATA(1)
          ENDIF
        ENDDO



        FORMAT='('' Enter the 6 denominator values for Ra'')'
        CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,0,0,
     '    NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,
     '    ICHAR,IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
     '    RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)


        DENOM_BOX(1,1)=3.1192d-07
        DENOM_BOX(2,1)=-1.0483d-08
        DENOM_BOX(3,1)=-3.17d-05
        DENOM_BOX(4,1)=-5.0848d-07
        DENOM_BOX(5,1)=2.2473d-09
        DENOM_BOX(6,1)=0.0010545d0

        DO nm=1,6
          RDEFLT(1)=DENOM_BOX(nm,1)
          WRITE(CHAR3,'(I1)') nm
          WRITE(CHAR4,'(E12.5)') DENOM_BOX(nm,1)
C LKC 6-DEC-2000 Zero length string
C          FORMAT='($,'' numerator value '//CHAR3(1:1)//''//
C     '      ' is ['//CHAR4(1:12)//']: '',G25.17)'
          FORMAT='($,'' numerator value '//CHAR3(1:1)//
     '      ' is ['//CHAR4(1:12)//']: '',G25.17)'
          IF(IOTYPE.EQ.3) THEN
            RDATA(1)=DENOM_BOX(nm,1)
          ENDIF
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,
     '      NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,
     '      ICHAR,IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
     '      RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            DENOM_BOX(nm,1)=RDATA(1)
          ENDIF
        ENDDO

        FORMAT='('' Enter the 6 numerator values for Rc'')'
        CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,0,0,
     '    NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,
     '    ICHAR,IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
     '    RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)

        NUM_BOX(1,2)=1.0d0
        NUM_BOX(2,2)=-0.11236d0
        NUM_BOX(3,2)=-0.083349d0
        NUM_BOX(4,2)=0.0055689d0
        NUM_BOX(5,2)=0.0082348d0
        NUM_BOX(6,2)=0.0069464d0

        DO nm=1,6
          RDEFLT(1)=NUM_BOX(nm,2)
          WRITE(CHAR3,'(I1)') nm
          WRITE(CHAR4,'(E12.5)') NUM_BOX(nm,2)
C LKC 6-DEC-2000 Zero length string
C          FORMAT='($,'' numerator value '//CHAR3(1:1)//''//
C     '      ' is ['//CHAR4(1:12)//']: '',G25.17)'
          FORMAT='($,'' numerator value '//CHAR3(1:1)//
     '      ' is ['//CHAR4(1:12)//']: '',G25.17)'
          IF(IOTYPE.EQ.3) THEN
            RDATA(1)=NUM_BOX(nm,2)
          ENDIF
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,
     '      NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,
     '      ICHAR,IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
     '      RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            NUM_BOX(nm,2)=RDATA(1)
          ENDIF
        ENDDO



        FORMAT='('' Enter the 6 denominator values for Rc'')'
        CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,0,0,
     '    NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,
     '    ICHAR,IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
     '    RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)

        DENOM_BOX(1,2)=0.00043352d0
        DENOM_BOX(2,2)=1.9563d-05
        DENOM_BOX(3,2)=3.2255d-05
        DENOM_BOX(4,2)=2.4147d-06
        DENOM_BOX(5,2)=9.9149d-06
        DENOM_BOX(6,2)=1.2144d-05

        DO nm=1,6
          RDEFLT(1)=DENOM_BOX(nm,2)
          WRITE(CHAR3,'(I1)') nm
          WRITE(CHAR4,'(E12.5)') DENOM_BOX(nm,2)
C LKC 6-DEC-2000 Zero length string
C          FORMAT='($,'' numerator value '//CHAR3(1:1)//''//
C     '      ' is ['//CHAR4(1:12)//']: '',G25.17)'
          FORMAT='($,'' numerator value '//CHAR3(1:1)//
     '      ' is ['//CHAR4(1:12)//']: '',G25.17)'
          IF(IOTYPE.EQ.3) THEN
            RDATA(1)=DENOM_BOX(nm,2)
          ENDIF
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,
     '      NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,
     '      ICHAR,IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
     '      RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            DENOM_BOX(nm,2)=RDATA(1)
          ENDIF
        ENDDO

        FORMAT='('' Enter the 6 numerator values for Rv'')'
        CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,0,0,
     '    NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,
     '    ICHAR,IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
     '    RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)

        NUM_BOX(1,3)=0.0034779d0
        NUM_BOX(2,3)=-0.00066418d0
        NUM_BOX(3,3)=-0.0047076d0
        NUM_BOX(4,3)=-0.0125d0
        NUM_BOX(5,3)=7.4907d-05
        NUM_BOX(6,3)=1.0d0

        DO nm=1,6
          RDEFLT(1)=NUM_BOX(nm,3)
          WRITE(CHAR3,'(I1)') nm
          WRITE(CHAR4,'(E12.5)') NUM_BOX(nm,3)
C LKC 6-DEC-2000 Zero length string
C          FORMAT='($,'' numerator value '//CHAR3(1:1)//''//
C     '      ' is ['//CHAR4(1:12)//']: '',G25.17)'
          FORMAT='($,'' numerator value '//CHAR3(1:1)//
     '      ' is ['//CHAR4(1:12)//']: '',G25.17)'
          IF(IOTYPE.EQ.3) THEN
            RDATA(1)=NUM_BOX(nm,3)
          ENDIF
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,
     '      NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,
     '      ICHAR,IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
     '      RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            NUM_BOX(nm,3)=RDATA(1)
          ENDIF
        ENDDO


        FORMAT='('' Enter the 6 denominator values for Rv'')'
        CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,0,0,
     '    NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,
     '    ICHAR,IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
     '    RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)

        DENOM_BOX(1,3)=7.6075d-06
        DENOM_BOX(2,3)=3.276d-07
        DENOM_BOX(3,3)=-0.00046939d0
        DENOM_BOX(4,3)=-2.8796d-05
        DENOM_BOX(5,3)=2.1416d-07
        DENOM_BOX(6,3)=0.0027147d0

        DO nm=1,6
          RDEFLT(1)=DENOM_BOX(nm,3)
          WRITE(CHAR3,'(I1)') nm
          WRITE(CHAR4,'(E12.5)') DENOM_BOX(nm,3)
C LKC 6-DEC-2000 Zero length string
C          FORMAT='($,'' numerator value '//CHAR3(1:1)//''//
C     '      ' is ['//CHAR4(1:12)//']: '',G25.17)'
          FORMAT='($,'' numerator value '//CHAR3(1:1)//
     '      ' is ['//CHAR4(1:12)//']: '',G25.17)'
          IF(IOTYPE.EQ.3) THEN
            RDATA(1)=DENOM_BOX(nm,3)
          ENDIF
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,
     '      NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,
     '      ICHAR,IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
     '      RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            DENOM_BOX(nm,3)=RDATA(1)
          ENDIF
        ENDDO

        FORMAT='('' Enter the 6 numerator values for C1'')'
        CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,0,0,
     '    NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,
     '    ICHAR,IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
     '    RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)

        NUM_BOX(1,4)=3.0104d-05
        NUM_BOX(2,4)=3.2133d-06
        NUM_BOX(3,4)=-1.2763d-05
        NUM_BOX(4,4)=-8.6688d-07
        NUM_BOX(5,4)=6.9337d-08
        NUM_BOX(6,4)=3.3129d-06

        DO nm=1,6
          RDEFLT(1)=NUM_BOX(nm,4)
          WRITE(CHAR3,'(I1)') nm
          WRITE(CHAR4,'(E12.5)') NUM_BOX(nm,4)
C LKC 6-DEC-2000 Zero length string
C          FORMAT='($,'' numerator value '//CHAR3(1:1)//''//
C     '      ' is ['//CHAR4(1:12)//']: '',G25.17)'
          FORMAT='($,'' numerator value '//CHAR3(1:1)//
     '      ' is ['//CHAR4(1:12)//']: '',G25.17)'
          IF(IOTYPE.EQ.3) THEN
            RDATA(1)=NUM_BOX(nm,4)
          ENDIF
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,
     '      NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,
     '      ICHAR,IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
     '      RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            NUM_BOX(nm,4)=RDATA(1)
          ENDIF
        ENDDO


        FORMAT='('' Enter the 6 denominator values for C1'')'
        CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,0,0,
     '    NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,
     '    ICHAR,IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
     '    RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)


        DENOM_BOX(1,4)=1.0d0
        DENOM_BOX(2,4)=0.069804d0
        DENOM_BOX(3,4)=-0.74405d0
        DENOM_BOX(4,4)=-0.023532d0
        DENOM_BOX(5,4)=0.0015796d0
        DENOM_BOX(6,4)=0.16945d0

        DO nm=1,6
          RDEFLT(1)=DENOM_BOX(nm,4)
          WRITE(CHAR3,'(I1)') nm
          WRITE(CHAR4,'(E12.5)') DENOM_BOX(nm,4)
C LKC 6-DEC-2000 Zero length string
C          FORMAT='($,'' numerator value '//CHAR3(1:1)//''//
C     '      ' is ['//CHAR4(1:12)//']: '',G25.17)'
          FORMAT='($,'' numerator value '//CHAR3(1:1)//
     '      ' is ['//CHAR4(1:12)//']: '',G25.17)'
          IF(IOTYPE.EQ.3) THEN
            RDATA(1)=DENOM_BOX(nm,4)
          ENDIF
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,
     '      NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,
     '      ICHAR,IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
     '      RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            DENOM_BOX(nm,4)=RDATA(1)
          ENDIF
        ENDDO


        FORMAT='('' Enter the 6 numerator values for C2'')'
        CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,0,0,
     '    NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,
     '    ICHAR,IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
     '    RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)

        NUM_BOX(1,5)=0.00020578d0
        NUM_BOX(2,5)=1.3136d-05
        NUM_BOX(3,5)=-4.428d-05
        NUM_BOX(4,5)=-2.3643d-06
        NUM_BOX(5,5)=5.2479d-06
        NUM_BOX(6,5)=8.2425d-06

        DO nm=1,6
          RDEFLT(1)=NUM_BOX(nm,5)
          WRITE(CHAR3,'(I1)') nm
          WRITE(CHAR4,'(E12.5)') NUM_BOX(nm,5)
C LKC 6-DEC-2000 Zero length string
C          FORMAT='($,'' numerator value '//CHAR3(1:1)//''//
C     '      ' is ['//CHAR4(1:12)//']: '',G25.17)'
          FORMAT='($,'' numerator value '//CHAR3(1:1)//
     '      ' is ['//CHAR4(1:12)//']: '',G25.17)'
          IF(IOTYPE.EQ.3) THEN
            RDATA(1)=NUM_BOX(nm,5)
          ENDIF
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,
     '      NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,
     '      ICHAR,IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
     '      RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            NUM_BOX(nm,5)=RDATA(1)
          ENDIF
        ENDDO


        FORMAT='('' Enter the 6 denominator values for C2'')'
        CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,0,0,
     '    NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,
     '    ICHAR,IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
     '    RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)

        DENOM_BOX(1,5)=1.0d0
        DENOM_BOX(2,5)=0.001729d0
        DENOM_BOX(3,5)=-0.36117d0
        DENOM_BOX(4,5)=-0.0062574d0
        DENOM_BOX(5,5)=0.0443551d0
        DENOM_BOX(6,5)=0.05926d0

        DO nm=1,6
          RDEFLT(1)=DENOM_BOX(nm,5)
          WRITE(CHAR3,'(I1)') nm
          WRITE(CHAR4,'(E12.5)') DENOM_BOX(nm,5)
C LKC 6-DEC-2000 Zero length string
C          FORMAT='($,'' numerator value '//CHAR3(1:1)//''//
C     '      ' is ['//CHAR4(1:12)//']: '',G25.17)'
          FORMAT='($,'' numerator value '//CHAR3(1:1)//
     '      ' is ['//CHAR4(1:12)//']: '',G25.17)'
          IF(IOTYPE.EQ.3) THEN
            RDATA(1)=DENOM_BOX(nm,5)
          ENDIF
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,
     '      NOQUES,FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,
     '      ICHAR,IDATA,IDEFLT,IMIN,IMAX,LDATA,LDEFLT,
     '      RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            DENOM_BOX(nm,5)=RDATA(1)
          ENDIF
        ENDDO

      ENDIF !coronaries

      IF(ITYP5(nr,nx).EQ.3.AND.ITYP2(nr,nx).EQ.3) THEN
C ***   Modal; vocal tract
        IF(ILP(1,1,nr,nx).EQ.2) NQT=NET(nr)
        IF(ILP(1,1,nr,nx).EQ.3) NQT=NPT(nr)
        FORMAT='($,'' Specify whether area parameters are fixed (T)'//
     '    ' or free [F]:'')'
        CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,0,ERROR,*9999)
        DO nq=1,NQT
          WRITE(CHAR5,'(I5)') nq
          CALL STRING_TRIM(CHAR5,IB5,IE5)
C         CHAR=TITL10(ILP(1,1,nr,nx)-1)
          CALL STRING_TRIM(CHAR,IB,IE)
          FORMAT='($,'' '//CHAR(IB:IE)//' '//CHAR5(IB5:IE5)//
     '      ' [T]: '',L1)'
          IF(IOTYPE.EQ.3) LDATA(1)=FIX(nq,4)
          CALL GINOUT(IOTYPE,IPLOGI,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LTRUE,RDATA,RDEFLT,RMIN,RMAX,0,ERROR,*9999)
          IF(iotype.ne.3) FIX(nq,4)=LDATA(1)
        ENDDO

        FORMAT='($,'' Specify whether vocal tract length is'//
     '    ' is fixed (t) or free [F]: '',L1)'
        IF(IOTYPE.EQ.3) LDATA(1)=FIX(NQT+1,4)
        CALL GINOUT(IOTYPE,IPLOGI,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LFALSE,RDATA,RDEFLT,RMIN,RMAX,0,ERROR,*9999)
        IF(iotype.ne.3) FIX(NQT+1,4)=LDATA(1)
        WRITE(CHAR2,'(I2)') KTYP17
        CALL STRING_TRIM(CHAR2,IB2,IE2)
        DO K=1,KTYP17
          IF(K.EQ.1) RDEFLT(1)=0.d0
          IF(K.GT.1) RDEFLT(1)=YP(K-1,4)
          RSAMP=RDEFLT(1)
          WRITE(CHAR11,'(E11.4)') RSAMP
          CALL STRING_TRIM(CHAR11,IB,IE)
          FORMAT='('' Enter '//CHAR2(IB2:IE2)//' observed '''//
     '      ' eigenvalue for the'''//
     '      '/$,'' least squares fit ['//CHAR11(IB:IE)//']: '',11.4)'
          IF(IOTYPE.EQ.3) RDATA(1)=YP(K,4)
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,0,ERROR,*9999)
          IF(iotype.ne.3) YP(K,4)=RDATA(1)
        ENDDO
      ENDIF

C *** DPN 05 October 1999 - hack for coupled models
      IF(COUPLED_MODEL_HACK) THEN
        ITYP19(nr,nx)=TEMP_ITYP19
        ITYP3(nr,nx)=TEMP_ITYP3
      ENDIF
C DMAL 06-JUN-2002 Time-dependent advection-diffusion for multiple species
C Needed to be changed for CPCG to calculate the material parameters at
C gauss point for all dependent variables.
      IF(ITYP5(nr,nx).EQ.2.AND.ITYP2(nr,nx).EQ.3) THEN
        IF(KTYP3A(nx).GT.1) THEN
          ILT(1,nr,nx)=ILT(1,nr,nx)*NH_LOC(0,nx)
        ENDIF
      ENDIF
C GBS  8-Jul-2003 Another hack for grid-based FEM/FVM with Laplace equation
      IF(ITYP5(nr,nx).EQ.1.AND.ITYP2(nr,nx).EQ.3.AND.
     '  (ITYP4(nr,nx).EQ.6.OR.ITYP4(nr,nx).EQ.7)) THEN
        DO nq=1,NQT
          CQ(1,nq)=1.0d0
          CQ(2,nq)=1.0d0
          CQ(3,nq)=1.0d0
          CQ(6,nq)=1.0d0
          CQ(7,nq)=1.0d0
          CQ(8,nq)=1.0d0
        ENDDO
      ENDIF


      IF(FILEIP.AND..NOT.ALL_REGIONS) CALL CLOSEF(IFILE,ERROR,*9999)
      CALL EXITS('IPMAT3')
      RETURN
 9999 CALL ERRORS('IPMAT3',ERROR)
      IF(FILEIP) CLOSE(UNIT=IFILE)
      CALL EXITS('IPMAT3')
      RETURN 1
      END



