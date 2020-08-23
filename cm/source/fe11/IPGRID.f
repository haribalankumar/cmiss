      SUBROUTINE IPGRID(TYPE,NBJ,NEELEM,NELIST,NENP,NENQ,NLL,NNB,
     &  NPNE,NQET,NQNE,NQS,NQSCNB,NQXI,NRLIST,NWQ,NXI,NXQ,
     &  NLATNE_PTR,NLATNQ_PTR,NLATPNQ_PTR,NQNLAT_PTR,DL,XIQ,ERROR,*)

C#### Subroutine: IPGRID
C###  Description:
C###    IPGRID inputs parameters for grid generation. The grids are
C###    tied directly to finite elements, the grid points are generated
C###    at equally spaced xi locations. There must be a separate grid
C###    scheme for each number of xi coordinates, ie. if there are
C###    1d and 2d elements then there must be at least 2 schemes,
C###    one for the 1d elements and one for the 2d elements. A different
C###    number of collocation points may be specified in each xi
C###    direction but care must be taken to ensure that the connectivity
C###    is consistent along element boundaries.
C###    The coronary option generates the finite difference grid
C###    positions and connectivity arrays within the finite element
C###    network mesh such that the delta x spacing between grid points
C###    is consistent between elements. Because of the solution method
C###    for coronary problems three grid points are placed at
C###    each bifurcation position.
C**** Written by Martin Buist, June 1997

C#### Variable: NQSCT
C###  Type: INTEGER
C###  Set_up: IPGRID
C###  Description:
C###    NQSCT is the total number of grid schemes.

C#### Variable: nqsc
C###  Type: INTEGER
C###  Description:
C###    nqsc is the grid scheme loop variable.

      IMPLICIT NONE

      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     &  NELIST(0:NEM),NENP(NPM,0:NEPM,0:NRM),NENQ(0:8,NQM),
     &  NLL(12,NEM),NNB(4,4,4,NBFM),NPNE(NNM,NBFM,NEM),NQET(NQSCM),
     &  NQNE(NEQM,NQEM),NQS(NEQM),NQSCNB(NQSCM),NQXI(0:NIM,NQSCM),
     &  NRLIST(0:NRM),NWQ(8,0:NQM,NAM),NXI(-NIM:NIM,0:NEIM,0:NEM),
     &  NXQ(-NIM:NIM,0:4,0:NQM,NAM)
      INTEGER*4 NLATNE_PTR,NLATNQ_PTR,NLATPNQ_PTR,NQNLAT_PTR
      REAL*8 DL(3,NLM),XIQ(NIM,NQM)
      CHARACTER ERROR*(*),TYPE*8

!     Local Variables
      INTEGER COR_SCHEME(1820),i,IBEG,ICHAR,IEND,IJ,IK,INFO,ISUM,mq,
     &  mq1,mq2,mq3,mq4,mq5,mq6,n,na,nb,ne,ni,nii,nij,nik,NITB,
     &  nl,noelem,NOQUES,nq,nq_ele,nq_adj,nqe,nqe1,nqe2,nqe3,
     &  NQFINISH(3),nqnew,nqq,nqsc,NQSTART(3),nr,nrr,nstep,
     &  NUM_NQ,SCHEME
      REAL*8 DELTA_X,LENGTH,
     &  SMALL_LENGTH                      
      CHARACTER CHAR1*4
      LOGICAL ALLSET,BACK_BIFUR,FOR_BIFUR,CHECKNXQ,FILEIP,
     &  GENER(-3:3),INTERNAL
!     External functions
      !INTEGER IDIGITS

      CALL ENTERS('IPGRID',*9999)


! Initialise variables
C      CALL ASSERT(NAM.GE.2,'>>Increase NAM to be >= 2',ERROR,*9999)
C Use NWQ(4..6,nq,1) for temp storage.

C!!! KAT: Is CALL_GRID set to false when rereading region 1?
C         What are UP_N.N. used for?  Shouldn't they be reset?
C 17-NOV-2004 JHC If-Statement is made so that when grid is 
C defined in region 2, variables are not initialised again
C      IF ((.NOT.UP_NQNP).AND.(.NOT.UP_NENQ)) THEN
      IF(.NOT.CALL_GRID)THEN
        UP_NQNP=.TRUE.
        UP_NENQ=.TRUE.
C KAT: Is this necessary?
        DO nq=0,NQM
          DO ni=-NIM,NIM
            DO i=0,4
              DO na=1,NAM
                NXQ(ni,i,nq,na)=0
              ENDDO !na
            ENDDO !i
          ENDDO !ni
        ENDDO !nq

C LKC 11-JAN-2011 I believe we need to initialise all NQS values
C (or the whole array) -> remove the initalisation for a specific
C region below ... 
        DO ne=1,NEQM !initialise NQS(ne)
          NQS(ne)=0
        ENDDO
      
      ENDIF !call grid

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999
      N_DISCRET=1

! Define the collocation schemes to be used
      IF(TYPE(1:7).EQ.'DEFAULT') THEN
        IDEFLT(1)=1
        IF(IOTYPE.EQ.3) IDATA(1)=NQSCT
        FORMAT='($,'' The number of grid schemes is [1]: '',I5)'
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     &    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NQSCM,
     &    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) NQSCT=IDATA(1)

        CALL ASSERT(NQSCT.LE.NQSCM,'>>Increase NQSCM',
     &    ERROR,*9999)
        CHECKNXQ=.FALSE.

        DO nqsc=1,NQSCT
          IF(IOTYPE.NE.3) THEN
            DO ni=1,NIM
C DAH 28/5/02 Now initialising all ni for each grid scheme as 1
C             NQXI(ni,nqsc)=0
              NQXI(ni,nqsc)=1
            ENDDO !ni
          ENDIF

          WRITE(CHAR1,'(I2)') nqsc
          CALL STRING_TRIM(CHAR1,IBEG,IEND)

          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=NQSCNB(nqsc)
          FORMAT='($,'' Enter the interpolating basis function for '
     &      // 'scheme '//CHAR1(IBEG:IEND)//' [1]: '',I2)'
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NBM,
     &      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NQSCNB(nqsc)=IDATA(1)

C KAT 2007-01-30: Could this be determined from the basis function?
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=NQXI(0,nqsc)
          FORMAT='($,'' Enter the number of xi coordinates for '
     &      // 'scheme '//CHAR1(IBEG:IEND)//' [1]: '',I1)'
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NIM,
     &      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NQXI(0,nqsc)=IDATA(1)

          DO ni=1,NQXI(0,nqsc)
            WRITE(CHAR1,'(I2)') ni
            CALL STRING_TRIM(CHAR1,IBEG,IEND)
            IDEFLT(1)=1
            IF(IOTYPE.EQ.3) IDATA(1)=NQXI(ni,nqsc)
            FORMAT='($,'' Enter the number of collocation points in '
     &        // 'xi '//CHAR1(IBEG:IEND)//' direction [1]: '',I1)'
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,
     &        NQEM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) NQXI(ni,nqsc)=IDATA(1)
C PJH 1Jun99 temporarily commented to check use of de grid;d
C            CALL ASSERT(NQXI(ni,nqsc).GT.2,
C     '        '>>Need at least 3 grid points',ERROR,*9999)
          ENDDO !ni


      
          IDEFLT(1)=1
          FORMAT='($,'' Enter the number of grid levels [1]: '',I2)'
          IF(IOTYPE.EQ.3) IDATA(1)=NMGT
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,11,
     &      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NMGT=IDATA(1)
          CALL ASSERT(NAM.GE.NMGT,'>>NAM too small',ERROR,*9999)

          IF(NQXI(0,nqsc).EQ.1) NQET(nqsc)=NQXI(1,nqsc)
          IF(NQXI(0,nqsc).GT.1) THEN
            IF(NQXI(1,nqsc).NE.NQXI(2,nqsc)) CHECKNXQ=.TRUE.
            NQET(nqsc)=NQXI(1,nqsc)*NQXI(2,nqsc)
          ENDIF
          IF(NQXI(0,nqsc).GT.2) THEN
            IF(NQXI(1,nqsc).NE.NQXI(3,nqsc)) CHECKNXQ=.TRUE.
            IF(NQXI(2,nqsc).NE.NQXI(3,nqsc)) CHECKNXQ=.TRUE.
            NQET(nqsc)=NQXI(1,nqsc)*NQXI(2,nqsc)*NQXI(3,nqsc)
          ENDIF
          WRITE(OP_STRING(1),
     &      '(''>>Increase NQEM >= '',I12)') NQET(nqsc)
          CALL STRING_TRIM(OP_STRING(1),IBEG,IEND)
          CALL ASSERT(NQET(nqsc).LE.NQEM,OP_STRING(1)(IBEG:IEND),ERROR,
     &      *9999)
        ENDDO !nqsc
        
        FORMAT='($,'' Use lattice based grid points'
     &    //' [N]? '',A)'      
        IF(IOTYPE.EQ.3) THEN
          IF(USE_LAT.EQ.1) THEN
            ADATA(1)='Y'
          ELSE
            ADATA(1)='N'
          ENDIF
        ENDIF
        CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     &    ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     &    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          IF(ADATA(1).EQ.'Y') THEN
            USE_LAT=1
            WRITE(OP_STRING,'(''WARNING: The lattice grid '
     &        //'scheme needs further testing'')')
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)            
          ELSE
            USE_LAT=0
          ENDIF
        ENDIF

        
! Warning if different number of points in each direction
        IF(CHECKNXQ) THEN
          WRITE(OP_STRING,'('' >>WARNING Check points on element '
     &      // 'boundaries are consistent'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

! Assign a collocation scheme to each specified element
        IF(IOTYPE.NE.3) THEN

C LKC 11-JAN-2011 - the whole NQS array is now initialised above.
C
C          DO nrr=1,NRLIST(0)
C            nr=NRLIST(nrr)
C            DO noelem=1,NEELEM(0,nr)
C              ne=NEELEM(noelem,nr)
C              NQS(ne)=0
C            ENDDO !ne
C          ENDDO !nr

 802      FORMAT='($,'' Enter element #s/name [EXIT]: '',I5)'
          CDATA(1)='ELEMENTS' !for use with group input
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,NET(0),
     &      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)

          IF(IDATA(1).NE.0) THEN !not default exit
            NELIST(0)=IDATA(0)
            DO n=1,IDATA(0)
              NELIST(n)=IDATA(n)
            ENDDO !n

            ne=NELIST(1) !rest of group is filled in afterwards
            FORMAT='($,'' Enter the grid scheme number for '
     &        // 'element(s) [1]: '',I2)'
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,
     &        NQSCM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &        *9999)
            NQS(ne)=IDATA(1)
            CALL ASSERT(NIT(NBJ(1,ne)).EQ.NQXI(0,NQS(ne)),
     &        '>>Scheme and element have different # xi dirns',
     &        ERROR,*9999)

            DO n=2,NELIST(0) !fill in rest of group
              ne=NELIST(n)
              NQS(ne)=NQS(NELIST(1))
            ENDDO
            GO TO 802

C PJH 1Jun99 add new default for de grid;d
          ELSE !default exit
            IF(NEELEM(0,0).EQ.1.AND.NQXI(0,1).EQ.1) THEN !one element & one grid pt
              NQS(1)=1 !is default grid scheme for one element
            ENDIF
          ENDIF !IDATA(1)

        ELSE IF(IOTYPE.EQ.3) THEN
          DO nrr=1,NRLIST(0)
            nr=NRLIST(nrr)
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              FORMAT='($,'' Enter element #s/name [EXIT]: '',I5)'
              IDATA(1)=ne
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,
     &          NET(0),LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &          *9999)
              FORMAT='($,'' Enter the grid scheme number for '
     &          // 'element(s) [1]: '',I2)'
              IDATA(1)=NQS(ne)
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,
     &          NQSCM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &          *9999)
            ENDDO !noelem
          ENDDO !nrr
          FORMAT='($,'' Enter element #s/name [EXIT]: '',I5)'
          IDATA(1)=0
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,NET(0),
     &      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        ENDIF !IOTYPE

        ALLSET=.TRUE.
        DO nrr=1,NRLIST(0)
          nr=NRLIST(nrr)
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            IF(NQS(ne).EQ.0) ALLSET=.FALSE.
          ENDDO !ne
        ENDDO !nr
        CALL ASSERT(ALLSET,'>>Must set grid scheme for all elements'
     &    //' in set regions',ERROR,*9999)

      ELSE IF(TYPE(1:8).EQ.'CORONARY') THEN
C Creating and assigning grid point schemes for coronary mesh
        DELTA_X=10.0d0
        FORMAT='('' Enter value for delta x [10.0]:'''//
     &    '/$,''    '',F4.2)'
        RDEFLT(1)=10.0d0
        IF(IOTYPE.EQ.3) RDATA(1)=DELTA_X
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     &    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     &    LDATA,LDEFLT,RDATA,RDEFLT,0.0D0,20.0D0,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) DELTA_X=RDATA(1)

C reads in the delta x spacing between grid points

C PM 27-NOV-01: New question added.
        FORMAT='('' Specify whether mesh discretisation is '//
     &    'with [1] : '''//
     &    '/''   (1) Uniform xi intervals   '''//
     &    '/''   (2) Uniform arc-length intervals  '''//
     &    '/$,''    '',I1)'
        IDEFLT(1)=1
        IF(IOTYPE.EQ.3) IDATA(1)=N_DISCRET
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &    1,2,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) N_DISCRET=IDATA(1)

        NQSCT=0
        SMALL_LENGTH=DELTA_X
        DO i=1,1820
          COR_SCHEME(i)=0
        ENDDO

        DO nrr=1,NRLIST(0) ! loops over regions
          nr=NRLIST(nrr)
          NJT=NJ_LOC(NJL_GEOM,0,nr)
          DO noelem=1,NEELEM(0,nr) ! loops over elements
C determines  length of each cornary element
            ne=NEELEM(noelem,nr)
            nb=NBJ(1,ne)
C PM 26-JUL-01 : length is now calculated using DL
C           np1=NPNE(1,nb,ne)
C           np2=NPNE(NNT(nb),nb,ne)

C           DO nj=1,NJT
C           NODE1(nj)=XP(1,1,nj,np1)
C           NODE2(nj)=XP(1,1,nj,np2)
C           ENDDO !nj

C           CALL COORD(ITYP10(nr),1,NODE1,CART_NODE1,ERROR,*9999)
C           CALL COORD(ITYP10(nr),1,NODE2,CART_NODE2,ERROR,*9999)

C calculates the cartesian position of the nodes in each element

C            LENGTH=0.0d0
C            DO nj=1,NJT
C              LENGTH=LENGTH+((CART_NODE1(nj)-CART_NODE2(nj))**2.0d0)
C            ENDDO !nj
C            LENGTH=LENGTH**0.5d0

            nl=NLL(1,ne)
            LENGTH=DL(3,nl)

            NUM_NQ=NINT(LENGTH/DELTA_X)
C calcultes the number of grid points in element ne
            NUM_NQ=NUM_NQ+1 !one more grid point than x steps
            IF(NUM_NQ.EQ.3) THEN !records smallest delta x step
              IF((LENGTH/2.0d0).LT.SMALL_LENGTH) THEN
                SMALL_LENGTH=LENGTH/2.0d0
              ENDIF
            ENDIF
            IF((NXI(-1,0,ne).LE.1).AND.(NXI(1,0,ne).LE.1)) THEN
c element does not have a bifurcation at either end and thus 2
c grid points is the minimum needed
              IF(NUM_NQ.LT.2) THEN
                NUM_NQ=2
                IF(LENGTH.LT.SMALL_LENGTH) THEN
                  SMALL_LENGTH=LENGTH
                ENDIF
              ENDIF
            ELSE IF(NENP(NPNE(1,nb,ne),0,nr).GT.2.AND.
     &          NENP(NPNE(2,nb,ne),0,nr).GT.2) THEN
c bifurcation at each end of element, needs 4 grid points minimum
              IF(NUM_NQ.LT.4) THEN
                NUM_NQ=4
                IF(LENGTH.LT.SMALL_LENGTH) THEN
                  SMALL_LENGTH=LENGTH
                ENDIF
              ENDIF
            ELSE !bifurcation at 1 end
c element does have a bifurcation at least one end and thus 3
c grid points is the minium needed
              IF(NUM_NQ.LT.3) THEN
                NUM_NQ=3
                IF((LENGTH/2.0d0).LT.SMALL_LENGTH) THEN
                  SMALL_LENGTH=LENGTH/2.0d0
                ENDIF
              ENDIF
            ENDIF !no bifurcation/bifurcation
            IF(COR_SCHEME(NUM_NQ).GT.0) THEN
              NQS(ne)=COR_SCHEME(NUM_NQ)
            ELSE ! create new scheme
              NQSCT=NQSCT+1
              CALL ASSERT(NQSCT.LE.NQSCM,'>>Increase NQSCM',
     &          ERROR,*9999)
              NQSCNB(NQSCT)=nb
              NQXI(0,NQSCT)=1
              NQXI(NQXI(0,NQSCT),NQSCT)=NUM_NQ
              NQET(NQSCT)=NQXI(1,NQSCT)
              WRITE(OP_STRING(1),
     &          '(''>>Increase NQEM >= '',I12)') NQET(NQSCT)
              CALL STRING_TRIM(OP_STRING(1),IBEG,IEND)
              CALL ASSERT(NQET(NQSCT).LE.NQEM,OP_STRING(1)(IBEG:IEND),
     &          ERROR,*9999)
              NQS(ne)=NQSCT
              COR_SCHEME(NUM_NQ)=NQSCT
            ENDIF !scheme
          ENDDO !noelem
        ENDDO !nrr
      ENDIF !TYPE

      IF(USE_LAT.EQ.0) THEN !use the standard method of grid generation
        
! Find element adjacency
        CALL NENXI(NBJ,NEELEM,NENP,NNB,NPNE,NXI,ERROR,*9999)
        
        nqnew=NQT+1 ! Generate the required number of grid points in each element

        WRITE(OP_STRING(1),
     &    '(''>>Increase NQM >= '',I12)') nqnew
        CALL STRING_TRIM(OP_STRING(1),IBEG,IEND)
        CALL ASSERT(nqnew.LE.NQM,OP_STRING(1)(IBEG:IEND),ERROR,*9999)
        DO nrr=1,NRLIST(0)
          nr=NRLIST(nrr)
          NQR(1,nr)=nqnew
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            
! Set variables for the current element
            SCHEME=NQS(ne)
            DO ni=1,NQXI(0,SCHEME)
              NQSTART(ni)=1
              NQFINISH(ni)=NQXI(ni,SCHEME)
              GENER(ni)=.FALSE.
              GENER(-ni)=.FALSE.
            ENDDO !ni
            
! If there is only one xi direction
            IF(NQXI(0,SCHEME).EQ.1) THEN
              nb=NQSCNB(SCHEME)
             IF(KTYP3B.NE.2) THEN !not Gauss point scheme
! Element in -xi 1
              IF(NXI(-1,0,ne).GT.0) THEN
! Tests for a coronary bifurcation in the -xi 1 direction
                IF(TYPE(1:8).EQ.'CORONARY'
     &            .AND.(NXI(-1,0,ne).GE.2)) THEN
                  BACK_BIFUR=.TRUE.
                ELSE
                  BACK_BIFUR=.FALSE.
                ENDIF
                DO nii=1,NXI(-1,0,ne)
                  IF(NXI(-1,nii,ne).LT.ne) THEN
                    GENER(-1)=.TRUE.
! If we are looking back onto the end of the previous element
                    IF(NPNE(1,nb,ne).EQ.NPNE(NNT(nb),nb,
     &                NXI(-1,nii,ne))) THEN
                      IF(.NOT.BACK_BIFUR) THEN
                        NQNE(ne,1)=NQNE(NXI(-1,nii,ne),
     &                    NQET(NQS(NXI(-1,nii,ne))))
                      ELSE
                        nq_ele=nqnew
                        nq_adj=NQNE(NXI(-1,nii,ne),
     &                    NQET(NQS(NXI(-1,nii,ne))))
                        WRITE(OP_STRING(1),
     &                    '(''>>Increase NQM >= '',I12)') MAX(nq_ele,
     &                    nq_adj)
                        CALL STRING_TRIM(OP_STRING(1),IBEG,IEND)
                        CALL ASSERT(nq_ele.LE.NQM.AND.nq_adj.LE.NQM,
     &                    OP_STRING(1)(IBEG:IEND),ERROR,*9999)
                        NXQ(-1,0,nq_ele,1)=NXQ(-1,0,nq_ele,1)+1
                        NXQ(-1,NXQ(-1,0,nq_ele,1),nq_ele,1)=nq_adj
                        NXQ(1,0,nq_adj,1)=NXQ(1,0,nq_adj,1)+1
                        NXQ(1,NXQ(1,0,nq_adj,1),nq_adj,1)=nq_ele
                      ENDIF
! If we are looking back onto the start of the previous element
                    ELSE IF(NPNE(1,nb,ne).EQ.NPNE(1,nb,NXI(-1,nii,ne)))
     &                  THEN
                      IF(.NOT.BACK_BIFUR) THEN
                        NQNE(ne,1)=NQNE(NXI(-1,nii,ne),1)
                      ELSE
                        nq_ele=nqnew
                        nq_adj=NQNE(NXI(-1,nii,ne),1)
                        WRITE(OP_STRING(1),
     &                    '(''>>Increase NQM >= '',I12)') MAX(nq_ele,
     &                    nq_adj)
                        CALL STRING_TRIM(OP_STRING(1),IBEG,IEND)
                        CALL ASSERT(nq_ele.LE.NQM.AND.nq_adj.LE.NQM,
     &                    OP_STRING(1)(IBEG:IEND),ERROR,*9999)
                        NXQ(-1,0,nq_ele,1)=NXQ(-1,0,nq_ele,1)+1
                        NXQ(-1,NXQ(-1,0,nq_ele,1),nq_ele,1)=nq_adj
                        NXQ(-1,0,nq_adj,1)=NXQ(-1,0,nq_adj,1)+1
                        NXQ(-1,NXQ(-1,0,nq_adj,1),nq_adj,1)=nq_ele
                      ENDIF
                    ENDIF
                  ENDIF
                ENDDO !nii
                IF(GENER(-1).AND.(.NOT.BACK_BIFUR))
     &            NQSTART(1)=NQSTART(1)+1
              ENDIF
! Element in +xi 1
              IF(NXI(1,0,ne).GT.0) THEN
                IF(TYPE(1:8).EQ.'CORONARY'
     &            .AND.(NXI(1,0,ne).GE.2)) THEN
                  FOR_BIFUR=.TRUE.
                ELSE
                  FOR_BIFUR=.FALSE.
                ENDIF
                DO nii=1,NXI(1,0,ne)
                  IF(NXI(1,nii,ne).LT.ne) THEN
                    GENER(1)=.TRUE.
! If we are looking forward onto the start of the next element
                    IF(NPNE(NNT(nb),nb,ne).EQ.NPNE(1,nb,NXI(1,nii,ne)))
     &                THEN
                      IF(.NOT.FOR_BIFUR) THEN
                        NQNE(ne,NQXI(1,SCHEME))=NQNE(NXI(1,nii,ne),1)
                      ELSE
                        nq_ele=nqnew+(NQFINISH(1)-NQSTART(1))
                        nq_adj=NQNE(NXI(1,nii,ne),1)
                        WRITE(OP_STRING(1),
     &                    '(''>>Increase NQM >= '',I12)') MAX(nq_ele,
     &                    nq_adj)
                        CALL STRING_TRIM(OP_STRING(1),IBEG,IEND)
                        CALL ASSERT(nq_ele.LE.NQM.AND.nq_adj.LE.NQM,
     &                    OP_STRING(1)(IBEG:IEND),ERROR,*9999)
                        NXQ(1,0,nq_ele,1)=NXQ(1,0,nq_ele,1)+1
                        NXQ(1,NXQ(1,0,nq_ele,1),nq_ele,1)=nq_adj
                        NXQ(-1,0,nq_adj,1)=NXQ(-1,0,nq_adj,1)+1
                        NXQ(-1,NXQ(-1,0,nq_adj,1),nq_adj,1)=nq_ele
                      ENDIF
! IF we are looking forward onto the end of the next element
                    ELSE IF(NPNE(NNT(nb),nb,ne).EQ.
     &                  NPNE(NNT(nb),nb,NXI(1,nii,ne))) THEN
                      IF(.NOT.FOR_BIFUR) THEN
                        NQNE(ne,NQXI(1,SCHEME))=NQNE(NXI(1,nii,ne),
     &                    NQET(NQS(NXI(1,nii,ne))))
                      ELSE
                        nq_ele=nqnew+(NQFINISH(1)-NQSTART(1))
                        nq_adj=NQNE(NXI(1,nii,ne),
     &                    NQET(NQS(NXI(1,nii,ne))))
                        WRITE(OP_STRING(1),
     &                    '(''>>Increase NQM >= '',I12)') MAX(nq_ele,
     &                    nq_adj)
                        CALL STRING_TRIM(OP_STRING(1),IBEG,IEND)
                        CALL ASSERT(nq_ele.LE.NQM.AND.nq_adj.LE.NQM,
     &                    OP_STRING(1)(IBEG:IEND),ERROR,*9999)
                        NXQ(1,0,nq_ele,1)=NXQ(1,0,nq_ele,1)+1
                        NXQ(1,NXQ(1,0,nq_ele,1),nq_ele,1)=nq_adj
                        NXQ(1,0,nq_adj,1)=NXQ(1,0,nq_adj,1)+1
                        NXQ(1,NXQ(1,0,nq_adj,1),nq_adj,1)=nq_ele
                      ENDIF
                    ENDIF
                  ENDIF
                ENDDO !nii
                IF(GENER(1).AND.(.NOT.FOR_BIFUR))
     &            NQFINISH(1)=NQFINISH(1)-1
              ENDIF
             ENDIF 
! Create new grid points where required
              DO nqe=NQSTART(1),NQFINISH(1)
                NQNE(ne,nqe)=nqnew
                IF(nqnew.LE.NQM) THEN
                  NXQ(0,1,nqnew,1)=nqnew
                  NENQ(0,nqnew)=0
                  NWQ(4,nqnew,1)=0
                ENDIF ! <NQM
                nqnew=nqnew+1
              ENDDO !nq
              
! Calculate grid point connectivity
              DO nqe=1,NQET(SCHEME)
                nq=NQNE(ne,nqe)
                mq1=0
                mq2=0
                IF(nqe.GT.1) mq1=NQNE(ne,nqe-1)
                IF(nqe.LT.NQET(SCHEME)) mq2=NQNE(ne,nqe+1)
                
! Connections in the xi 1 direction
                IF(mq1.EQ.0) NWQ(4,nq,1)=1
                IF(mq2.GT.0) THEN
                  NXQ(1,0,nq,1)=NXQ(1,0,nq,1)+1
                  WRITE(OP_STRING(1),
     &              '(''>>Increase NQM >= '',I12)') MAX(mq2,nq)
                  CALL STRING_TRIM(OP_STRING(1),IBEG,IEND)
                  CALL ASSERT(nq.LE.NQM,OP_STRING(1)(IBEG:IEND),ERROR,
     &              *9999)      
                  NXQ(1,NXQ(1,0,nq,1),nq,1)=mq2
                  CALL ASSERT(mq2.LE.NQM,OP_STRING(1)(IBEG:IEND),ERROR,
     &              *9999)
                  NXQ(-1,0,mq2,1)=NXQ(-1,0,mq2,1)+1
                  NXQ(-1,NXQ(-1,0,mq2,1),mq2,1)=nq
                ELSE
                  NWQ(4,nq,1)=-1
                ENDIF
              ENDDO !nqe
              
! If there are two xi directions
            ELSE IF(NQXI(0,SCHEME).EQ.2) THEN
! Element in -xi 1
              IF((NXI(-1,1,ne).LT.ne).AND.(NXI(-1,1,ne).GT.0)) THEN
                GENER(-1)=.TRUE.
                DO nqe2=1,NQXI(2,SCHEME)
                  NQNE(ne,((nqe2-1)*NQXI(1,SCHEME))+1)=NQNE(NXI(-1,
     &              1,ne),(nqe2*NQXI(1,NQS(NXI(-1,1,ne)))))
                ENDDO !nqe2
              ENDIF
! Element in +xi 1
              IF((NXI(1,1,ne).LT.ne).AND.(NXI(1,1,ne).GT.0)) THEN
                GENER(1)=.TRUE.
                DO nqe2=1,NQXI(2,SCHEME)
                  NQNE(ne,(nqe2*NQXI(1,SCHEME)))=NQNE(NXI(1,1,ne),
     &              ((nqe2-1)*NQXI(1,NQS(NXI(1,1,ne))))+1)
                ENDDO !nqe2
              ENDIF
! Element in -xi 2
              IF((NXI(-2,1,ne).LT.ne).AND.(NXI(-2,1,ne).GT.0)) THEN
                GENER(-2)=.TRUE.
                DO nqe1=1,NQXI(1,SCHEME)
                  NQNE(ne,nqe1)=NQNE(NXI(-2,1,ne),(nqe1+(NQXI(1,
     &              NQS(NXI(-2,1,ne)))*(NQXI(2,NQS(NXI(-2,1,ne)))-1))))
                ENDDO !nqe1
              ENDIF
! Element in +xi 2
              IF((NXI(2,1,ne).LT.ne).AND.(NXI(2,1,ne).GT.0)) THEN
                GENER(2)=.TRUE.
                DO nqe1=1,NQXI(1,SCHEME)
                  NQNE(ne,(nqe1+(NQXI(1,SCHEME)*(NQXI(2,SCHEME)-1))))=
     &              NQNE(NXI(2,1,ne),nqe1)
                ENDDO !nqe1
              ENDIF

! Update new grid point start/finish counters
              DO ni=1,NQXI(0,SCHEME)
                IF(GENER(-ni)) NQSTART(ni)=NQSTART(ni)+1
                IF(GENER(ni)) NQFINISH(ni)=NQFINISH(ni)-1
              ENDDO !ni
              
! Create new grid points where required
              DO nqe2=NQSTART(2),NQFINISH(2)
                DO nqe1=NQSTART(1),NQFINISH(1)
                  nqe=nqe1+((nqe2-1)*NQXI(1,SCHEME))
                  NQNE(ne,nqe)=nqnew
                  IF(nqnew.LE.NQM) THEN
                    NXQ(0,1,nqnew,1)=nqnew
                    NENQ(0,nqnew)=0
                    DO nstep=4,5
                      NWQ(nstep,nqnew,1)=0
                    ENDDO !ni
                  ENDIF ! <NQM
                  nqnew=nqnew+1
                ENDDO !nqe2
              ENDDO !nqe
              
! Calculate grid point connectivity
              DO nqe2=1,NQXI(2,SCHEME)
                DO nqe1=1,NQXI(1,SCHEME)
                  nqe=nqe1+(nqe2-1)*NQXI(1,SCHEME)
                  nq=NQNE(ne,nqe)
                  mq1=0
                  mq2=0
                  mq3=0
                  mq4=0
                  IF(nqe1.GT.1) mq1=NQNE(ne,nqe-1)
                  IF(nqe1.LT.NQXI(1,SCHEME)) mq2=NQNE(ne,nqe+1)
                  IF(nqe2.GT.1) mq3=NQNE(ne,nqe-NQXI(1,SCHEME))
                  IF(nqe2.LT.NQXI(2,SCHEME)) mq4=NQNE(ne,nqe+NQXI(1,
     &              SCHEME))
                  
                  IF(nq.LE.NQM) THEN
! Connections in the xi 1 direction
                    IF(mq1.EQ.0) NWQ(4,nq,1)=NWQ(4,nq,1)+1
                    IF(mq2.GT.0.AND.mq2.LE.NQM) THEN
C                   NXQ(1,0,nq,1)=1
C                   NXQ(1,1,nq,1)=mq2
C                   NXQ(-1,0,mq2,1)=1
C                   NXQ(-1,1,mq2,1)=nq
                      IF(mq2.NE.NXQ(1,1,nq,1)) THEN
                        NXQ(1,0,nq,1)=NXQ(1,0,nq,1)+1
                        NXQ(1,NXQ(1,0,nq,1),nq,1)=mq2
                        NXQ(-1,0,mq2,1)=NXQ(-1,0,mq2,1)+1
                        NXQ(-1,NXQ(-1,0,mq2,1),mq2,1)=nq
                      ENDIF
                    ELSE
                      NWQ(4,nq,1)=NWQ(4,nq,1)-1
                    ENDIF

! Connections in the  xi 2 direction
                    IF(mq3.EQ.0) NWQ(5,nq,1)=NWQ(5,nq,1)+1
                    IF(mq4.GT.0) THEN
C                   NXQ(2,0,nq,1)=1
C                   NXQ(2,1,nq,1)=mq4
C                   NXQ(-2,0,mq4,1)=1
C                   NXQ(-2,1,mq4,1)=nq
                      IF(mq4.NE.NXQ(2,1,nq,1).AND.mq4.LE.NQM) THEN
                        NXQ(2,0,nq,1)=NXQ(2,0,nq,1)+1
                        NXQ(2,NXQ(2,0,nq,1),nq,1)=mq4
                        NXQ(-2,0,mq4,1)=NXQ(-2,0,mq4,1)+1
                        NXQ(-2,NXQ(-2,0,mq4,1),mq4,1)=nq
                      ENDIF
                    ELSE
                      NWQ(5,nq,1)=NWQ(5,nq,1)-1
                    ENDIF
                  ENDIF ! nq <= NQM
                ENDDO !nqe1
              ENDDO !nqe2

! If there are three xi directions
            ELSE IF(NQXI(0,SCHEME).EQ.3) THEN
              IF(KTYP3B.NE.2) THEN
! Element in -xi 1
                IF((NXI(-1,1,ne).LT.ne).AND.(NXI(-1,1,ne).GT.0)) THEN
                  GENER(-1)=.TRUE.
                  DO nqe3=1,NQXI(3,SCHEME)
                    DO nqe2=1,NQXI(2,SCHEME)
                      nqe=1+((nqe2-1)*NQXI(1,SCHEME))+((nqe3-1)*
     &                  NQXI(1,SCHEME)*NQXI(2,SCHEME))
                      NQNE(ne,nqe)=NQNE(NXI(-1,1,ne),(nqe2*NQXI(1,
     &                  NQS(NXI(-1,1,ne))))+((nqe3-1)*NQXI(1,NQS(NXI(-1,
     &                  1,ne)))*NQXI(2,NQS(NXI(-1,1,ne)))))
                    ENDDO !nqe2
                  ENDDO !nqe3
                ENDIF
! Element in +xi 1
                IF((NXI(1,1,ne).LT.ne).AND.(NXI(1,1,ne).GT.0)) THEN
                  GENER(1)=.TRUE.
                  DO nqe3=1,NQXI(3,SCHEME)
                    DO nqe2=1,NQXI(2,SCHEME)
                      nqe=(nqe2*NQXI(1,SCHEME))+((nqe3-1)*
     &                  NQXI(1,SCHEME)*NQXI(2,SCHEME))
                      NQNE(ne,nqe)=NQNE(NXI(1,1,ne),1+((nqe2-1)*
     &                  NQXI(1,NQS(NXI(1,1,ne))))+((nqe3-1)*
     &                  NQXI(1,NQS(NXI(1,1,ne)))*NQXI(2,NQS(NXI(1,1,
     &                  ne)))))
                    ENDDO !nqe2
                  ENDDO !nqe3
                ENDIF
! Element in -xi 2
                IF((NXI(-2,1,ne).LT.ne).AND.(NXI(-2,1,ne).GT.0)) THEN
                  GENER(-2)=.TRUE.
                  DO nqe3=1,NQXI(3,SCHEME)
                    DO nqe1=1,NQXI(1,SCHEME)
                      nqe=nqe1+((nqe3-1)*NQXI(1,SCHEME)*NQXI(2,SCHEME))
                      NQNE(ne,nqe)=NQNE(NXI(-2,1,ne),nqe1+
     &                  ((NQXI(2,NQS(NXI(-2,1,ne)))-1)*NQXI(1,
     &                  NQS(NXI(-2,1,ne))))+
     &                  ((nqe3-1)*NQXI(1,NQS(NXI(-2,1,ne)))*NQXI(2,
     &                  NQS(NXI(-2,1,ne)))))
                    ENDDO !nqe1
                  ENDDO !nqe3
                ENDIF
! Element in +xi 2
                IF((NXI(2,1,ne).LT.ne).AND.(NXI(2,1,ne).GT.0)) THEN
                  GENER(2)=.TRUE.
                  DO nqe3=1,NQXI(3,SCHEME)
                    DO nqe1=1,NQXI(1,SCHEME)
                      nqe=nqe1+((NQXI(2,SCHEME)-1)*NQXI(1,SCHEME))+
     &                  ((nqe3-1)*NQXI(1,SCHEME)*NQXI(2,SCHEME))
                      NQNE(ne,nqe)=NQNE(NXI(2,1,ne),nqe1+((nqe3-1)*
     &                  NQXI(1,NQS(NXI(2,1,ne)))*NQXI(2,NQS(NXI(2,1,
     &                  ne)))))
                    ENDDO !nqe1
                  ENDDO !nqe3
                ENDIF
! Element in -xi 3
                IF((NXI(-3,1,ne).LT.ne).AND.(NXI(-3,1,ne).GT.0)) THEN
                  GENER(-3)=.TRUE.
                  DO nqe2=1,NQXI(2,SCHEME)
                    DO nqe1=1,NQXI(1,SCHEME)
                      nqe=nqe1+((nqe2-1)*NQXI(1,SCHEME))
                      NQNE(ne,nqe)=NQNE(NXI(-3,1,ne),nqe1+((nqe2-1)*
     &                  NQXI(1,NQS(NXI(-3,1,ne))))+((NQXI(3,
     &                  NQS(NXI(-3,1,ne)))-1)*NQXI(1,NQS(NXI(-3,1,ne)))*
     &                  NQXI(2,NQS(NXI(-3,1,ne)))))
                    ENDDO !nqe1
                  ENDDO !nqe2
                ENDIF
! Element in +xi 3
                IF((NXI(3,1,ne).LT.ne).AND.(NXI(3,1,ne).GT.0)) THEN
                  GENER(3)=.TRUE.
                  DO nqe2=1,NQXI(2,SCHEME)
                    DO nqe1=1,NQXI(1,SCHEME)
                      nqe=nqe1+((nqe2-1)*NQXI(1,SCHEME))+((NQXI(3,
     &                  SCHEME)-1)*NQXI(1,SCHEME)*NQXI(2,SCHEME))
                      NQNE(ne,nqe)=NQNE(NXI(3,1,ne),nqe1+((nqe2-1)*
     &                  NQXI(1,NQS(NXI(3,1,ne)))))
                    ENDDO !nqe1
                  ENDDO !nqe2
                ENDIF
              ENDIF
              
! Update new grid point start/finish counters
              DO ni=1,NQXI(0,SCHEME)
                IF(GENER(-ni)) NQSTART(ni)=NQSTART(ni)+1
                IF(GENER(ni)) NQFINISH(ni)=NQFINISH(ni)-1
              ENDDO !ni
              
! Create new grid points where required
              DO nqe3=NQSTART(3),NQFINISH(3)
                DO nqe2=NQSTART(2),NQFINISH(2)
                  DO nqe1=NQSTART(1),NQFINISH(1)
                    nqe=nqe1+((nqe2-1)*NQXI(1,SCHEME))+((nqe3-1)*
     &                NQXI(1,SCHEME)*NQXI(2,SCHEME))
                    NQNE(ne,nqe)=nqnew
                    IF(nqnew.LE.NQM) THEN
                      NXQ(0,1,nqnew,1)=nqnew
                      NENQ(0,nqnew)=0
                      DO nstep=4,6
                        NWQ(nstep,nqnew,1)=0
                      ENDDO !ni
                    ENDIF ! <NQM
                    nqnew=nqnew+1
                  ENDDO !nqe1
                ENDDO !nqe2
              ENDDO !nqe3
              
! Calculate grid point connectivity
              DO nqe3=1,NQXI(3,SCHEME)
                DO nqe2=1,NQXI(2,SCHEME)
                  DO nqe1=1,NQXI(1,SCHEME)
                    nqe=nqe1+((nqe2-1)*NQXI(1,SCHEME))+((nqe3-1)*
     &                NQXI(1,SCHEME)*NQXI(2,SCHEME))
                    nq=NQNE(ne,nqe)
                    mq1=0
                    mq2=0
                    mq3=0
                    mq4=0
                    mq5=0
                    mq6=0
                    IF(nqe1.GT.1) mq1=NQNE(ne,nqe-1)
                    IF(nqe1.LT.NQXI(1,SCHEME)) mq2=NQNE(ne,nqe+1)
                    IF(nqe2.GT.1) mq3=NQNE(ne,nqe-NQXI(1,SCHEME))
                    IF(nqe2.LT.NQXI(2,SCHEME)) mq4=NQNE(ne,nqe+NQXI(1,
     &                SCHEME))
                    IF(nqe3.GT.1) mq5=NQNE(ne,nqe-(NQXI(1,SCHEME)*
     &                NQXI(2,SCHEME)))
                    IF(nqe3.LT.NQXI(3,SCHEME)) mq6=NQNE(ne,nqe+(NQXI(1,
     &                SCHEME)*NQXI(2,SCHEME)))
                    
                    IF(nq.LE.NQM) THEN
! Connections in the xi 1 direction
                      IF(mq1.EQ.0) NWQ(4,nq,1)=NWQ(4,nq,1)+1
                      IF((mq2.GT.0).AND.(mq2.LE.NQM)) THEN
C                       IF(mq2.GT.0) THEN
C                       NXQ(1,0,nq,1)=1
C                       NXQ(1,1,nq,1)=mq2
C                       IF(mq2.LE.NQM) THEN
C                       NXQ(-1,0,mq2,1)=1
C                       NXQ(-1,1,mq2,1)=nq
C                       ENDIF
                        IF((mq2.NE.NXQ(1,1,nq,1)).AND.
     &                    (mq2.NE.NXQ(1,2,nq,1))) THEN
                          NXQ(1,0,nq,1)=NXQ(1,0,nq,1)+1
                          NXQ(1,NXQ(1,0,nq,1),nq,1)=mq2
                          NXQ(-1,0,mq2,1)=NXQ(-1,0,mq2,1)+1
                          NXQ(-1,NXQ(-1,0,mq2,1),mq2,1)=nq
                        ENDIF
                      ELSE
                        NWQ(4,nq,1)=NWQ(4,nq,1)-1
                      ENDIF
                      
! Connections in the xi 2 direction
                      IF(mq3.EQ.0) NWQ(5,nq,1)=NWQ(5,nq,1)+1
                      IF((mq4.GT.0).AND.(mq4.LE.NQM)) THEN
C                       IF(mq4.GT.0) THEN
C                         NXQ(2,0,nq,1)=1
C                         NXQ(2,1,nq,1)=mq4
C                       IF(mq4.LE.NQM) THEN
C                         NXQ(-2,0,mq4,1)=1
C                         NXQ(-2,1,mq4,1)=nq
C                       ENDIF
                        IF((mq4.NE.NXQ(2,1,nq,1)).AND.
     &                    (mq4.NE.NXQ(2,2,nq,1))) THEN
                          NXQ(2,0,nq,1)=NXQ(2,0,nq,1)+1
                          NXQ(2,NXQ(2,0,nq,1),nq,1)=mq4
                          NXQ(-2,0,mq4,1)=NXQ(-2,0,mq4,1)+1
                          NXQ(-2,NXQ(-2,0,mq4,1),mq4,1)=nq
                        ENDIF
                      ELSE
                        NWQ(5,nq,1)=NWQ(5,nq,1)-1
                      ENDIF

! Connections in the xi 3 direction
                      IF(mq5.EQ.0) NWQ(6,nq,1)=NWQ(6,nq,1)+1
                      IF((mq6.GT.0).AND.(mq6.LE.NQM)) THEN
C                    IF(mq6.GT.0) THEN
C                      NXQ(3,0,nq,1)=1
C                      NXQ(3,1,nq,1)=mq6
C                      IF(mq6.LE.NQM) THEN
C                        NXQ(-3,0,mq6,1)=1
C                        NXQ(-3,1,mq6,1)=nq
C                      ENDIF
                        IF((mq6.NE.NXQ(3,1,nq,1)).AND.
     &                    (mq6.NE.NXQ(3,2,nq,1))) THEN
                          NXQ(3,0,nq,1)=NXQ(3,0,nq,1)+1
                          NXQ(3,NXQ(3,0,nq,1),nq,1)=mq6
                          NXQ(-3,0,mq6,1)=NXQ(-3,0,mq6,1)+1
                          NXQ(-3,NXQ(-3,0,mq6,1),mq6,1)=nq
                        ENDIF
                      ELSE
                        NWQ(6,nq,1)=NWQ(6,nq,1)-1
                      ENDIF
                    ENDIF ! nq <= NQM
                  ENDDO !nqe1
                ENDDO !nqe2
              ENDDO !nqe3
            ELSE
              IF(KTYP3B.NE.2) THEN
                CALL ASSERT(.FALSE.,'>>Invalid number of xi directions',
     &            ERROR,*9999)
              ENDIF
            ENDIF

! Create the NENQ array
            DO nqq=1,NQET(SCHEME)
              nq=NQNE(ne,nqq)
              IF(nq.LE.NQM) THEN
                NENQ(0,nq)=NENQ(0,nq)+1
                NENQ(NENQ(0,nq),nq)=ne
              ENDIF
            ENDDO !nqq
            
            IF(DOP) THEN
              WRITE(OP_STRING,'(''Element '',I4,'' Region '',I4)') ne,nr
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              DO nq=1,NQET(NQS(ne))
                WRITE(OP_STRING,'(''Grid number '',I6,I6)') nq,NQNE(ne,
     &            nq)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDDO
            ENDIF
          ENDDO !element

          NQR(2,nr)=nqnew-1

! Check for and flag external grid points
          NITB=NQXI(0,NQS(NEELEM(1,nr)))
          IK=MAX(0,NITB-2) !zero for 1,2-D, one for 3-D
          IJ=MIN(NITB-1,1) !zero for 1-D, one for 2,3-D
          DO nq=NQR(1,nr),NQR(2,nr)
            IF(nq.LE.NQM) THEN
              INTERNAL=.TRUE.
              ISUM=0
              DO nstep=4,3+NITB
                IF(NWQ(nstep,nq,1).NE.0)
     &            NWQ(nstep,nq,1)=SIGN(1,NWQ(nstep,nq,1))
              ENDDO !ni
              DO nik=-IK,IK
                DO nij=-IJ,IJ
                  DO nii=-1,1
                    IF(NXQ(nik*3,1,nq,1).LE.NQM) THEN
                      IF(NXQ(nij*2,1,NXQ(nik*3,1,nq,1),1).LE.NQM) THEN
                        mq=NXQ(nii,1,NXQ(nij*2,1,NXQ(nik*3,1,nq,1),1),1)
                        IF(mq.EQ.0) THEN
                          INTERNAL=.FALSE.
                          ISUM=ISUM+1
                        ENDIF
                      ENDIF
                    ENDIF
                  ENDDO !nii
                ENDDO !nij
              ENDDO !nik
              IF(INTERNAL) THEN !internal g.p.
                DO ni=1,2
                  NWQ(ni,nq,1)=0
                ENDDO
              ELSE !on external bdy
                IF(NITB.EQ.2) THEN
                  IF(ISUM.EQ.1) THEN
! Special case for "270 degree" corner (internal corner)
! Find which g.p. is missing, and head away from that.
                    DO nij=-1,1
                      DO nii=-1,1
                        mq=NXQ(nii,1,NXQ(nij*2,1,nq,1),1)
                        IF(mq.EQ.0) THEN
                          NWQ(4,nq,1)=-nii
                          NWQ(5,nq,1)=-nij
                        ENDIF
                      ENDDO
                    ENDDO
                  ENDIF
                ELSE IF(NITB.EQ.3) THEN
                ENDIF
! Now find g.p.s to use for no-flux b.c.s
                mq=nq
                DO nstep=1,2
! Step twice in dirn indicated from NWQ(ni+3,nq,1)
                  DO ni=1,NITB
                    IF(mq.LE.NQM) mq=NXQ(ni*NWQ(ni+3,nq,1),1,mq,1)
                  ENDDO
                  NWQ(nstep,nq,1)=mq
                ENDDO
              ENDIF !int/ext
            ENDIF
          ENDDO !nq
C        DO nq=NQR(1,nr),NQR(2,nr)
C          IF(NITB.EQ.2) THEN !2 xi directions
C            !Check for cusp which should be external
C            DO ni=1,NITB
C              IF(NXQ(ni,0,nq,1).GT.1) THEN
C                NWQ(1,nq,1)=NXQ(-ni,1,nq,1)
C                NWQ(2,nq,1)=NXQ(-ni,1,NWQ(1,nq,1),1)
CC                DO nii=1,NXQ(ni,0,nq,1)
CC                  NXQ(ni,nii,nq,1)=0
CC                ENDDO
CC                NXQ(ni,0,nq,1)=0
C              ELSEIF(NXQ(-ni,0,nq,1).GT.1) THEN
C                NWQ(1,nq,1)=NXQ(ni,1,nq,1)
C                NWQ(2,nq,1)=NXQ(ni,1,NWQ(1,nq,1),1)
CC                DO nii=1,NXQ(-ni,0,nq,1)
CC                  NXQ(-ni,nii,nq,1)=0
CC                ENDDO
CC                NXQ(-ni,0,nq,1)=0
C              ENDIF !nxq
C            ENDDO !ni
C          ENDIF !2d
C        ENDDO  !nq
        ENDDO !region

        NQT=nqnew-1 !Update total number of grid points
        WRITE(OP_STRING(1),
     &    '(''>>Increase NQM >= '',I12)') NQT
        CALL STRING_TRIM(OP_STRING(1),IBEG,IEND)
        CALL ASSERT(nqnew-1.LE.NQM,OP_STRING(1)(IBEG:IEND),ERROR,*9999)

        WRITE(OP_STRING(1),'(''Minimum NQM needed: '',I12)') NQT
        CALL STRING_TRIM(OP_STRING(1),IBEG,IEND)
        CALL WRITES(IODI,OP_STRING(1),ERROR,*9999)
        
        IF(TYPE(1:8).EQ.'CORONARY') THEN
          WRITE(OP_STRING,
     &      '(''smallest length'',F12.8)') SMALL_LENGTH
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

        IF(DOP) THEN
          IF(TYPE(1:8).EQ.'CORONARY') THEN
            WRITE(OP_STRING,
     &        '(''smallest length'',F12.8)') SMALL_LENGTH
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ELSE
            DO nq=1,NQT
              WRITE(OP_STRING,'(''Grid point '',I6)') nq
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'(''Element numbers'',8I4)')
     &          (NENQ(ne,nq),ne=1,NENQ(0,nq))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO !nq
            DO nq=1,NQT
              WRITE(OP_STRING,
     &          '(''nq,NWQ(1..2,nq,1)'',3I6)') nq,NWQ(1,nq,1),
     &          NWQ(2,nq,1)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO !nq
            DO ni=-NIM,NIM
              DO i=0,4
                DO nq=1,NQT
                  WRITE(OP_STRING,'(''ni,i,nq,nxq'',4I5)') ni,i,nq,
     &              NXQ(ni,i,nq,1)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDDO !ni
              ENDDO !i
            ENDDO !nq
          ENDIF
        ENDIF

      ELSE IF(USE_LAT.EQ.1) THEN
C DAH 17Jun02: Adding in lattice grid points. If lattice grid points
C       are used, then the following code is skipped.

C       Allocate memory if needed for lattice grid points
C GBS 16-Jul-2003 Changed to external routine for allocation
C       and dimensions only as large as needed:
C         NEQM is set to NEELEM(0,0) in DEGRID
        CALL ALLOCATE_LATTICE(NLATNE_PTR,NLATNQ_PTR,NLATPNQ_PTR,
     &    NQNLAT_PTR,ERROR,*9999)
      
        NLATXIM=0
        DO nqsc=1,NQSCT
          DO ni=1,NQXI(0,nqsc)
            IF(NQXI(ni,nqsc).GT.NLATXIM) THEN
              NLATXIM=NQXI(ni,nqsc)
            ENDIF
          ENDDO
        ENDDO
        CALL CREATE_LATTICE(NBJ,NEELEM,NENP,NENQ,%VAL(NLATNE_PTR),
     &    %VAL(NLATNQ_PTR),%VAL(NLATPNQ_PTR),NNB,NPNE,NQNE,
     &    %VAL(NQNLAT_PTR),NQS,NQSCNB,NQXI,NRLIST,NWQ,XIQ,ERROR,*9999)
      ENDIF

C KAT 2005-02-18:  Was there a reason for this?  What is it used for?
C       !Reset temp work memory to zero
C       DO nq=1,NQT
C         NWQ(4,nq,1)=0
C         NWQ(5,nq,1)=0
C         NWQ(6,nq,1)=0
C       ENDDO

      CALL EXITS('IPGRID')
      RETURN
 9999 CALL ERRORS('IPGRID',ERROR)
      CALL EXITS('IPGRID')
      RETURN 1
      END
