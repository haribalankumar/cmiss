      SUBROUTINE BASIS6(IBT,IDO,INP,nb,NDET,NGAP,DET,PG,WG,XIG,ERROR,*)

C#### Subroutine: BASIS6
C###  Description:
C###    BASIS6 inputs parameters for BE simplex/serendipity/sector
C###    basis functions.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  nb,NDET(NBFM,0:NNM),NGAP(NIM,NBM)
      REAL*8 DET(NBFM,0:NNM,NGM,6),PG(NSM,NUM,NGM,NBM),
     '  WG(NGM,NBM),XIG(NIM,NGM,NBM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER COLLAPSED(3),i,ICHAR,INFO,k,NBTOP,ni,ni1,ni2,NI2_TOP,
     '  nk,nkinc,nn,nnmax,NOQUES,numnn,NUM_COLLAPSED,NUMHERMXI,
     '  NUMNODES(3),POSITION(4)
      CHARACTER CHAR1*100,CHAR2*3,CHAR3*3
      LOGICAL ATCOLLAPSE,FILEIP,ISATCOLLAPSE

      CALL ENTERS('BASIS6',*9999)

      FILEIP=.FALSE.
      NOQUES=0
      INFO=2
      ICHAR=999

      DO K=1,100
        CHAR1(K:K)=' '
      ENDDO

C cpb 31/10/96 Can only have 2 xi direcitons for BEM sector type bases
C      FORMAT='($,'' Enter the number of Xi-coordinates [1]: '',I1)'
C      IF(IOTYPE.EQ.3) IDATA(1)=NIT(nb)
C      CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NJT,
C     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C      IF(iotype.ne.3) NIT(nb)=IDATA(1)

      NIT(nb)=2

      IF(iotype.ne.3) THEN
        NUT(nb)=NIT(nb)*NIT(nb)+2
        CALL ASSERT(NUT(nb).LE.NUM,'>>Increase NUM',ERROR,*9999)
        NGT(nb)=1
        NNT(nb)=1
        DO ni=1,NIT(nb)
          NGAP(ni,nb)=1
          IBT(1,ni,nb)=1
C news MPN 9-Jul-96 Initialise IBT
          IBT(2,ni,nb)=0
          IBT(3,ni,nb)=0
        ENDDO !ni
      ENDIF

      FORMAT='('' Specify the type of basis [1]:'''//
     '  '/''   (1) Simplex'''//
     '  '/''   (2) Serendipity'''//
     '  '/''   (3) Sector'''//
     '  '/$,''    '',I1)'
      IF(IOTYPE.EQ.3) IDATA(1)=NBSC(1,nb)
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,3,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(iotype.ne.3) NBSC(1,nb)=IDATA(1)

C cpb 31/10/96 Adding general sectors
      IF(NBSC(1,nb).LE.2) THEN
        FORMAT='('' Specify the basis degree [1]:'''//
     '    '/''   (1) Linear'''//
     '    '/''   (2) Quadratic'''//
     '    '/''   (3) Cubic '''//
     '    '/''   (4) Hermite'''//
     '    '/$,''    '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=NBSC(2,nb)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,4,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(iotype.ne.3) NBSC(2,nb)=IDATA(1)
      ELSE
C       Find the basis type in each direction
        NUM_COLLAPSED=0
        DO ni=1,NIT(nb)
          WRITE(CHAR1,'(I1)') ni
          FORMAT='(/'' Specify the basis type in the Xi('//CHAR1(1:1)
     '      //') direction [1]:'''//
     '      '/''   (1) Linear Lagrange'''//
     '      '/''   (2) Quadratic Lagrange'''//
     '      '/''   (3) Cubic Lagrange'''//
     '      '/''   (4) Hermite'''//
     '      '/$,''    '',I1)'
          IF(IOTYPE.EQ.3) THEN
            IF(IBT(1,ni,nb).EQ.2) THEN
              IDATA(1)=4
            ELSE
              IDATA(1)=IBT(2,ni,nb)
            ENDIF
          ENDIF
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,4,
     &      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) IBT(2,ni,nb)=IDATA(1)
          IF(nb.GT.1) THEN
            IF(NIT(nb).EQ.NIT(nb-1).AND.
     '        (NBC(nb-1).NE.5.OR.NBC(nb-1).EQ.6)) THEN
              IDEFLT(1)=NGAP(ni,nb-1)
            ENDIF
          ELSE
            IF(IBT(2,ni,nb).EQ.1) THEN
              IDEFLT(1)=2
            ELSE IF(IBT(2,ni,nb).GT.1) THEN
              IDEFLT(1)=3
            ENDIF
          ENDIF
          IDEFLT(1)=3
          IDEFLT(2)=6
          WRITE(CHAR2,'(I3)') IDEFLT(1)
          WRITE(CHAR3,'(I3)') IDEFLT(2)
          FORMAT='('' Enter the number of Gauss points in the Xi('
     '      //CHAR1(1:1)//') direction for the'''
     '      //'/$,'' low and high order schemes '
     '      //'['//CHAR2(1:3)//','//CHAR3(1:3)//']: '',I3,I3)'
          IF(IOTYPE.EQ.3) THEN
            IDATA(1)=NGAP(ni,NBASEF(nb,NBASEF(nb,0)-1))
            IDATA(2)=NGAP(ni,NBASEF(nb,1))
          ENDIF
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      2,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,IMAX,
     &      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(iotype.ne.3) THEN
            NGLIMITS(ni,nb,1)=IDATA(1) !Low order
            NGLIMITS(ni,nb,2)=IDATA(2) !High order
            IF(NGLIMITS(ni,nb,1).GT.NGLIMITS(ni,nb,2)) THEN
C             Entered in the wrong order
              NGLIMITS(ni,nb,1)=IDATA(2)
              NGLIMITS(ni,nb,2)=IDATA(1)
            ENDIF
            NGAP(ni,nb)=NGLIMITS(ni,nb,2) !High order scheme
          ENDIF
          FORMAT='($,'' Is the Xi('//CHAR1(1:1)//') direction '
     '      //'collapsed [N]? '',A)'
          IF(IOTYPE.EQ.3) THEN
            IF((IBT(1,ni,nb).EQ.5).OR.(IBT(1,ni,nb).EQ.6)) THEN
              ADATA(1)='Y'
              NUM_COLLAPSED=NUM_COLLAPSED+1
              COLLAPSED(NUM_COLLAPSED)=ni
            ELSE
              ADATA(1)='N'
            ENDIF
          ENDIF
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,3,
     &      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            IF(ADATA(1).EQ.'Y') THEN
              NUM_COLLAPSED=NUM_COLLAPSED+1
              COLLAPSED(NUM_COLLAPSED)=ni
            ELSE
              IF(IBT(2,ni,nb).EQ.4) THEN
                IBT(1,ni,nb)=2
                IBT(2,ni,nb)=1
              ELSE
                IBT(1,ni,nb)=1
              ENDIF
            ENDIF
          ENDIF
        ENDDO !ni
        CALL ASSERT(NUM_COLLAPSED.LE.(NIT(nb)-1),
     '    '>>Invalid collapsing of element',ERROR,*9999)
        CALL ASSERT(NUM_COLLAPSED.GT.0,
     '    '>>Use normal tensor product basis function',ERROR,*9999)
        IF(COLLAPSED(1).EQ.1) THEN
          IBT(3,COLLAPSED(1),nb)=2
        ELSE
          IBT(3,COLLAPSED(1),nb)=1
        ENDIF
        WRITE(CHAR1,'(I1)') IBT(3,COLLAPSED(1),nb)
        FORMAT='(/'' Specify the end of collapse [1]:'''//
     '    '/''   (1) Xi=0 end of Xi('//CHAR1(1:1)//')'''//
     '    '/''   (2) Xi=1 end of Xi('//CHAR1(1:1)//')'''//
     '    '/$,''    '',I1)'
        IF(IOTYPE.EQ.3) THEN
          IF(IBT(1,COLLAPSED(1),nb).EQ.5) THEN
            IDATA(1)=1
          ELSE
            IDATA(1)=2
          ENDIF
        ENDIF
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '    1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,2,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          IF(IDATA(1).EQ.1) THEN
            DO i=1,NUM_COLLAPSED
              IBT(1,COLLAPSED(i),nb)=5
            ENDDO !i
            IF(IBT(1,IBT(3,COLLAPSED(1),nb),nb).EQ.2)
     '        IBT(2,IBT(3,COLLAPSED(1),nb),nb)=2
          ELSE
            DO i=1,NUM_COLLAPSED
              IBT(1,COLLAPSED(i),nb)=6
            ENDDO !i
            IF(IBT(1,IBT(3,COLLAPSED(1),nb),nb).EQ.2)
     '        IBT(2,IBT(3,COLLAPSED(1),nb),nb)=3
          ENDIF
        ENDIF
      ENDIF

C cpb 11/9/95 Adding zero cross derivatives.
      NUMHERMXI=0
      DO ni=1,NIT(nb)
        IF(((IBT(1,ni,nb).EQ.5.OR.IBT(1,ni,nb).EQ.6).AND.
     '    IBT(2,ni,nb).EQ.4).OR.IBT(1,ni,nb).EQ.2) THEN
          NUMHERMXI=NUMHERMXI+1
        ENDIF
      ENDDO !ni
      IF(NUMHERMXI.GE.2) THEN
        FORMAT='($,'' Do you want to set cross derivatives to zero '
     '    //'[N]? '',A)'
        IF(IOTYPE.EQ.3) THEN
          IF(NBCD(nb).EQ.1) THEN
            ADATA(1)='Y'
          ELSE
            ADATA(1)='N'
          ENDIF
        ENDIF
        CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '    1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          IF(ADATA(1).EQ.'Y') THEN
            NBCD(nb)=1
          ENDIF
        ENDIF
      ENDIF

      IF(NBSC(1,nb).EQ.1) THEN
        nn=0
        IF(NBSC(2,nb).EQ.4)THEN !AJP 29-5-93
          NI2_TOP=1
        ELSE
          NI2_TOP=NBSC(2,nb)
        ENDIF
        DO ni2=1,NI2_TOP+1
          DO ni1=1,NI2_TOP+2-ni2
            nn=nn+1
            INP(nn,1,nb)=ni1
            INP(nn,2,nb)=ni2
            IF(NIT(nb).EQ.3) INP(nn,3,nb)=1
          ENDDO
        ENDDO
      ELSE IF(NBSC(1,nb).EQ.2) THEN
        CALL ASSERT(NBSC(2,nb).LT.4,'>>Hermite serendipity not set up',
     '    ERROR,*9999)
      ELSE IF(NBSC(1,nb).EQ.3) THEN
C       Determine the maximum number of nodes
        nnmax=1
        DO ni=1,NIT(nb)
          IF((IBT(1,ni,nb).EQ.2).OR.(IBT(1,ni,nb).EQ.5.AND.
     '      IBT(2,ni,nb).EQ.4).OR.(IBT(1,ni,nb).EQ.6.AND.
     '      IBT(2,ni,nb).EQ.4)) THEN
            numnn=2
          ELSE
            numnn=IBT(2,ni,nb)+1
          ENDIF
          NUMNODES(ni)=numnn
          nnmax=nnmax*numnn
          POSITION(ni)=1
        ENDDO
C       Determine the actual number of nodes in the collapsed tensor
C       product and their positions (INP)
        nn=0
        DO i=1,nnmax
          ATCOLLAPSE=.FALSE.
          DO ni=1,NIT(nb)
            IF((IBT(1,ni,nb).EQ.5.AND.
     '        POSITION(IBT(3,COLLAPSED(1),nb)).EQ.1).OR.
     '        (IBT(1,ni,nb).EQ.6.AND.
     '        POSITION(IBT(3,COLLAPSED(1),nb)).EQ.
     '        NUMNODES(IBT(3,COLLAPSED(1),nb)))) ATCOLLAPSE=.TRUE.
          ENDDO !ni
          IF(ATCOLLAPSE) THEN
            IF(POSITION(COLLAPSED(1)).EQ.1) THEN
              nn=nn+1
              IF(nn.LE.NNM) THEN
                DO ni=1,NIT(nb)
                  INP(nn,ni,nb)=POSITION(ni)
                ENDDO !ni
              ENDIF
            ENDIF
          ELSE
            nn=nn+1
            IF(nn.LE.NNM) THEN
              DO ni=1,NIT(nb)
                INP(nn,ni,nb)=POSITION(ni)
              ENDDO !ni
            ENDIF
          ENDIF
          POSITION(1)=POSITION(1)+1
          DO ni=1,NIT(nb)
            IF(POSITION(ni).GT.NUMNODES(ni)) THEN
              POSITION(ni)=1
              POSITION(ni+1)=POSITION(ni+1)+1
            ENDIF
          ENDDO !ni
        ENDDO !i
        NNT(nb)=nn
        CALL ASSERT(NNT(nb).LE.NNM,'>>Increase NNM',ERROR,*9999)
      ENDIF

      IF(NBSC(1,nb).LE.2) THEN !Simplex/Serendipity
        IBT(1,1,nb)=NBSC(1,nb)+2
        IBT(1,2,nb)=NBSC(1,nb)+2
        IBT(2,1,nb)=NBSC(2,nb)
        IBT(2,2,nb)=NBSC(2,nb)
        NNT(nb)=nn
        CALL ASSERT(nn.LE.NNM,'>>Increase NNM',ERROR,*9999)
        CALL ASSERT(NBSC(2,nb).EQ.4,
     '    '>>This type of basis function is not implemented',
     '    ERROR,*9999)
C       Only set up for 2 dimensions initially
        CALL ASSERT(NIT(nb).LE.2,'>>Hermite simplex not implemented'
     '    //' for 3 dimensions',ERROR,*9999)
        IDEFLT(1)=3
        DO ni=1,NIT(nb)
          WRITE(CHAR1,'(I1)') ni
          IDEFLT(1)=3
          IDEFLT(2)=6
          WRITE(CHAR2,'(I3)') IDEFLT(1)
          WRITE(CHAR3,'(I3)') IDEFLT(2)
          FORMAT='('' Enter the number of Gauss points in the Xi('
     '      //CHAR1(1:1)//') direction for the'''
     '      //'/$,'' low and high order schemes '
     '      //'['//CHAR2(1:3)//','//CHAR3(1:3)//']: '',I3,I3)'
          IF(IOTYPE.EQ.3) THEN
            IDATA(1)=NGAP(ni,NBASEF(nb,NBASEF(nb,0)))
            IDATA(2)=NGAP(ni,NBASEF(nb,1))
          ENDIF
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      2,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,7,
     &      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(iotype.ne.3) THEN
            NGLIMITS(ni,nb,1)=IDATA(1) !Low order
            NGLIMITS(ni,nb,2)=IDATA(2) !High order
            IF(NGLIMITS(ni,nb,1).GT.NGLIMITS(ni,nb,2)) THEN
 !           !Entered in the wrong order
              NGLIMITS(ni,nb,1)=IDATA(2)
              NGLIMITS(ni,nb,2)=IDATA(1)
            ENDIF
            NGAP(ni,nb)=NGLIMITS(ni,nb,2) !High order scheme
            IBT(1,ni,nb)=3
            IBT(2,ni,nb)=4 !Special hermite
            NNT(nb)=3
          ENDIF
        ENDDO
        FORMAT='($,'' Is the apex node local node 1 or 3 [1]: '',I1)'
        IF(IOTYPE.EQ.3) THEN
          nn=1
          DO WHILE((NKT(nn,nb).GT.1).AND.(nn.LE.NNT(nb)))
            nn=nn+1
          ENDDO
          CALL ASSERT(nn.LE.NNT(nb),'>>Error in hermite simplex write',
     '      ERROR,*9999)
          IDATA(1)=nn
        ENDIF
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,3,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(iotype.ne.3) THEN
          NST(nb)=0
          NKT(IDATA(1),nb)=1
          DO nn=1,NNT(nb)
            IF(nn.ne.IDATA(1)) NKT(nn,nb)=4
            NST(nb)=NST(nb)+NKT(nn,nb)
          ENDDO
          NKT(0,nb)=4
        ENDIF
        IF(iotype.ne.3)THEN
 !       !Set up standard IDO
c cpb 14/7/95 Temporary nn loop
          DO nn=1,NNT(nb)
            nk=0
            DO ni2=1,2
              DO ni1=1,2
                nk=nk+1
                IDO(nk,nn,1,nb)=ni1
                IDO(nk,nn,2,nb)=ni2
              ENDDO
            ENDDO
            IDO(1,nn,0,nb)=1 !Standard mapping for zeroth component
            IDO(2,nn,0,nb)=2
            IDO(3,nn,0,nb)=4
            IDO(4,nn,0,nb)=7
          ENDDO !nn
        ENDIF
        IF(IDATA(1).ne.3)THEN
 !       !Need to correct INP for apex at node 1.
          INP(1,1,nb)=1
          INP(1,2,nb)=1
          INP(2,1,nb)=1
          INP(2,2,nb)=2
          INP(3,1,nb)=2
          INP(3,2,nb)=2
        ENDIF
        IF(iotype.ne.3)THEN !Don't recalculate PG if writing the file
          NBTOP=NBT !Highest basis function number so far
c PJH 9Sep95     IOD=.FALSE.
          CALL GAUSS6(IBT,IDO,INP,nb,NBTOP,NDET,NGAP,DET,PG,WG,XIG,
     '      ERROR,*9999)
        ENDIF
      ELSE
C cpb 16/7/95 For the moment don't allow any derivatives at the
C collapsed node.
        NST(nb)=0
        NKT(0,nb)=1
        DO nn=1,NNT(nb)
          NKT(nn,nb)=1
          ATCOLLAPSE=ISATCOLLAPSE(IBT(1,1,nb),INP(1,1,nb),nb,nn)
          IF(ATCOLLAPSE) THEN
            DO ni=1,NIT(nb)
              IDO(1,nn,ni,nb)=1
            ENDDO !ni
          ELSE
            DO ni=1,NIT(nb)
              nkinc=1
              IF((IBT(1,ni,nb).EQ.2.AND.(IBT(2,ni,nb).EQ.1.OR.
     '          INP(nn,ni,nb).NE.(IBT(2,ni,nb)-1))).OR.
     '          (IBT(1,ni,nb).EQ.5.AND.IBT(2,ni,nb).EQ.4).OR.
     '          (IBT(1,ni,nb).EQ.6.AND.IBT(2,ni,nb).EQ.4))THEN
                nkinc=2
              ENDIF
              IF(nkinc.EQ.2) THEN
                DO nk=1,NKT(nn,nb)
                  IDO(nk,nn,ni,nb)=1
                  IF(nk+NKT(nn,nb).LE.NKM) THEN
                    DO ni2=1,ni-1
                      IDO(nk+NKT(nn,nb),nn,ni2,nb)=IDO(nk,nn,ni2,nb)
                    ENDDO !ni2
                  ENDIF
                ENDDO !nk
                IF(NKT(nn,nb)*2.LE.NKM) THEN
                  DO nk=NKT(nn,nb)+1,NKT(nn,nb)*2
                    IDO(nk,nn,ni,nb)=2
                  ENDDO !nk
                ENDIF
              ELSE
                IF(NKT(nn,nb).LE.NKM) THEN
                  DO nk=1,NKT(nn,nb)
                    IDO(nk,nn,ni,nb)=1
                  ENDDO !nk
                ENDIF
              ENDIF
              NKT(nn,nb)=NKT(nn,nb)*nkinc
            ENDDO !ni
          ENDIF
C cpb 11/9/95 Adding zero cross derivatives. This is not the best way
C to do this but it is only temporary.
          IF(NBCD(nb).EQ.1) THEN
            IF(NKT(nn,nb).EQ.4) THEN
              NKT(nn,nb)=3
              DO ni=1,NIT(nb)
                IDO(4,nn,ni,nb)=0
              ENDDO !ni
            ENDIF
          ENDIF
          IF(NKT(nn,nb).GT.NKT(0,nb)) NKT(0,nb)=NKT(nn,nb)
          NST(nb)=NST(nb)+NKT(nn,nb)
        ENDDO !nn
        CALL ASSERT(NKT(0,nb).LE.NKM,'>>Increase NKM',ERROR,*9999)
        CALL ASSERT(NST(nb).LE.NSM,'>>Increase NSM',ERROR,*9999)

        IF(IOTYPE.NE.3) THEN
          NBTOP=NBT
          CALL GAUSS6(IBT,IDO,INP,nb,NBTOP,NDET,NGAP,DET,PG,WG,XIG,
     '      ERROR,*9999)
        ENDIF

      ENDIF

      CALL EXITS('BASIS6')
      RETURN
 9999 CALL ERRORS('BASIS6',ERROR)
      CALL EXITS('BASIS6')
      RETURN 1
      END

C 24/2/97 LC archive entire routine : Old BASIS6 - replaced AJP 31-5-93

