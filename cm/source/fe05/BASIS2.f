      SUBROUTINE BASIS2(IBT,IDO,INP,nb,NGAP,PG,WG,XIG,ERROR,*)

C#### Subroutine: BASIS2
C###  Description:
C###    BASIS2 inpust parameters for 2D or 3D Simplex/Serendipity/
C###    Sector basis functions and Gauss-Legendre or triangular
C###    geometry quadrature.

C**** Note:
C**** NBSC(1,nb) is 1,2,3 for simplex,serendipity,sector   in 1-2 plane
C**** NBSC(2,nb)  " 1,2,3,4 " linear,quadr,cubic,hermite   in 1-2 plane
C**** NBSC(3,nb)  " 1,2,3  "  lagrange,simplex,serendipity in 3-dir.n
C**** NBSC(4,nb)  " 1,2,3  "  linear,quadratic,cubic        "     "

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM),IDO(NKM,NNM,0:NIM),INP(NNM,NIM),nb,
     '  NGAP(NIM,NBM)
      REAL*8 PG(NSM,NUM,NGM),WG(NGM),XIG(NIM,NGM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,ICHAR,INCR,INFO,LIMIT,mn,ni,ni1,ni2,NI2_TOP,ni3,
     '  nk,nn,NOQUES,COLLAPSED(3),IDIR(2),NNMAX,NUM_COLLAPSED,
     '  NUMHERMXI,NUMNODES(3),NUMNN,POSITION(4)
      CHARACTER CHAR1*100,CHAR2*100
      LOGICAL ATCOLLAPSE,FILEIP,ISATCOLLAPSE

      CALL ENTERS('BASIS2',*9999)
      FILEIP=.FALSE.
      NOQUES=0
      INFO=2
      ICHAR=999 !RGB 14/12/97

C KAT 14Dec99: Only set if IOTYPE.NE.3
CC cpb 11/9/95 Adding zero cross derivatives.
C      NBCD(nb)=0

      FORMAT='($,'' Enter the number of Xi-coordinates [2]: '',I1)'
      IDEFLT(1)=2
      IF(IOTYPE.EQ.3) IDATA(1)=NIT(nb)
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,2,3,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NIT(nb)=IDATA(1)
      CALL ASSERT(NIT(nb).LE.NIM,'>>Increase NIM',ERROR,*9999)
      IF(IOTYPE.EQ.1.AND.NIT(nb).EQ.3) THEN
        WRITE(OP_STRING,'('' >>Remember to define 2D face basis '','
     '    //'''functions'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      IF(IOTYPE.NE.3) THEN
        NUT(nb)=NIT(nb)*NIT(nb)+2
        CALL ASSERT(NUT(nb).LE.NUM,'>>Increase NUM',ERROR,*9999)
        NGT(nb)=1
        NNT(nb)=1
        DO ni=1,NIT(nb)
          NGAP(ni,nb)=1
          IBT(1,ni)=1
C news MPN 9-Jul-96 Initialise IBT
          IBT(2,ni)=0
          IBT(3,ni)=0
        ENDDO !ni
      ENDIF

C cpb 15/7/95 Generalising sectors
C      FORMAT='('' Specify the type of basis in the 1-2 plane [1]:'''//
C     '  '/''   (1) Simplex'''//
C     '  '/''   (2) Serendipity'''//
C     '  '/''   (3) Cubic sector'''//
C     '  '/$,''    '',I1)'
      FORMAT='('' Specify the type of basis [1]:'''//
     '  '/''   (1) Simplex'''//
     '  '/''   (2) Serendipity'''//
     '  '/''   (3) Sector'''//
     '  '/$,''    '',I1)'
      IF(IOTYPE.EQ.3) IDATA(1)=NBSC(1,nb)
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,3,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NBSC(1,nb)=IDATA(1)

      IF(NBSC(1,nb).LE.2) THEN
        FORMAT='('' Specify the basis degree in the 1-2 plane [1]:'''//
     '    '/''   (1) Linear'''//
     '    '/''   (2) Quadratic'''//
     '    '/''   (3) Cubic '''//
     '    '/''   (4) Hermite'''//
     '    '/$,''    '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=NBSC(2,nb)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,4,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(iotype.ne.3) NBSC(2,nb)=IDATA(1)

C CS 14/8/97 new parametric coordinate type option
        IF((NBSC(1,nb).EQ.1).AND.(NBSC(2,nb).LE.2)) THEN
          FORMAT='('' Specify the parametric coordinate type [1]:'''//
     '      '/''   (1) Area'''//
     '      '/''   (2) Xi'''//
     '      '/$,''    '',I1)'
          IF(IOTYPE.EQ.3) IDATA(1)=IBT(3,1)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,2,
     &      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(iotype.ne.3) IBT(3,1)=IDATA(1)
        ENDIF

C!!! CS 8-DEC-96: Initialising NNT(nb) for simplexes
        IF(NBSC(1,nb).EQ.1) THEN !Simplex
          IF(NBSC(2,nb).EQ.1) THEN !Linear
            NNT(nb)=NIT(nb) + 1
          ELSEIF(NBSC(2,nb).EQ.2) THEN !Quadratic
            NNT(nb)=NIT(nb)*3
            IF(NIT(nb).EQ.3) NNT(nb) = NNT(nb) + 1
          ENDIF
        ENDIF

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
            IF(IBT(1,ni).EQ.2) THEN
              IDATA(1)=4
            ELSE
              IDATA(1)=IBT(2,ni)
            ENDIF
          ENDIF
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,4,
     &      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) IBT(2,ni)=IDATA(1)
          IF(nb.GT.1) THEN
            IF(NIT(nb).EQ.NIT(nb-1).AND.
     '        (NBC(nb-1).NE.5.OR.NBC(nb-1).EQ.6)) THEN
              IDEFLT(1)=NGAP(ni,nb-1)
            ENDIF
          ELSE
            IF(IBT(2,ni).EQ.1) THEN
              IDEFLT(1)=2
            ELSE IF(IBT(2,ni).GT.1) THEN
C cpb 19/2/98 If using bicubic then use four Gauss points
              IF(NIT(nb).EQ.1) THEN
                IDEFLT(1)=3
              ELSE
                IDEFLT(1)=4
              ENDIF
            ENDIF
          ENDIF
          WRITE(CHAR2,'(I1)') IDEFLT(1)
          FORMAT='($,'' Enter the number of Gauss points in the Xi('//
     '      CHAR1(1:1)//') direction ['//CHAR2(1:1)//']: '',I1)'
          IF(IOTYPE.EQ.3) IDATA(1)=NGAP(ni,nb)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,7,
     &      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            NGAP(ni,nb)=IDATA(1)
            NGT(nb)=NGT(nb)*NGAP(ni,nb)
          ENDIF
          FORMAT='($,'' Is the Xi('//CHAR1(1:1)//') direction '
     '      //'collapsed [N]? '',A)'
          IF(IOTYPE.EQ.3) THEN
            IF((IBT(1,ni).EQ.5).OR.(IBT(1,ni).EQ.6)) THEN
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
              IF(IBT(2,ni).EQ.4) THEN
                IBT(1,ni)=2
                IBT(2,ni)=1
              ELSE
                IBT(1,ni)=1
              ENDIF
            ENDIF
          ENDIF
        ENDDO !ni
        CALL ASSERT(NGT(nb).LE.NGM,'>>Increase NGM',ERROR,*9999)
        CALL ASSERT(NUM_COLLAPSED.LE.(NIT(nb)-1),
     '    '>>Invalid collapsing of element',ERROR,*9999)
        CALL ASSERT(NUM_COLLAPSED.GT.0,
     '    '>>Use normal tensor product basis function',ERROR,*9999)
        IF(NUM_COLLAPSED.EQ.1.AND.NIT(nb).EQ.3) THEN
          IF(COLLAPSED(1).EQ.1) THEN
            IDIR(1)=2
            IDIR(2)=3
          ELSE IF(COLLAPSED(1).EQ.2) THEN
            IDIR(1)=1
            IDIR(2)=3
          ELSE
            IDIR(1)=1
            IDIR(2)=2
          ENDIF
          WRITE(CHAR1,'(I1)') IDIR(1)
          WRITE(CHAR2,'(I1)') IDIR(2)
          FORMAT='(/'' Enter the face of collapse [1]:'''//
     '      '/''   (1) Xi('//CHAR1(1:1)//')=0 face'''//
     '      '/''   (2) Xi('//CHAR1(1:1)//')=1 face'''//
     '      '/''   (3) Xi('//CHAR2(1:1)//')=0 face'''//
     '      '/''   (4) Xi('//CHAR2(1:1)//')=1 face'''//
     '      '/$,''    '',I1)'
          IF(IOTYPE.EQ.3) THEN
            IF(IBT(1,COLLAPSED(1)).EQ.5) THEN
              IF(IBT(3,COLLAPSED(1)).EQ.IDIR(1)) THEN
                IDATA(1)=1
              ELSE
                IDATA(1)=3
              ENDIF
            ELSE IF(IBT(1,COLLAPSED(1)).EQ.6) THEN
              IF(IBT(3,COLLAPSED(1)).EQ.IDIR(1)) THEN
                IDATA(1)=2
              ELSE
                IDATA(1)=4
              ENDIF
            ENDIF
          ENDIF
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,4,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            IF(IDATA(1).EQ.1) THEN
              IBT(1,COLLAPSED(1))=5
              IBT(3,COLLAPSED(1))=IDIR(1)
            ELSE IF(IDATA(1).EQ.2) THEN
              IBT(1,COLLAPSED(1))=6
              IBT(3,COLLAPSED(1))=IDIR(1)
            ELSE IF(IDATA(1).EQ.3) THEN
              IBT(1,COLLAPSED(1))=5
              IBT(3,COLLAPSED(1))=IDIR(2)
            ELSE
              IBT(1,COLLAPSED(1))=6
              IBT(3,COLLAPSED(1))=IDIR(2)
            ENDIF
C cpb 19/7/95 If the interpolation is Hermite along the non-collapsed
C direction the adjust the Xi interpolation type to be the appropriate
C quadratic Hermite (i.e. cannot have a derivative at a collapsed node)
            IF(IBT(1,COLLAPSED(1)).EQ.5) THEN
              IF(IBT(1,IBT(3,COLLAPSED(1))).EQ.2)
     '          IBT(2,IBT(3,COLLAPSED(1)))=2
            ELSE
              IF(IBT(1,IBT(3,COLLAPSED(1))).EQ.2)
     '          IBT(2,IBT(3,COLLAPSED(1)))=3
            ENDIF
          ENDIF
        ELSE
          IF(NIT(nb).EQ.2) THEN
            IF(COLLAPSED(1).EQ.1) THEN
              IBT(3,COLLAPSED(1))=2
            ELSE
              IBT(3,COLLAPSED(1))=1
            ENDIF
          ELSE
            IF(COLLAPSED(1).EQ.1.AND.COLLAPSED(2).EQ.2) THEN
              IBT(3,COLLAPSED(1))=3
              IBT(3,COLLAPSED(2))=3
            ELSE IF(COLLAPSED(1).EQ.1.AND.COLLAPSED(2).EQ.3) THEN
              IBT(3,COLLAPSED(1))=2
              IBT(3,COLLAPSED(2))=2
            ELSE IF(COLLAPSED(1).EQ.2.AND.COLLAPSED(2).EQ.3) THEN
              IBT(3,COLLAPSED(1))=1
              IBT(3,COLLAPSED(2))=1
            ENDIF
          ENDIF
          WRITE(CHAR1,'(I1)') IBT(3,COLLAPSED(1))
          FORMAT='(/'' Specify the end of collapse [1]:'''//
     '      '/''   (1) Xi=0 end of Xi('//CHAR1(1:1)//')'''//
     '      '/''   (2) Xi=1 end of Xi('//CHAR1(1:1)//')'''//
     '      '/$,''    '',I1)'
          IF(IOTYPE.EQ.3) THEN
            IF(IBT(1,COLLAPSED(1)).EQ.5) THEN
              IDATA(1)=1
            ELSE
              IDATA(1)=2
            ENDIF
          ENDIF
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,2,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            IF(IDATA(1).EQ.1) THEN
              DO i=1,NUM_COLLAPSED
                IBT(1,COLLAPSED(i))=5
              ENDDO !i
              IF(IBT(1,IBT(3,COLLAPSED(1))).EQ.2)
     '          IBT(2,IBT(3,COLLAPSED(1)))=2
            ELSE
              DO i=1,NUM_COLLAPSED
                IBT(1,COLLAPSED(i))=6
              ENDDO !i
              IF(IBT(1,IBT(3,COLLAPSED(1))).EQ.2)
     '          IBT(2,IBT(3,COLLAPSED(1)))=3
            ENDIF
          ENDIF
        ENDIF
      ENDIF

C cpb 11/9/95 Adding zero cross derivatives.
      NUMHERMXI=0
      DO ni=1,NIT(nb)
        IF(((IBT(1,ni).EQ.5.OR.IBT(1,ni).EQ.6).AND.IBT(2,ni).EQ.4).OR.
     '    IBT(1,ni).EQ.2) THEN
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
            INP(nn,1)=ni1
            INP(nn,2)=ni2
            IF(NIT(nb).EQ.3) INP(nn,3)=1
          ENDDO
        ENDDO

      ELSE IF(NBSC(1,nb).EQ.2) THEN
        nn=0
        CALL ASSERT(NBSC(2,nb).LT.4,'>>Hermite serendipity not set up',
     '    ERROR,*9999)
        DO ni2=1,NI2_TOP+1
          IF(ni2.GE.2.AND.ni2.LE.NBSC(2,nb)) THEN
            INCR=NBSC(2,nb)
          ELSE
            INCR=1
          ENDIF
          DO ni1=1,NBSC(2,1)+1,INCR
            nn=nn+1
            INP(nn,1)=ni1
            INP(nn,2)=ni2
            IF(NIT(nb).EQ.3) INP(nn,3)=1
          ENDDO
        ENDDO
        NNT(nb)=nn

      ELSE IF(NBSC(1,nb).EQ.3) THEN

C LC 24/2/97 LC archived section : cpb 15/7/95 Generalising sectors
C       Determine the maximum number of nodes
        nnmax=1
        DO ni=1,NIT(nb)
          IF((IBT(1,ni).EQ.2).OR.(IBT(1,ni).EQ.5.AND.IBT(2,ni).EQ.4)
     '      .OR.(IBT(1,ni).EQ.6.AND.IBT(2,ni).EQ.4)) THEN
            numnn=2
          ELSE
            numnn=IBT(2,ni)+1
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
            IF((IBT(1,ni).EQ.5.AND.POSITION(IBT(3,COLLAPSED(1))).EQ.1)
     '        .OR.(IBT(1,ni).EQ.6.AND.POSITION(IBT(3,COLLAPSED(1))).EQ.
     '        NUMNODES(IBT(3,COLLAPSED(1))))) ATCOLLAPSE=.TRUE.
          ENDDO !ni
          IF(ATCOLLAPSE) THEN
            IF(POSITION(COLLAPSED(1)).EQ.1) THEN
              nn=nn+1
              IF(nn.LE.NNM) THEN
                DO ni=1,NIT(nb)
                  INP(nn,ni)=POSITION(ni)
                ENDDO !ni
              ENDIF
            ENDIF
          ELSE
            nn=nn+1
            IF(nn.LE.NNM) THEN
              DO ni=1,NIT(nb)
                INP(nn,ni)=POSITION(ni)
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
        IBT(1,1)=NBSC(1,nb)+2
        IBT(1,2)=NBSC(1,nb)+2
        IBT(2,1)=NBSC(2,nb)
        IBT(2,2)=NBSC(2,nb)
      ENDIF
      CALL ASSERT(NNT(nb).LE.NNM,'>>Need to increase NNM',ERROR,*9999)

      IF(NBSC(2,nb).NE.4.AND.NBSC(1,nb).NE.3) THEN !Not a Herm Simplex or Sector
        FORMAT='('' The number of Gauss points in the 1-2 '
     '    //'plane is [4]:'''//
     '    '/''   ( if Simplex [1,4 or 7] )'''//
     '    '/''   ( if Serendipity [1,4 or 9] )'''//
     '    '/$,''   '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=NGAP(1,nb)
        IDEFLT(1) = 4
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,9,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(iotype.ne.3) NGAP(1,nb)=IDATA(1)

        NGAP(2,nb)=1

        IF(NIT(nb).EQ.3) THEN
          FORMAT='('' Specify the type of basis in the 3-dir.n [1]:'''//
     '      '/''   (1) Lagrange'''//
     '      '/''   (2) Simplex'''//
     '      '/''   (3) Serendipity'''//
     '      '/$,''    '',I1)'
          IF(IOTYPE.EQ.3) IDATA(1)=NBSC(3,nb)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,3,LDATA,
     &      LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) NBSC(3,nb)=IDATA(1)

          FORMAT='(/'' Specify the basis degree in the 3-dir.n [1]:'''//
     '      '/''   (1) Linear'''//
     '      '/''   (2) Quadratic'''//
     '      '/''   (3) Cubic'''//
     '      '/''   (4) Cubic Hermite'''//
     '      '/$,''    '',I1)'
          IF(IOTYPE.EQ.3) IDATA(1)=NBSC(4,nb)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,4,
     &      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(iotype.ne.3) NBSC(4,nb)=IDATA(1)

          IF(NBSC(1,nb).ne.3) THEN
            CALL ASSERT(NBSC(4,nb).LT.4,'>>Hermite only for sector.',
     '        ERROR,*9999)
          ENDIF

          nn=NNT(nb)
          IF(NBSC(3,nb).EQ.1) THEN   ! Lagrange/Hermite in Xi3
            IF(NBSC(4,nb).LE.3) THEN ! Lgrnge/Hrmte,smplx,srndpty in Xi3
              DO ni3=2,NBSC(4,nb)+1
                DO mn=1,NNT(nb)
                  nn=nn+1
                  INP(nn,1)=INP(mn,1)
                  INP(nn,2)=INP(mn,2)
                  INP(nn,3)=ni3
                ENDDO !mn
              ENDDO ! ni3
            ELSE                     ! Hermite in Xi3
              ni3=2
              DO mn=1,NNT(nb)
                nn=nn+1
                INP(nn,1)=INP(mn,1)
                INP(nn,2)=INP(mn,2)
                INP(nn,3)=ni3
              ENDDO
            ENDIF
            NNT(nb)=nn
          ELSE IF(NBSC(3,nb).EQ.2) THEN       ! Simplex in Xi3
            nn = nn - 1 !TEMPORARY RGB - fixing simplexes in 3D
            DO ni3=2,NBSC(4,nb)+1
              DO ni2=1,NBSC(2,nb)-ni3+2
                INCR=1
                IF(NBSC(1,nb).EQ.1) THEN      ! Simplex in Xi1-2
                  LIMIT=NBSC(2,nb)+3-ni2-ni3
                ELSE IF(NBSC(1,nb).EQ.2) THEN ! Serendipity in Xi1-2
                  LIMIT=NBSC(2,nb)+2-ni3
                  IF(ni2.GE.2.AND.ni2.LE.NBSC(2,nb)+1-ni3) THEN
                    INCR=NBSC(2,nb)+1-ni3
                  ENDIF
                ELSE IF(NBSC(1,nb).EQ.3) THEN ! Sector in Xi1-2
                ENDIF
                DO ni1=1,LIMIT,INCR
                  nn=nn+1
                  INP(nn,1)=ni1
                  INP(nn,2)=ni2
                  INP(nn,3)=ni3
                ENDDO ! ni1
              ENDDO ! ni2
            ENDDO ! ni3
          ELSE IF(NBSC(3,nb).EQ.3) THEN       ! Serendipity in Xi3
          ENDIF

          IBT(2,3)=NBSC(4,nb)
          IF(NBSC(3,nb).EQ.1) THEN
            IF(NBSC(4,nb).LE.3) THEN
              IBT(1,3)=1
            ELSE !Hermite
              IBT(1,3)=2
              IBT(2,3)=1
            ENDIF
          ELSE IF(NBSC(3,nb).GT.1) THEN
            IBT(1,3)=NBSC(3,nb)+1
          ENDIF
          FORMAT='($,'' The number of Gauss points in the'//
     '      ' 3-direction is [1]: '',I1)'
          IF(IOTYPE.EQ.3) IDATA(1)=NGAP(2,nb)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,9,
     &      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(iotype.ne.3) NGAP(2,nb)=IDATA(1)

        ENDIF
        NST(nb)=NNT(nb)
        IF(NBSC(1,nb).LE.2) THEN !Simplex/Serendipity
          NGT(nb)=NGAP(1,nb)*NGAP(2,nb)
        ELSE IF(NBSC(1,nb).EQ.3) THEN !Sector
          NGT(nb)=NGAP(1,nb)*NGAP(2,nb)
        ENDIF
        CALL ASSERT(NGT(nb).LE.NGM,' >>Need to increase NGM',
     '    ERROR,*9999)

C LC 24/2/97 archived section : PJH   9Aug91 - use default INP

        NKT(0,nb)=1
        DO nn=1,NNT(nb)
          NKT(nn,nb)=1
          DO ni=1,NIT(nb)
            IDO(1,nn,ni)=1
          ENDDO !ni
        ENDDO !nn

C LC 24/2/97 archived section except call line following:
        CALL GAUSS2(IBT,INP,nb,NGAP(1,nb),PG,WG,XIG,ERROR,*9999)

      ELSE IF(NBSC(1,nb).EQ.3) THEN !Sector

C KAT 3May98: modified to allow deriv in `non-sector' direction at
C       collapsed nodes.

        NST(nb)=0
        NKT(0,nb)=1
        DO nn=1,NNT(nb)
          NKT(nn,nb)=1
          ATCOLLAPSE=ISATCOLLAPSE(IBT,INP,nb,nn)
          DO ni=1,NIT(nb)
C cpb 02/02/00 Bug fix
C            IF((.NOT.ATCOLLAPSE.OR.IBT(3,ni).NE.0) !not collapsed in ni dirn
            IF((.NOT.ATCOLLAPSE.OR.IBT(3,ni).EQ.0) !not collapsed in ni dirn
     '        .AND.((IBT(1,ni).EQ.2.AND.(IBT(2,ni).EQ.1
     '        .OR.INP(nn,ni).NE.(IBT(2,ni)-1)))
     '        .OR.(IBT(1,ni).EQ.5.AND.IBT(2,ni).EQ.4)
     '        .OR.(IBT(1,ni).EQ.6.AND.IBT(2,ni).EQ.4))) THEN !deriv in ni dirn
C             Deriv in this direction
              DO nk=1,NKT(nn,nb)
                IDO(nk,nn,ni)=1
                IF(nk+NKT(nn,nb).LE.NKM) THEN
                  DO ni2=1,ni-1
                    IDO(nk+NKT(nn,nb),nn,ni2)=IDO(nk,nn,ni2)
                  ENDDO !ni2
                ENDIF
              ENDDO !nk
              IF(NKT(nn,nb)*2.LE.NKM) THEN
                DO nk=NKT(nn,nb)+1,NKT(nn,nb)*2
                  IDO(nk,nn,ni)=2
                ENDDO !nk
              ENDIF
              NKT(nn,nb)=NKT(nn,nb)*2
            ELSE
              IF(NKT(nn,nb).LE.NKM) THEN
                DO nk=1,NKT(nn,nb)
                  IDO(nk,nn,ni)=1
                ENDDO !nk
              ENDIF
            ENDIF
          ENDDO !ni
C cpb 11/9/95 Adding zero cross derivatives. This is not the best way
C to do this but it is only temporary.
          IF(NBCD(nb).EQ.1) THEN
            IF(NKT(nn,nb).EQ.4) THEN
              NKT(nn,nb)=3
              DO ni=1,NIT(nb)
                IDO(4,nn,ni)=0
              ENDDO !ni
            ENDIF
          ENDIF
          IF(NKT(nn,nb).GT.NKT(0,nb)) NKT(0,nb)=NKT(nn,nb)
          NST(nb)=NST(nb)+NKT(nn,nb)
        ENDDO !nn
        CALL ASSERT(NKT(0,nb).LE.NKM,'>>Increase NKM',ERROR,*9999)
        CALL ASSERT(NST(nb).LE.NSM,'>>Increase NSM',ERROR,*9999)

        CALL GAUSS1(IBT,IDO,INP,nb,NGAP(1,nb),PG,WG,XIG,ERROR,*9999)

      ELSE !Hermite simplex

C       Only set up for 2 dimensions initially
        CALL ASSERT(NIT(nb).LE.2,' >>Hermite simplex not implemented'
     '    //' for 3 dimensions',ERROR,*9999)

        IDEFLT(1)=3
        DO ni=1,NIT(nb)
          WRITE(CHAR1,'(I1)') ni
          WRITE(CHAR2,'(I1)') IDEFLT(1)
          FORMAT='($,'' Enter the number of Gauss points in the Xi('//
     '      CHAR1(1:1)//') direction ['//CHAR2(1:1)//']: '',I1)'
          IF(IOTYPE.EQ.3) IDATA(1)=NGAP(ni,nb)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,7,
     &      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(iotype.ne.3) THEN
            NGAP(ni,nb)=IDATA(1)
            NGT(nb)=NGT(nb)*NGAP(ni,nb)
            CALL ASSERT(NGT(nb).LE.NGM,'>>Need to increase NGM',
     '        ERROR,*9999)
            IBT(1,ni)=3
            IBT(2,ni)=4 !Special hermite
            NNT(nb)=3
          ENDIF
        ENDDO
        FORMAT='($,'' Is the apex node local node 1 or 3 [1]: '',I1)'
        IF(IOTYPE.EQ.3) THEN
          nn=1
          DO WHILE((NKT(nn,nb).GT.1).AND.(nn.LE.NNT(nb)))
            nn=nn+1
          ENDDO
          CALL ASSERT(nn.LE.NNT(nb),'>>Error in Hermite simplex write',
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
!         !Set up standard IDO
          DO nn=1,NNT(nb)
            nk=0
            DO ni2=1,2
              DO ni1=1,2
                nk=nk+1
                IDO(nk,nn,1)=ni1
                IDO(nk,nn,2)=ni2
              ENDDO
            ENDDO
            IDO(1,nn,0)=1 !Standard mapping for zeroth component
            IDO(2,nn,0)=2
            IDO(3,nn,0)=4
            IDO(4,nn,0)=6 !GMH 2/11/95 cross derivative
!            IDO(4,nn,0)=7
          ENDDO !nn
        ENDIF
        IF(IDATA(1).ne.3)THEN
!         !Need to correct INP for apex at node 1.
          INP(1,1)=1
          INP(1,2)=1
          INP(2,1)=1
          INP(2,2)=2
          INP(3,1)=2
          INP(3,2)=2
        ENDIF

        CALL GAUSS2_HERMITE(IDO,INP,nb,NGAP(1,nb),PG,WG,XIG,
     '    ERROR,*9999)

      ENDIF !End of Hermite Simplex loop

      CALL EXITS('BASIS2')
      RETURN
 9999 CALL ERRORS('BASIS2',ERROR)
      CALL EXITS('BASIS2')
      RETURN 1
      END


