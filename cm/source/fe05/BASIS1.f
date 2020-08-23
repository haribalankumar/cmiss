      SUBROUTINE BASIS1(IBT,IDO,INP,nb,NGAP,PG,WG,XIG,ERROR,*)

C#### Subroutine: BASIS1
C###  Description:
C###    BASIS1 inputs parameters for Lagrange/Hermite tensor
C###    product basis functions and Gauss-Legendre or Gauss-Lobatto
C###    quadrature.

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
      INTEGER IB,IB1,IBEG,ICHAR,IDO1,IDO2,IDO3,IEND,IMAX_L,INFO,k,ng,
     '  ni,ni2,nk,NKT_OLD,nn,NOQUES,NUMHERMXI,NUMNODES(3),POSITION(4)
      CHARACTER CHAR*1,CHAR1*192,CHAR2*1,CHAR3*2
      LOGICAL DIFFIDO,FILEIP
!     External functions
      INTEGER IDIGITS

C CHAR1 must be 192 to account for tri-linear lagrange (64*3)

      CALL ENTERS('BASIS1',*9999)
      FILEIP=.FALSE.
      NOQUES=0
      INFO=2
      ICHAR=999 !RGB 14/12/97

C KAT 14Dec99: Only set if IOTYPE.NE.3
CC cpb 11/9/95 Adding zero cross derivatives.
C      NBCD(nb)=0

      DO k=1,100
        CHAR1(k:k)=' '
      ENDDO
      FORMAT='($,'' Enter the number of Xi-coordinates [1]: '',i1)'
      IF(IOTYPE.EQ.3) IDATA(1)=NIT(nb)
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NJT,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NIT(nb)=IDATA(1)
C     CALL ASSERT(NIT(nb).LE.NIM,'>>Increase NIM',ERROR,*9999)
      IF(NIT(nb).GT.NIM) THEN
        WRITE(CHAR,'(I1)') IDIGITS(NIT(nb))
        WRITE(ERROR,'(''>>Increase NIM to '',I'//CHAR//')') NIT(nb)
        GO TO 9999
      ENDIF
      IF(IOTYPE.EQ.1.AND.NIT(nb).EQ.3) THEN
        WRITE(OP_STRING,'('' >>Remember to define 2D face basis '','
     '    //'''functions'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      IF(IOTYPE.NE.3) THEN
        NUT(nb)=NIT(nb)*NIT(nb)+2
C       CALL ASSERT(NUT(nb).LE.NUM,'>>Increase NUM',ERROR,*9999)
        IF(NUT(nb).GT.NUM) THEN
          WRITE(CHAR,'(I1)') IDIGITS(NUT(nb))
          WRITE(ERROR,'(''>>Increase NUM to '',I'//CHAR//')') NUT(nb)
          GO TO 9999
        ENDIF
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

      NOQUES=0

      NUMHERMXI=0
      DO ni=1,NIT(nb)
        IB1=0
        WRITE(CHAR1,'(I1)') ni
C CPB 13/7/95 Adding quadratic Hermite interpolation
        FORMAT='(/'' The interpolant in the Xi('//
     '    CHAR1(1:1)//') direction is [1]: '''//
     '    '/''   (1) Linear Lagrange'''//
     '    '/''   (2) Quadratic Lagrange'''//
     '    '/''   (3) Cubic Lagrange'''//
     '    '/''   (4) Quadratic Hermite'''//
     '    '/''   (5) Cubic Hermite'''//
     '    '/$,''    '',I1)'
        IF(IOTYPE.EQ.3) THEN
          IF(IBT(1,ni).EQ.1) THEN
            IB=IBT(2,ni)
          ELSE IF(IBT(1,ni).EQ.2) THEN
            IF(IBT(2,ni).EQ.1) THEN
              IB=5
            ELSE IF((IBT(2,ni).EQ.2).OR.(IBT(2,ni).EQ.3)) THEN
              IB=4
            ENDIF
          ENDIF
          IDATA(1)=IB
        ENDIF
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,5,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(iotype.ne.3) IB=IDATA(1)

C MPN 16-Nov-95: adding tri-cubic Hermite interpolation
CC MPN 17/7/95 Flaging warning if try to use tri-cubic Hermite fns
C        IF(IB.NE.5) ALLCUBHERM=.FALSE.
C        IF(ni.EQ.3.AND.ALLCUBHERM) THEN !user chose tri-cubic Hermite
C          ERROR='>>Tri-cubic Hermite basis fns not implemented'
C          GO TO 9999
C        ENDIF

c cpb 13/7/95 Adding quadratic Hermite interpolation
        IF(IB.EQ.4) THEN !Quadratic Hermite
          FORMAT='($,'' Enter the local node # with no derivative '
     '      //'(1 or 2) [1]: '',I1)'
          IF(IOTYPE.EQ.3) IDATA(1)=IBT(2,ni)-1
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,2,
     &      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) IB1=IDATA(1)
        ENDIF

        IF(nb.GT.1) THEN
          IF(NIT(nb).EQ.NIT(nb-1).AND.
     '      (NBC(nb-1).ne.5.OR.NBC(nb-1).EQ.6)) THEN
            IDEFLT(1)=NGAP(ni,nb-1)
          ELSE
            IF(IB.EQ.1) THEN
              IDEFLT(1)=2
            ELSE IF(IB.GT.1) THEN
              IDEFLT(1)=3
            ENDIF
          ENDIF
        ELSE
          IF(IB.EQ.1) THEN
            IDEFLT(1)=2
          ELSE IF(IB.GT.1) THEN
C cpb 28/9/98 Just use 3 Gauss points and trade of accuracy of the
C elemental stiffness integrals with amount of computational work
CC cpb 19/2/98 If using bicubic then use four Gauss points
C            IF(NIT(nb).EQ.1) THEN
C              IDEFLT(1)=3
C            ELSE
C              IDEFLT(1)=4
C            ENDIF
            IDEFLT(1)=3
          ENDIF
        ENDIF
        WRITE(CHAR2,'(I1)') IDEFLT(1)
        FORMAT='($,'' Enter the number of Gauss points in the Xi('//
     '    CHAR1(1:1)//') direction ['//CHAR2(1:1)//']: '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=NGAP(ni,nb)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,7,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          NGAP(ni,nb)=IDATA(1)
          NGT(nb)=NGT(nb)*NGAP(ni,nb)
          IF(IB.LE.3) THEN !Lagrange
            IBT(1,ni)=1
            IBT(2,ni)=IB
            NUMNODES(ni)=IB+1
            NNT(nb)=NNT(nb)*(IB+1)
          ELSE IF(IB.EQ.4) THEN !quadratic Hermite
            IBT(1,ni)=2
            IBT(2,ni)=1+IB1
            NUMNODES(ni)=2
            NNT(nb)=NNT(nb)*2
          ELSE !cubic Hermite
            IBT(1,ni)=2
            IBT(2,ni)=1
            NUMNODES(ni)=2
            NNT(nb)=NNT(nb)*2
          ENDIF
        ENDIF
        IF(IB.GE.4) !quadratic/cubic Hermite
     '    NUMHERMXI=NUMHERMXI+1
        POSITION(ni)=1
      ENDDO !ni
C     CALL ASSERT(NGT(nb).LE.NGM,'>>Increase NGM',ERROR,*9999)
      IF(NGT(nb).GT.NGM) THEN
        WRITE(CHAR,'(I1)') IDIGITS(NGT(nb))
        WRITE(ERROR,'(''>>Increase NGM to '',I'//CHAR//')') NGT(nb)
        GO TO 9999
      ENDIF
C     CALL ASSERT(NNT(nb).LE.NNM,'>>Increase NNM',ERROR,*9999)
      IF(NNT(nb).GT.NNM) THEN
        WRITE(CHAR,'(I1)') IDIGITS(NNT(nb))
        WRITE(ERROR,'(''>>Increase NNM to '',I'//CHAR//')') NNT(nb)
        GO TO 9999
      ENDIF
      
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
C 24/2/97 LC archived section :
C      cpb 13/7/95 Adding quadratic Hermite interpolation and this now needs
C 24/2/97 LC archived section :
C       cpb 13/7/95 Need to calc NKT (and hence NST) after INP

      IF(IOTYPE.NE.3) THEN
C       Determine the node positions (INP)
        DO nn=1,NNT(nb)
          ni=1
          DO WHILE(POSITION(ni).GT.NUMNODES(ni))
            POSITION(ni)=1
            ni=ni+1
            POSITION(ni)=POSITION(ni)+1
          ENDDO !ni
          DO ni=1,NIT(nb)
            INP(nn,ni)=POSITION(ni)
          ENDDO !ni
          POSITION(1)=POSITION(1)+1
        ENDDO !nn
      ENDIF

      k=0
      DO nn=1,NNT(nb)
        DO ni=1,NIT(nb)
          k=k+1
          IDEFLT(k)=INP(nn,ni)
          WRITE(CHAR1(k:k),'(I1)') IDEFLT(k)
        ENDDO
      ENDDO
C cpb 21/3/96 There is a problem when we have a large number of local
C nodes. In this case the INP string is very large and runs of the end
C of the record. Hence if the number of local nodes is over 27 (which
C corresponds to a tri-quadratic Lagrange basis prompt the user for
C the indicies for each local node
      IF(K.LE.27) THEN
        CALL STRING_TRIM(CHAR1,IBEG,IEND)
        FORMAT='($,'' Enter the node position indices ['//
     '    CHAR1(IBEG:IEND)//']: '',27I2)'
        IMAX_L=NNT(nb)
        IF(IOTYPE.EQ.3) THEN
          k=0
C          IMAX_L=10
          DO nn=1,NNT(nb)
            DO ni=1,NIT(nb)
              k=k+1
              IDATA(k)=INP(nn,ni)
C              IMAX_L=IMAX_L*10
            ENDDO
          ENDDO
        ENDIF
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,k,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,IMAX_L,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.EQ.3) THEN
          k=0
          DO nn=1,NNT(nb)
            DO ni=1,NIT(nb)
              k=k+1
              INP(nn,ni)=IDATA(k)
            ENDDO
          ENDDO
        ENDIF
      ELSE
        CHAR1=' ' !rgb - initialising string
        DO nn=1,NNT(nb)
          WRITE(CHAR3,'(I2)') nn
          DO ni=1,NIT(nb)
            IDEFLT(ni)=INP(nn,ni)
            WRITE(CHAR1(ni:ni),'(I1)') IDEFLT(ni)
          ENDDO !ni
          CALL STRING_TRIM(CHAR1,IBEG,IEND)
          FORMAT='($,'' Enter the node position indices for local '
     '      //'node '//CHAR3//' ['//CHAR1(IBEG:IEND)//']: '',27I2)'
          IF(IOTYPE.EQ.3) THEN
            DO ni=1,NIT(nb)
              IDATA(ni)=INP(nn,ni)
            ENDDO !ni
          ENDIF
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     '      NIT(nb),ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,4,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            DO ni=1,NIT(nb)
              INP(nn,ni)=IDATA(ni)
            ENDDO !ni
          ENDIF
        ENDDO !nn
      ENDIF

C 24/2/97 LC archived section : cpb 14/7/95 Old NKT,IDO calculation
C cpb 14/7/95 New NKT,NST,IDO calculations for different numbers of
C derivatives at each local node (e.g. quadratic Hermite)

      IF(IOTYPE.NE.3) THEN
        NST(nb)=0
        NKT(0,nb)=1
        DO nn=1,NNT(nb)
          NKT(nn,nb)=1
          DO ni=1,NIT(nb)
            IF(IBT(1,ni).EQ.2.AND.
     '        (IBT(2,ni).EQ.1.OR.INP(nn,ni).NE.IBT(2,ni)-1)) THEN
              NKT_OLD=NKT(nn,nb)
              NKT(nn,nb)=NKT(nn,nb)*2
              IF(NKT(nn,nb).LE.NKM) THEN
                DO nk=1,NKT_OLD
                  IDO(nk,nn,ni)=1
                  IDO(NKT_OLD+nk,nn,ni)=2
                  DO ni2=1,ni-1
                    IDO(NKT_OLD+nk,nn,ni2)=IDO(nk,nn,ni2)
                  ENDDO !ni2
                ENDDO !nk
              ENDIF !NKM
            ELSE
              IF(NKT(nn,nb).LE.NKM) THEN
                DO nk=1,NKT(nn,nb)
                  IDO(nk,nn,ni)=1
                ENDDO !nk
              ENDIF !NKM
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
      ENDIF

      IF(NKT(0,nb).GT.NKM) THEN
        IEND=0
        CALL APPENDC(IEND,'Increase NKM to ',ERROR)
        CALL APPENDI(IEND,NKT(0,nb),ERROR)
        GOTO 9999
      ENDIF

      IF(NKT(0,nb).GT.1) THEN
        DIFFIDO=.FALSE.
        DO nn=1,NNT(nb)
          IF(NKT(nn,nb).NE.NKT(0,nb)) DIFFIDO=.TRUE.
        ENDDO !nn
        IF(DIFFIDO) THEN
          DO nn=1,NNT(nb)
            CHAR1=' '
            k=0
            DO nk=1,NKT(nn,nb)
              DO ni=1,NIT(nb)
                k=k+1
                IDEFLT(k)=IDO(nk,nn,ni)
                WRITE(CHAR1(k:k),'(I1)') IDEFLT(k)
              ENDDO !ni
            ENDDO !nk
            CALL STRING_TRIM(CHAR1,IBEG,IEND)
            WRITE(CHAR2,'(I1)') nn
            FORMAT='($,'' Enter the derivative order indices for '
     '        //'local node '//CHAR2(1:1)//' ['//CHAR1(IBEG:IEND)
     '        //']: '',40I2)'
            IMAX_L=NKT(0,nb)
            IF(IOTYPE.EQ.3) THEN
              k=0
C              IMAX_L=10
              DO nk=1,NKT(nn,nb)
                DO ni=1,NIT(nb)
                  k=k+1
                  IDATA(k)=IDO(nk,nn,ni)
C                  IMAX_L=10
                ENDDO !ni
              ENDDO !nk
            ENDIF
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,k,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     &        IMAX_L,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) THEN
              k=0
              DO nk=1,NKT(NN,nb)
                DO ni=1,NIT(nb)
                  k=k+1
                  IDO(nk,NN,ni)=IDATA(k)
                ENDDO !ni
              ENDDO !nk
            ENDIF
          ENDDO !nn
        ELSE
          CHAR1=' '
          k=0
          DO nk=1,NKT(0,nb)
            DO ni=1,NIT(nb)
              k=k+1
              IDEFLT(k)=IDO(nk,1,ni)
              WRITE(CHAR1(k:k),'(I1)') IDEFLT(k)
            ENDDO !ni
          ENDDO !nk
          CALL STRING_TRIM(CHAR1,IBEG,IEND)
          FORMAT='($,'' Enter the derivative order indices ['//
     '      CHAR1(IBEG:IEND)//']: '',40I2)'
          IMAX_L=NKT(0,nb)
          IF(IOTYPE.EQ.3) THEN
            k=0
C            IMAX_L=10
            DO nk=1,NKT(0,nb)
              DO ni=1,NIT(nb)
                k=k+1
                IDATA(k)=IDO(nk,1,ni)
C                IMAX_L=IMAX_L*10
              ENDDO !ni
            ENDDO !nk
          ENDIF
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      k,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,IMAX_L,
     &      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            k=0
            DO nk=1,NKT(0,nb)
              DO ni=1,NIT(nb)
                k=k+1
                IDO(nk,1,ni)=IDATA(k)
              ENDDO !ni
            ENDDO !nk
            DO nn=2,NNT(nb)
              DO nk=1,NKT(0,nb)
                DO ni=1,NIT(nb)
                  k=k+1
                  IDO(nk,nn,ni)=IDO(nk,1,ni)
                ENDDO !ni
              ENDDO !nk
            ENDDO !nn
          ENDIF
        ENDIF
      ENDIF

      IF(IOTYPE.NE.3) THEN
        DO nn=1,NNT(nb)
          DO nk=1,NKT(nn,nb)
            IDO1=IDO(nk,nn,1)
            IDO2=1
            IDO3=1
            IF(NIT(nb).GE.2) IDO2=IDO(nk,nn,2)
            IF(NIT(nb).EQ.3) IDO3=IDO(nk,nn,3)
            IF(IDO1.EQ.1.AND.IDO2.EQ.1.AND.IDO3.EQ.1) IDO(nk,nn,0)=1
            IF(IDO1.EQ.2.AND.IDO2.EQ.1.AND.IDO3.EQ.1) IDO(nk,nn,0)=2
            IF(IDO1.EQ.1.AND.IDO2.EQ.2.AND.IDO3.EQ.1) IDO(nk,nn,0)=4
            IF(IDO1.EQ.2.AND.IDO2.EQ.2.AND.IDO3.EQ.1) IDO(nk,nn,0)=6
            IF(IDO1.EQ.1.AND.IDO2.EQ.1.AND.IDO3.EQ.2) IDO(nk,nn,0)=7
            IF(IDO1.EQ.2.AND.IDO2.EQ.1.AND.IDO3.EQ.2) IDO(nk,nn,0)=9
            IF(IDO1.EQ.1.AND.IDO2.EQ.2.AND.IDO3.EQ.2) IDO(nk,nn,0)=10
            IF(IDO1.EQ.2.AND.IDO2.EQ.2.AND.IDO3.EQ.2) IDO(nk,nn,0)=11
          ENDDO !nk
        ENDDO !nn

        CALL GAUSS1(IBT,IDO,INP,nb,NGAP(1,nb),PG,WG,XIG,ERROR,*9999)

      ENDIF

      IF(DOP) THEN
        DO ng=1,NGT(nb)
          WRITE(OP_STRING,'('' XIG(ni,'',I2,''): '',3D11.3)')
     '       ng,(XIG(ni,ng),ni=1,NIT(nb))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

      CALL EXITS('BASIS1')
      RETURN
 9999 CALL ERRORS('BASIS1',ERROR)
      CALL EXITS('BASIS1')
      RETURN 1
      END


