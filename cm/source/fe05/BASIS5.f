      SUBROUTINE BASIS5(IBT,IDO,INP,nb,NDET,NGAP,DET,PG,WG,XIG,ERROR,*)

C#### Subroutine: BASIS5
C###  Description:
C###    BASIS5 is the Boundary Element basis function routine.
C###    Input of parameters for Lagrange/Hermite tensor product basis
C###    functions and Gauss-Legendre quadrature.
C**** INP(nn,ni,nb) gives the index for element node nn in each Xi-dir.
C**** Thus: INP(nn,1..,nb) = 1,2,2 indicates that node nn is the first
C**** node in the Xi(1) dir and second in the Xi(2) & Xi(3) directions.
C**** The basis number, nb, passed to this routine is interpreted as
C**** the parent basis number of the boundary element family.

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
      INTEGER IB,IBEG,ICHAR,IDO1,IDO2,IDO3,IEND,INFO,k,n1,n2,
     '  NBTOP,ni,ni1,ni2,nk,nn,NOQUES,NUMHERMXI
      CHARACTER CHAR1*100,CHAR2*3,CHAR3*3
      LOGICAL FILEIP

      CALL ENTERS('BASIS5',*9999)
      FILEIP=.FALSE.
      NOQUES=0
      INFO=2
      ICHAR=999

      DO k=1,100
        CHAR1(k:k)=' '
      ENDDO
      FORMAT='($,'' Enter the number of Xi-coordinates [1]: '',I1)'
      IF(IOTYPE.EQ.3) IDATA(1)=NIT(nb)
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NJT,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(iotype.ne.3) NIT(nb)=IDATA(1)

      IF(IOTYPE.NE.3) THEN
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

      NOQUES=0

      NUMHERMXI=0
      DO ni=1,NIT(nb)
        WRITE(CHAR1,'(I1)') ni
        FORMAT='(/'' The interpolant in the Xi('//
     '    CHAR1(1:1)//') direction is [1]: '''//
     '    '/''   (1) Linear Lagrange'''//
     '    '/''   (2) Quadr. Lagrange'''//
     '    '/''   (3) Cubic  Lagrange'''//
     '    '/''   (4) Cubic  Hermite'''//
     '    '/$,''    '',I1)'
        IF(IOTYPE.EQ.3) THEN
          IF(IBT(1,ni,nb).EQ.1) THEN
            IB=IBT(2,ni,nb)
          ELSE IF(IBT(1,ni,nb).EQ.2 ) THEN
            IB=4
          ENDIF
          IDATA(1)=IB
        ENDIF
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,4,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(iotype.ne.3) IB=IDATA(1)

        IF(IB.EQ.1) THEN
          IDEFLT(1)=3
          IDEFLT(2)=6
        ELSE IF(IB.EQ.2.OR.IB.EQ.3) THEN
          IDEFLT(1)=3
          IDEFLT(2)=20
        ELSE IF(IB.EQ.4) THEN
          IDEFLT(1)=3
          IDEFLT(2)=6
        ENDIF
        WRITE(CHAR2,'(I3)') IDEFLT(1)
        WRITE(CHAR3,'(I3)') IDEFLT(2)
        FORMAT='('' Enter the number of Gauss points in the Xi('
     '    //CHAR1(1:1)//') direction for the'''
     '    //'/$,'' low and high order schemes '
     '    //'['//CHAR2(1:3)//','//CHAR3(1:3)//']: '',I3,I3)'
        IF(IOTYPE.EQ.3) THEN
          IDATA(1)=NGAP(ni,NBASEF(nb,NBASEF(nb,0)-1))
          IDATA(2)=NGAP(ni,NBASEF(nb,1))
        ENDIF
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,2,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,64,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          NGLIMITS(ni,nb,1)=IDATA(1) !Low order
          NGLIMITS(ni,nb,2)=IDATA(2) !High order
          IF(NGLIMITS(ni,nb,1).GT.NGLIMITS(ni,nb,2)) THEN
!           !Entered in the wrong order
            NGLIMITS(ni,nb,1)=IDATA(2)
            NGLIMITS(ni,nb,2)=IDATA(1)
          ENDIF
          NGAP(ni,nb)=NGLIMITS(ni,nb,2) !High order scheme

          IF(IB.LE.3) THEN
            IBT(1,ni,nb)=1
            IBT(2,ni,nb)=IB
            NNT(nb)=NNT(nb)*(IB+1)
          ELSE
            IBT(1,ni,nb)=2
            IBT(2,ni,nb)=1
            NNT(nb)=NNT(nb)*2
          ENDIF
          CALL ASSERT(NNT(nb).LE.NNM,'>>Need to increase NNM',
     '      ERROR,*9999)
        ENDIF
        IF(IB.GE.4) !quadratic/cubic Hermite
     '    NUMHERMXI=NUMHERMXI+1
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

      IF(IOTYPE.NE.3) THEN
c cpb 20/7/95 Temporary nn loop
        DO nn=1,NNT(nb)
          nk=0
          IF(NIT(nb).EQ.1) THEN
            DO ni1=1,IBT(1,1,nb)
              nk=nk+1
              IDO(nk,nn,1,nb)=ni1
            ENDDO
          ELSE IF(NIT(nb).EQ.2) THEN

C LKC 24-APR-1998 added assert
            CALL ASSERT(NKM.GE.IBT(1,2,nb)*IBT(1,1,nb),
     '        '>>Increase NKM',ERROR,*9999)

            DO ni2=1,IBT(1,2,nb)
              DO ni1=1,IBT(1,1,nb)
                nk=nk+1
                IDO(nk,nn,1,nb)=ni1
                IDO(nk,nn,2,nb)=ni2
              ENDDO
            ENDDO
          ENDIF
C KAT 14Dec99: cpb 11/9/95 Adding zero cross derivatives. This is not
C         the best way to do this but it is only temporary.
          IF(NBCD(nb).EQ.1) THEN
            IF(nk.EQ.4) THEN
              nk=3
              DO ni=1,NIT(nb)
                IDO(4,nn,ni,nb)=0
              ENDDO !ni
            ENDIF
          ENDIF
        ENDDO !nn
        NKT(0,nb)=nk
        NST(nb)=0
        DO nn=1,NNT(nb) !AJP 25-5-93  All nkt(nn,nb) the same for now
          NKT(nn,nb)=NKT(0,nb)
          NST(nb)=NST(nb)+NKT(nn,nb)
        ENDDO
      ENDIF

      IF(NST(nb).GT.0) THEN
        IF(iotype.ne.3) THEN
          nn=0
          IF(NIT(nb).EQ.1) THEN
            DO n1=1,IBT(2,1,nb)+1
              nn=nn+1
              INP(nn,1,nb)=n1
            ENDDO
          ELSE IF(NIT(nb).EQ.2) THEN
            DO n2=1,IBT(2,2,nb)+1
              DO n1=1,IBT(2,1,nb)+1
                nn=nn+1
                INP(nn,1,nb)=n1
                INP(nn,2,nb)=n2
              ENDDO
            ENDDO
          ENDIF
          CALL ASSERT(nn.LE.NNM,'>>Increase NNM',ERROR,*9999)
        ENDIF

        k=0
        DO nn=1,NNT(nb)
          DO ni=1,NIT(nb)
            k=k+1
            IDEFLT(k)=INP(nn,ni,nb)
            WRITE(CHAR1(k:k),'(I1)') IDEFLT(k)
          ENDDO
        ENDDO
        CALL STRING_TRIM(CHAR1,IBEG,IEND)
        FORMAT='($,'' Enter the node position indices ['//
     '    CHAR1(IBEG:IEND)//']: '',40I2)'
        IF(IOTYPE.EQ.3) THEN
          k=0
          DO nn=1,NNT(nb)
            DO ni=1,NIT(nb)
              k=k+1
              IDATA(k)=INP(nn,ni,nb)
            ENDDO
          ENDDO
        ENDIF
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,k,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,4,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(iotype.ne.3) THEN
          k=0
          DO nn=1,NNT(nb)
            DO ni=1,NIT(nb)
              k=k+1
              INP(nn,ni,nb)=IDATA(k)
            ENDDO
          ENDDO
        ENDIF

        IF(NKT(0,nb).GT.1) THEN
          k=0
          DO nk=1,NKT(0,nb)
            DO ni=1,NIT(nb)
              k=k+1
              IDEFLT(k)=IDO(nk,1,ni,nb)
              WRITE(CHAR1(k:k),'(I1)') IDEFLT(k)
            ENDDO
          ENDDO
          CALL STRING_TRIM(CHAR1,IBEG,IEND)
          FORMAT='($,'' Enter the derivative order indices ['//
     '      CHAR1(IBEG:IEND)//']: '',40I2)'
          IF(IOTYPE.EQ.3) THEN
            k=0
            DO nk=1,NKT(0,nb)
              DO ni=1,NIT(nb)
                k=k+1
                IDATA(k)=IDO(nk,1,ni,nb)
              ENDDO
            ENDDO
          ENDIF
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      k,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     &      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(iotype.ne.3) THEN
            k=0
            DO nk=1,NKT(0,nb)
              DO ni=1,NIT(nb)
                k=k+1
                IDO(nk,1,ni,nb)=IDATA(k)
              ENDDO !ni
            ENDDO !nk
            DO nn=2,NNT(nb)
              DO nk=1,NKT(0,nb)
                DO ni=1,NIT(nb)
                  k=k+1
                  IDO(nk,nn,ni,nb)=IDO(nk,1,ni,nb)
                ENDDO !ni
              ENDDO !nk
            ENDDO !nn
          ENDIF
        ENDIF
        DO nn=1,NNT(nb)
          DO nk=1,NKT(nn,nb)
            IDO1=IDO(nk,nn,1,nb)
            IDO2=1
            IDO3=1
            IF(NIT(nb).GE.2) IDO2=IDO(nk,nn,2,nb)
            IF(NIT(nb).EQ.3) IDO3=IDO(nk,nn,3,nb)
            IF(IDO1.EQ.1.AND.IDO2.EQ.1.AND.IDO3.EQ.1) IDO(nk,nn,0,nb)=1
            IF(IDO1.EQ.2.AND.IDO2.EQ.1.AND.IDO3.EQ.1) IDO(nk,nn,0,nb)=2
            IF(IDO1.EQ.1.AND.IDO2.EQ.2.AND.IDO3.EQ.1) IDO(nk,nn,0,nb)=4
            IF(IDO1.EQ.2.AND.IDO2.EQ.2.AND.IDO3.EQ.1) IDO(nk,nn,0,nb)=6
            IF(IDO1.EQ.1.AND.IDO2.EQ.1.AND.IDO3.EQ.2) IDO(nk,nn,0,nb)=7
            IF(IDO1.EQ.2.AND.IDO2.EQ.1.AND.IDO3.EQ.2) IDO(nk,nn,0,nb)=9
            IF(IDO1.EQ.1.AND.IDO2.EQ.2.AND.IDO3.EQ.2) IDO(nk,nn,0,nb)=10
            IF(IDO1.EQ.2.AND.IDO2.EQ.2.AND.IDO3.EQ.2) IDO(nk,nn,0,nb)=11
          ENDDO !nk
        ENDDO !nn
      ENDIF

      IF(iotype.ne.3)THEN !Don't recalculate PG if writing the file
        NBTOP=NBT !Highest basis function number so far
c PJH 9Sep95     IOD=.FALSE.
        CALL GAUSS5(IBT,IDO,INP,nb,NBTOP,NDET,NGAP,DET,PG,WG,XIG,
     '    ERROR,*9999)
      ENDIF

      CALL EXITS('BASIS5')
      RETURN
 9999 CALL ERRORS('BASIS5',ERROR)
      CALL EXITS('BASIS5')
      RETURN 1
      END


