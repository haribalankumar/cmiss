      SUBROUTINE BASIS4(IBT,IDO,INP,nb,NGAP,PG,WG,XIG,ERROR,*)

C#### Subroutine: BASIS4
C###  Description:
C###    BASIS4 is Fourier Series tensor product basis. Similar to
C###    BASIS1 except Fourier series is always in the last ni position
C###    and number of nodes is not doubled beyond that required by the
C###    Lagrange/Hermite part. INP is only defined for ni=1,..,NIT-1.
C**** IBT(1,ni,nb) = 9 for Fourier basis
C**** IBT(2,ni,nb) = # Fourier coefficients

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'four00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM),IDO(NKM,NNM,0:NIM),INP(NNM,NIM),nb,
     '  NGAP(NIM,NBM)
      REAL*8 PG(NSM,NUM,NGM),WG(NGM),XIG(NIM,NGM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IB,IBEG,ICHAR,IDO1,IDO2,IDO3,IEND,INFO,k,n1,n2,n3,
     '  ni,ni1,ni2,ni3,ni4,nk,nn,NOQUES
      CHARACTER CHAR1*100,CHAR2*1
      LOGICAL FILEIP

      CALL ENTERS('BASIS4',*9999)
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
c     CALL ASSERT(NIT(nb).GT.1,'>>NIT must be 2 or higher',ERROR,*9999)
      IF(IOTYPE.EQ.1.AND.NIT(nb).EQ.3) THEN
        WRITE(OP_STRING,'('' >>Remember to define 2D face basis '','
     '    //'''functions'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF
      IF(NIT(nb).GT.1) THEN
        WRITE(OP_STRING,
     '    '('' >>Fourier series must be defined in last xi position'')')
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

      IF(iotype.ne.3) THEN
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

      NOQUES=0
      DO ni=1,NIT(nb)
        WRITE(CHAR1,'(I1)') ni
        FORMAT='(/'' The interpolant in the Xi('//
     '    CHAR1(1:1)//') direction is [1]: '''//
     '    '/''   (1) Linear Lagrange'''//
     '    '/''   (2) Quadr. Lagrange'''//
     '    '/''   (3) Cubic  Lagrange'''//
     '    '/''   (4) Cubic  Hermite'''//
     '    '/''   (5) Fourier series'''//
     '    '/$,''    '',I1)'
        IF(IOTYPE.EQ.3) THEN
          IF(IBT(1,ni).EQ.1) THEN  !Lagrange
            IB=IBT(2,ni)
          ELSE IF(IBT(1,ni).EQ.2) THEN  !Hermite
            IB=4
          ELSE IF(IBT(1,ni).EQ.9) THEN  !Fourier
            IB=5
          ENDIF
          IDATA(1)=IB
        ENDIF
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,5,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(iotype.ne.3) IB=IDATA(1)

        IF(IB.EQ.5) THEN !Fourier
          FORMAT='($,'' Enter number of harmonics [1]: '',I2)'
          IDEFLT(1)=1
          IF(IOTYPE.EQ.3) IDATA(1)=(IBT(2,ni)-1)/2
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NKM,
     &      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(iotype.ne.3) THEN
            IBT(2,ni)=1+2*IDATA(1)
          ENDIF
          FORMAT='($,'' Enter ang. freq. of fundamental [1.0]: '','
     '      //'E12.4)'
          RDEFLT(1)=1.0D0
          IF(IOTYPE.EQ.3) RDATA(1)=OMEGA
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     &      LDATA,LDEFLT,RDATA,RDEFLT,0.0D0,RMAX,INFO,ERROR,*9999)
          IF(iotype.ne.3) OMEGA=RDATA(1)
        ENDIF

        IF(nb.GT.1) THEN
          IF(NIT(nb).EQ.NIT(nb-1)) THEN
            IDEFLT(1)=NGAP(ni,nb-1)
          ENDIF
        ELSE
          IF(IB.EQ.1) THEN
            IDEFLT(1)=2
          ELSE IF(IB.GT.1) THEN
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
        IF(iotype.ne.3) THEN
          NGAP(ni,nb)=IDATA(1)
          NGT(nb)=NGT(nb)*NGAP(ni,nb)
          CALL ASSERT(NGT(nb).LE.NGM,'>>Need to increase NGM',
     '      ERROR,*9999)
          IF(IB.LE.3) THEN !Lagrange
            IBT(1,ni)=1
            IBT(2,ni)=IB
            NNT(nb)=NNT(nb)*(IB+1)
          ELSE IF(IB.EQ.4) THEN !Hermite
            IBT(1,ni)=2
            IBT(2,ni)=1
            NNT(nb)=NNT(nb)*2
          ELSE IF(IB.EQ.5) THEN !Fourier
            IBT(1,ni)=9
          ENDIF
          CALL ASSERT(NNT(nb).LE.NNM,'>>Need to increase NNM',
     '      ERROR,*9999)
        ENDIF
      ENDDO !ni

      IF(iotype.ne.3) THEN
c cpb 14/7/95 Temporary nn loop
        DO nn=1,NNT(nb)
          nk=0
          IF(NIT(nb).EQ.2) THEN
            DO ni2=1,IBT(2,2) !#terms in Fourier series
              DO ni1=1,IBT(1,1) !1 for Lagrange, 2 for Hermite
                nk=nk+1
                IDO(nk,nn,1)=ni1
                IDO(nk,nn,2)=ni2
              ENDDO
            ENDDO
          ELSE IF(NIT(nb).EQ.3) THEN
            DO ni3=1,IBT(2,3)
              DO ni2=1,IBT(1,2)
                DO ni1=1,IBT(1,1)
                  nk=nk+1
                  IDO(nk,nn,1)=ni1
                  IDO(nk,nn,2)=ni2
                  IDO(nk,nn,3)=ni3
                ENDDO
              ENDDO
            ENDDO
          ELSE IF(NIT(nb).EQ.4) THEN
            DO ni4=1,IBT(2,3)
              DO ni3=1,IBT(1,3)
                DO ni2=1,IBT(1,2)
                  DO ni1=1,IBT(1,1)
                    nk=nk+1
                    IDO(nk,nn,1)=ni1
                    IDO(nk,nn,2)=ni2
                    IDO(nk,nn,3)=ni3
                    IDO(nk,nn,4)=ni4
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
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
          IF(NIT(nb).EQ.2) THEN
            DO n1=1,IBT(2,1)+1
              nn=nn+1
              INP(nn,1)=n1
            ENDDO
          ELSE IF(NIT(nb).EQ.3) THEN
            DO n2=1,IBT(2,2)+1
              DO n1=1,IBT(2,1)+1
                nn=nn+1
                INP(nn,1)=n1
                INP(nn,2)=n2
              ENDDO
            ENDDO
          ELSE IF(NIT(nb).EQ.4) THEN
            DO n3=1,IBT(2,3)+1
              DO n2=1,IBT(2,2)+1
                DO n1=1,IBT(2,1)+1
                  nn=nn+1
                  INP(nn,1)=n1
                  INP(nn,2)=n2
                  INP(nn,3)=n3
                ENDDO
              ENDDO
            ENDDO
          ENDIF
          CALL ASSERT(nn.LE.NNM,'>>Need to increase NNM',ERROR,*9999)
        ENDIF

        K=0
        DO nn=1,NNT(nb)
          DO ni=1,NIT(nb)
            K=K+1
            IDEFLT(K)=INP(nn,ni)
          ENDDO
        ENDDO
        CALL STRING_TRIM(CHAR1,IBEG,IEND)
        FORMAT='($,'' Enter the node position indices [default]: '','
     '    //'40I2)'
        IF(IOTYPE.EQ.3) THEN
          k=0
          DO nn=1,NNT(nb)
            DO ni=1,NIT(nb)
              k=k+1
              IDATA(k)=INP(nn,ni)
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
              INP(nn,ni)=IDATA(k)
            ENDDO
          ENDDO
        ENDIF

        DO nn=1,NNT(nb)
          DO nk=1,NKT(nn,nb)
            IDO1=IDO(nk,nn,1)
            IDO2=1
            IDO3=1
            IF(NIT(nb).GE.3) IDO2=IDO(nk,nn,2)
            IF(NIT(nb).EQ.4) IDO3=IDO(nk,nn,3)
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

      ENDIF

      CALL GAUSS4(IBT,IDO,INP,nb,NGAP(1,nb),PG,WG,XIG,ERROR,*9999)

      CALL EXITS('BASIS4')
      RETURN
 9999 CALL ERRORS('BASIS4',ERROR)
      CALL EXITS('BASIS4')
      RETURN 1
      END


