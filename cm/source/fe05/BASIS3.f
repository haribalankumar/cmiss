      SUBROUTINE BASIS3(IBT,INP,nb,NGAP,ERROR,*)

C#### Subroutine: BASIS3
C###  Description:
C###    BASIS3 inputs parameters for B-spline tensor product basis
C###    functions and Gauss-Legendre or Gauss-Lobatto quadrature.
C**** The PG array for a Spline basis contains the Gauss point
C**** evaluations of the polynomial powers of Xi(1),Xi(2) and Xi(3).
C**** NST(nb) is the number of dofs coupled to each element.
C**** NKT(0,nb) "   "   "   " polynomial terms in the Spline basis
C****
C**** 23rd March 1989: this basis is now defined using a knot vector.
C**** This is defined by the user and stored in the common block BSPLINE:
C**** NTKN(ni) is the number of knots in the ni direction
C**** BKNOT(nknot,ni) is the xi(ni) position of the knot (non-decreasing)
C**** MBSPL(ni) is the degree of the spline in the ni direction.
C**** NST(nb) now holds the number of d.o.f.s per element
C**** NNT(np)  "   "     "   "     "     "     "     "
C**** NAT(nb) is now zero
C**** NKT(0,nb) is now one
C**** INP(nn,ni,nb) is the xi index for each control point (cf Lagrange basis)

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'bspln00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM),INP(NNM,NIM),nb,NGAP(NIM,NBM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,IBEG,ICHAR,IEND,INFO,j,n,n1,n2,n3,ni,nkn,NKN2,nn,
     '  NOQUES
      REAL*8 T(100)
      CHARACTER CHAR1*100,CHAR2*3,CHAR3*12
      LOGICAL FILEIP

      CALL ENTERS('BASIS3',*9999)
      FILEIP=.FALSE.
      NOQUES=0
      INFO=2
      ICHAR=999

      FORMAT='($,'' Enter the number of Xi-coordinates [1]: '',I1)'
      IF(IOTYPE.EQ.3) IDATA(1)=NIT(nb)
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,3,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(iotype.ne.3) NIT(nb)=IDATA(1)

      NUT(nb)=NIT(nb)*NIT(nb)+2
      CALL ASSERT(NUT(nb).LE.NUM,'>>Increase NUM',ERROR,*9999)
C     NAT(nb)=1
      NGT(nb)=1
      NKT(0,nb)=1
      NST(nb)=1
      NNT(nb)=1
      DO ni=1,NIT(nb)
        NGAP(ni,nb)=1
        IBT(1,ni)=1
C news MPN 9-Jul-96 Initialise IBT
        IBT(2,ni)=0
        IBT(3,ni)=0
      ENDDO !ni
      DO ni=1,NIT(nb)
        WRITE(CHAR1,'(I1)') ni

        FORMAT='($,'' Enter number of knots in the Xi('//CHAR1(1:1)//
     '    ') direction [4]:'',I3)'
        IF(IOTYPE.EQ.3) IDATA(1)=NTKN(ni)
        IDEFLT(1)=4
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,20,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(iotype.ne.3) NTKN(ni)=IDATA(1)
        RDEFLT(1)=0.D0
        DO nkn=1,NTKN(ni)
          WRITE(CHAR2,'(I3)')nkn
C          CHAR3=CFROMR(RDEFLT(1),'(E12.5)')
          WRITE(CHAR3,'(E12.5)') RDEFLT(1)
          CALL STRING_TRIM(CHAR2,IBEG,IEND)
          FORMAT='($,'' For knot '//CHAR2(IBEG:IEND)//
     '      ' enter xi position ['//CHAR3(1:12)//']:'',E12.5)'
          IF(IOTYPE.EQ.3) RDATA(1)=RDEFLT(1)
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     &      LDATA,LDEFLT,RDATA,RDEFLT,RDEFLT(1),RMAX,INFO,ERROR,*9999)
          IF(iotype.ne.3) BKNOT(nkn,ni)=RDATA(1)
          RDEFLT(1)=BKNOT(nkn,ni)
        ENDDO

        FORMAT='(/'' Specify whether interpolant in the Xi('
     '    //CHAR1(1:1)//') direction is [1]:'''//
     '    '/''   (1) Linear B-spline'''//
     '    '/''   (2) Quadr. B-Spline'''//
     '    '/''   (3) Cubic  B-Spline'''//
     '    '/$,''    '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=MBSPL(ni)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,3,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(iotype.ne.3) MBSPL(ni)=IDATA(1)

        IF(nb.GT.1) THEN
          IF(NIT(nb).EQ.NIT(nb-1)) THEN
            IDEFLT(1)=NGAP(ni,nb-1)
          ENDIF
        ELSE
          IF(MBSPL(ni).EQ.1) THEN
            IDEFLT(1)=2
          ELSE IF(MBSPL(ni).GT.1) THEN
            IDEFLT(1)=3
          ENDIF
        ENDIF
        WRITE(CHAR2,'(I1)') IDEFLT(1)
        FORMAT='($,'' Enter the number of Gauss points in the Xi('//
     '    CHAR1(1:1)//') direction ['//CHAR2(1:1)//']: '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=NGAP(ni,nb)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,9,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(iotype.ne.3) NGAP(ni,nb)=IDATA(1)

        NGT(nb)=NGT(nb)*NGAP(ni,nb)
        CALL ASSERT(NGT(nb).LE.NGM,'>>Need to increase NGM',
     '    ERROR,*9999)
        NST(nb)=NST(nb)*(MBSPL(ni)+1)
        NNT(nb)=NNT(nb)*(MBSPL(ni)+1)
        CALL ASSERT(NNT(nb).LE.NNM,'>>Need to increase NNM',
     '    ERROR,*9999)
        IBT(1,ni)=5
        IBT(2,ni)=MBSPL(ni)

C       calculate discrete b-splines over refined mesh for rendering
        NKN2=1
        T(1)=BKNOT(1,ni)
        DO nkn=1,NTKN(ni)-1
          IF(BKNOT(nkn,ni).LT.BKNOT(nkn+1,ni)) THEN
            DO n=1,6
              NKN2=NKN2+1
              T(NKN2)=BKNOT(nkn,ni)+N*(BKNOT(nkn+1,ni)
     '          -BKNOT(nkn,ni))/6.0d0
            ENDDO
          ELSE
            NKN2=NKN2+1
            T(NKN2)=BKNOT(nkn,ni)
          ENDIF
        ENDDO
        NTKN2(ni)=NKN2
        IF(DOP) THEN
          WRITE(OP_STRING,*)' T=',(T(J),J=1,NTKN2(ni))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        DO j=1,NTKN2(ni)-(MBSPL(ni)+1)
          CALL OSLO(j,NTKN(ni),BKNOT(1,ni),T,MBSPL(ni)+1,ALPHA(1,j,ni))
        ENDDO
        IF(DOP) THEN
          DO i=1,NTKN(ni)-(MBSPL(ni)+1)
            WRITE(OP_STRING,'('' ALPHA('',I2,'',..,'',I2,'')='',/,'
     '        //'20(5F11.4,/))')
     '        i,ni,(ALPHA(i,j,ni),j=1,NTKN2(ni)-(MBSPL(ni)+1))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDIF

      ENDDO !ni

C     need to redefine this later
C     CALL GAUSS3(IBT,nb,NGAP(1,nb),PG,WG,XIG,ERROR,*9999)

      IF(DOP) THEN
        WRITE(OP_STRING,'('' ntkn(..)='',3I4)') (NTKN(ni),ni=1,3)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' mbspl(..)='',3I4)') (MBSPL(ni),ni=1,3)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO ni=1,NIT(nb)
          WRITE(CHAR1,'(I1)') ni
          WRITE(OP_STRING,'('' bknot(..,'//CHAR1(1:1)//')='','
     '      //'(/,8E12.5))') (BKNOT(nkn,ni),nkn=1,NTKN(ni))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

C     calculate inp
      IF(NST(nb).GT.0) THEN
        nn=0
        IF(NIT(nb).EQ.1) THEN
          DO n1=1,IBT(2,1)+1
            nn=nn+1
            INP(nn,1)=n1
          ENDDO
        ELSE IF(NIT(nb).EQ.2) THEN
          DO n2=1,IBT(2,2)+1
            DO n1=1,IBT(2,1)+1
              nn=nn+1
              INP(nn,1)=n1
              INP(nn,2)=n2
            ENDDO
          ENDDO
        ELSE IF(NIT(nb).EQ.3) THEN
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

      CALL EXITS('BASIS3')
      RETURN
 9999 CALL ERRORS('BASIS3',ERROR)
      CALL EXITS('BASIS3')
      RETURN 1
      END


