      SUBROUTINE OPBASE1(IBT,IDO,INP,NAN,nb,NFAM,NGAP,NKEF,NNF,NNL,nu,
     '  PG,XIG,ERROR,*)

C#### Subroutine: OPBASE1
C###  Description:
C###    OPBASE1 outputs basis function nb.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'four00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'jtyp00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),nb,NFAM,NGAP(NIM,NBM),NKEF(0:4,16,6,NBFM),
     '  NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),nu
      REAL*8 PG(NSM,NUM,NGM,NBM),XIG(NIM,NGM,NBM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IB,na,nbb,nf,ni,ng,nk,nl,nn,ns
      LOGICAL DIFFIDO
      CHARACTER CHAR*6,FORMAT2*100,
     '  NAME3(5)*7,NAME4(9)*12,NAME6(7)*26,NAME7(0:7)*34,FORMAT1*100
C     '  NAME1(3)*6,NAME2(4)*10,
C     '  NAME5(2)*8

C      DATA NAME1(1)/'Length'/,NAME1(2)/'  Area'/,NAME1(3)/'Volume'/,
C     '     NAME2(1)/'l.Lagrange'/,NAME2(2)/'q.Lagrange'/,
C     '     NAME2(3)/'c.Lagrange'/,NAME2(4)/'c.Hermite '/
      DATA NAME3 /'linear ','quadr. ','cubic  ','quadr 1','quadr 2'/
      DATA NAME4 /'Lagrange    ','Hermite     ',
     '            'Simplex     ','Serendipity ',
     '            'Sector      ','Sector      ',
     '            'Transition  ','Singular    ','Fourier    '/
C      DATA NAME5 /'External','Internal'/
C cpb 5/12/96 re-adding arc-length scaling
c cpb swapping over nbi = 4/5
C      DATA NAME6 /'unit scale factors      ',
C     '            'element factors read in ',
C     '            'global factors read in  ',
C     '            'calc.d from arc length  ',
C     '            'calc.d from angle change'/
C KAT 28Aug98: Adding harmonic mean
C      DATA NAME6 /'unit scale factors         ',
C     '            'element factors read in    ',
C     '            'global factors read in     ',
C     '            'calc.d from angle change   ',
C     '            'calc.d from ave. arc length',
C     '            'calc.d from arc length     '/
C                  123456789012345678901234567890
      DATA NAME6 /'unit                      ',
     '            'element based (read in)   ',
     '            'node based (read in)      ',
     '            'calc.d from angle change  ',
     '            'calc.d from arc length    ',
     '            'arithmetic mean arc length',
     '            'harmonic mean arc length  '/
      DATA NAME7 /'Auxiliary basis only             ',
     '            'Lagrange/Hermite tensor product  ',
     '            'Simplex/Serendipity/Sector       ',
     '            'B-spline tensor product          ',
     '            'Fourier Series basis             ',
     '            'Boundary element Lagrange/Hermite',
     '            'Boundary element Simplex/Sector  ',
     '            'Extended Lagrange                '/

      CALL ENTERS('OPBASE1',*9999)
      nbb=nb !the global basis function which was passed to opbase1
      nb=NFBASE(1,nbb) !family name of basis passed to opbase1
      IF(JTYP8.NE.2) THEN
        FORMAT='(/''  Basis function type '',I2,'': '',A,'
     '     //'/3X,'' The no. of Xi-coordinates  =  '',I1,'
     '     //'/3X,'' The no. of element nodes   = '',I2)'
        WRITE(OP_STRING,FORMAT) nb,NAME7(NBC(nb)),NIT(nb),NNT(nb)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ELSE IF(JTYP8.EQ.2) THEN !list basis family
        IF(NFAM.EQ.0) THEN
          FORMAT='(/''  Basis function type '',I2,'': '',A,'
     '       //'/3X,'' The no. of Xi-coordinates  =  '',I1,'
     '       //'/3X,'' The no. of element nodes   = '',I2)'
          WRITE(OP_STRING,FORMAT) nb,NAME7(NBC(nb)),NIT(nb),NNT(nb)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE
          FORMAT='(/''  Basis function type '',I2,'': '',A)'
          WRITE(OP_STRING,FORMAT) nb,NAME7(NBC(nb))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          FORMAT='(/''  Family basis number '','
     '    //'I2,'' Child basis number '
     '    //''',I2,'' Global basis number '',I2)'
          WRITE(OP_STRING,FORMAT) nb,NFBASE(2,nbb),nbb
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          FORMAT='(/3X,'' The no. of Xi-coordinates  =  '',I1,'
     '       //'/3X,'' The no. of element nodes   = '',I2)'
          WRITE(OP_STRING,FORMAT) NIT(nb),NNT(nb)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
      ENDIF
      IF(NKT(0,nb).GT.0) THEN
        FORMAT='(3X,'' The Xi('',I1,'') basis function   =  '',A,1X,A)'
        DO ni=1,NIT(nb)
          IF(IBT(1,ni,nb).EQ.2) THEN
            IF(IBT(2,ni,nb).EQ.1) THEN
              IB=3
            ELSE IF(IBT(2,ni,nb).EQ.2) THEN
              IB=4
            ELSE IF(IBT(2,ni,nb).EQ.3) THEN
              IB=5
            ENDIF
          ELSE
            IB=IBT(2,ni,nb)
          ENDIF
          IF(IBT(1,ni,nb).ne.9) THEN
            IF(IBT(1,ni,nb).EQ.5.OR.IBT(1,ni,nb).EQ.6) THEN
              IF(IB.EQ.4) THEN
                WRITE(OP_STRING,'(3X,'' The Xi('',I1,'') basis '
     '            //'function   =  '',A,1X,A,'' (Xi('',I1,'')='',I1,'
     '            //''')'')') ni,'Hermite',NAME4(IBT(1,ni,nb)),
     '            IBT(3,ni,nb),(IBT(1,ni,nb)-5)
              ELSE
                WRITE(OP_STRING,'(3X,'' The Xi('',I1,'') basis '
     '            //'function   =  '',A,1X,A,'' (Xi('',I1,'')='',I1,'
     '            //''')'')') ni,NAME3(IB),NAME4(IBT(1,ni,nb)),
     '            IBT(3,ni,nb),(IBT(1,ni,nb)-5)
              ENDIF
            ELSE
              WRITE(OP_STRING,FORMAT) ni,NAME3(IB),NAME4(IBT(1,ni,nb))
            ENDIF
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE
            WRITE(CHAR,'(I6)') IB
            WRITE(OP_STRING,FORMAT) ni,CHAR,
     '        NAME4(IBT(1,ni,nb))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF
      IF(NIT(nb).EQ.1) THEN
        FORMAT='(3X,'' The no. of Gauss  points   = '',I2)'
      ELSE IF(NIT(nb).EQ.2) THEN
        FORMAT='(3X,'' The no. of Gauss  points   = '',I2,'' *'',I2)'
      ELSE IF(NIT(nb).EQ.3) THEN
        FORMAT='(3X,'' The no. of Gauss  points   = '',2(I2,'' *''),I2)'
      ELSE IF(NIT(nb).EQ.4) THEN
        FORMAT='(3X,'' The no. of Gauss  points   = '',3(I2,'' *''),I2)'
      ENDIF
      WRITE(OP_STRING,FORMAT) (NGAP(ni,nbb),ni=1,NIT(nb))
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      IF(NKT(0,nb).GT.0.AND.NNT(nb).GT.0) THEN
        FORMAT='(3X,'' The no. of derivs/variable = '',I2,'
     '      //'/,3X,'' Indices for nodal position = '',1X,50(I1))'
        WRITE(OP_STRING,FORMAT) NKT(0,nb),
     '    ((INP(nn,ni,nb),ni=1,NIT(nb)),nn=1,NNT(nb))
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DIFFIDO=.FALSE.
        DO nn=1,NNT(nb)
          IF(NKT(nn,nb).NE.NKT(0,nb)) DIFFIDO=.TRUE.
        ENDDO !ni
        IF(DIFFIDO) THEN
          DO nn=1,NNT(nb)
            FORMAT='(3X,'' Ind. for der. ord. (nn='',I2,'') = '',1X,'
     '        //'50(I1))'
            WRITE(OP_STRING,FORMAT) nn,((IDO(nk,nn,ni,nb),
     '        ni=1,NIT(nb)),nk=1,NKT(nn,nb))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO !nn
        ELSE
          FORMAT='(3X,'' Indices for deriv. order   = '',1X,50(I1))'
          WRITE(OP_STRING,FORMAT) ((IDO(nk,1,ni,nb),ni=1,NIT(nb)),
     '      nk=1,NKT(0,nb))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
      ELSE IF(NKT(0,nb).GT.0.AND.NNT(nb).EQ.0) THEN
        FORMAT='(3X,'' The no. of spline functions= '',I2,'
     '       //'/3X,'' The no. of polynomial terms= '',I2)'
        WRITE(OP_STRING,FORMAT) NST(nb),NKT(0,nb)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF

      IF(NBC(nb).EQ.7) THEN !Extended Lagrange (mg collocation)
        FORMAT='(3X,'' # multigrid colloc. levels = '',I2)'
        WRITE(OP_STRING,FORMAT) NMGT
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      ELSE IF(NBC(nb).ne.4) THEN
        FORMAT='(3X,'' # auxiliary parameters     = '',I2)'
        WRITE(OP_STRING,FORMAT) NAT(nb)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(NABTYP(nb).EQ.1) THEN      !Legendre aux basis
          FORMAT='(3X,'' Legendre auxiliary basis'')'
          WRITE(OP_STRING,FORMAT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          FORMAT='(3X,'' Degree  of aux. param. ('',I1,'') =  '','
     '      //'2(I2,''*''),I2)'
        ELSE IF(NABTYP(nb).EQ.2) THEN !Fourier  aux basis
          FORMAT='(3X,'' Fourier  auxiliary basis'')'
          WRITE(OP_STRING,FORMAT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          FORMAT='(3X,'' 1/2 wave# of aux param ('',I1,'') =  '','
     '      //'2(I2,''*''),I2)'
        ELSE IF(NABTYP(nb).EQ.3) THEN !Pressure aux basis
          FORMAT='(3X,'' Pressure auxiliary basis'')'
          WRITE(OP_STRING,FORMAT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          FORMAT='(3X,'' Degree  of aux. param. ('',I1,'') =  '','
     '      //'2(I2,''*''),I2)'
        ENDIF
        IF(NABTYP(nb).GT.0) THEN
          DO na=1,NAT(nb)
            WRITE(OP_STRING,FORMAT) na,(NAN(ni,na,nb),ni=1,3)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDIF
      ENDIF

      IF(NBC(nb).LE.2.OR.NBC(nb).EQ.5.OR.NBC(nb).EQ.6) THEN !Lagrange/Hermite tensor prod basis or simplex
        WRITE(OP_STRING,'(3X,'' Line scaling factor type   =  '',A)')
     '    NAME6(NBI(nb))
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ELSE IF(NBC(nb).EQ.4) THEN !Fourier Series basis
        WRITE(OP_STRING,'(3X,'' Number of harmonics        = '',I2)')
     '    (IBT(2,NIT(nb),nb)-1)/2
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(3X,'' Number of coefficients     = '',I2)')
     '    IBT(2,NIT(nb),nb)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(3X,'' Ang. freq. of fundamental  = '','
     '    //'E11.4)') OMEGA
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF

      IF(JTYP8.EQ.1) THEN

        WRITE(OP_STRING,'(/'' Number of lines in the basis: '',I2)')
     '    NLE(nb)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO nl=1,NLE(nb)
          WRITE(OP_STRING,'('' Node #s for line '',I2,'' : '',4I3)')
     '      nl,(NNL(nn,nl,nb),nn=1,4)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO !nl
        IF(NIT(nb).GT.2) THEN
          WRITE(OP_STRING,'(/'' Number of faces in the basis: '',I1)')
     '      NFE(nb)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO nf=1,NFE(nb)
            WRITE(OP_STRING,'('' Xi direction normal to face : '',I1)')
     '        NNF(1,nf,nb)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' Node #s for face '',I1,'' : '',16I3)')
     '        nf,(NNF(nn,nf,nb),nn=1,NNF(0,nf,nb))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            DO nn=1,NNF(0,nf,nb)
              WRITE(OP_STRING,'(''   Deriv #s for face nn='',I2,'
     '          //''' : '',10I2)') nn,(NKEF(nk,nn,nf,nb),
     '          nk=1,NKEF(0,nn,nf,nb))
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDDO !nn
          ENDDO !nf
        ENDIF

        DO ng=1,NGT(nbb)
          WRITE(OP_STRING,'(''    Xi positions at Gauss point '',I3,'
     '      //''': '',3D12.4)') ng,(XIG(ni,ng,nbb),ni=1,NIT(nb))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
        IF(NNT(nb).GT.0) THEN
          FORMAT1='(/3X,'' PG(ns,nu,ng,nb='',I2,'
     '      //'''),ng=1,2,3,..... :''/)'
          FORMAT2='(3X,'' ns='',I2,'' nu='',I2,6D11.3,/(15X,6D11.3))'
        ELSE IF(NAT(nb).GT.0) THEN
          FORMAT1='(/3X,'' PG(na,nu,ng,nb='',I2,'
     '      //'''),ng=1,2,3,..... :''/)'
          FORMAT2='(3X,'' na='',I2,'' nu='',I2,6D11.3,/(15X,6D11.3))'
        ENDIF
        WRITE(OP_STRING,FORMAT1) nb
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C MPN 7/4/93 - Separating out NST(nb) and NAT(nb)
c        DO ns=1,NST(nb)
        DO ns=1,NST(nb)+NAT(nb)
          WRITE(OP_STRING,FORMAT2) ns,nu,
     '      (PG(ns,nu,ng,nbb),ng=1,NGT(nbb))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

      CALL EXITS('OPBASE1')
      RETURN
 9999 CALL ERRORS('OPBASE1',ERROR)
      CALL EXITS('OPBASE1')
      RETURN 1
      END


