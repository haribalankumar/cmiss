      SUBROUTINE IPMAT5(GRNGLIST,ICQS_SPATIAL,ILPIN,ILTIN,
     '  IRCQS_SPATIAL,NBJ,NEELEM,NELIST,NGLIST,NMBIN,NPLIST,NPNODE,NQET,
     '  NQNE,NQS,NQXI,nr,nx,CE,CELL_RCQS_VALUE,CGE,CIN,CP,RCQS_SPATIAL,
     '  XIG,ERROR,*)

C#### Subroutine: IPMAT5
C###  Description:
C###    IPMAT5 inputs material parameters and determines constitutive
C###    law parameters from them (a straight except fot the pole-zero
C###    law which uses fibre distribution functionals to relate the
C###    parameters).

C**** For a polynomial in the principal strain invariants I1,I2,I3, or
C**** extension ratios L1,L2,L3, IT(1,nr),IT(2,nr) & IT(3,nr) record the
C**** maximum powers of I1-3,I2-3 & I3-1 or L1-1,L2-1 & L3-1,respec.
C**** ILP(il,1,nr,nx),il=1,ILT(1,nr,nx),
C****   where ILT(1,nr,nx)=(IT(1,nr)+1)*(IT(2,nr)+1)*(IT(3,nr)+1),
C****   records whether constitutive law parameters are:
C****   (1) Constant spatially: value  in CE(il,ne)
C****   (2) Piecewise constant: values in CE(il,ne)
C****   (3) Piecewise linear:   values in CP(il,np) basis NMB(il,nr,nx)
C****   (4) Defined by Gauss pts- values in CGE(nm,ng,ne,nx)
! CS 17/5/99 implementing, but putting in CGE
! C****   (4) Defined by Gauss pts- values in YG(ng,ne) -Not implemented
! new CS Feb 24 2000
C****   (5) Defined by Grid pts- values in CGE(nm,ng,ne,nx) interpolated from grid point fields
C**** and if parameter il is time varying ILP(il,1,nr,nx) is negative.
C**** NOTE: CIN, ILPIN, ILTIN and NMBIN store the corresponding
C**** information for the entered parameters (from which the
C**** constitutive parameters are derived).

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ipma50.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'titl50.cmn'
!     Parameter List
      INTEGER GRNGLIST(0:NEGM),ICQS_SPATIAL(NQISVM,NQM),
     '  ILPIN(NMM),ILTIN,IRCQS_SPATIAL(0:NQRSVM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),NGLIST(0:NGM),
     '  NMBIN(NMM),NPLIST(0:NPM),NPNODE(0:NP_R_M,0:NRM),NQET(NQSCM),
     '  NQNE(NEQM,NQEM),NQS(NEQM),NQXI(0:NIM,NQSCM),nr,nx
      REAL*8 CE(NMM,NEM),CELL_RCQS_VALUE(NQRM,NQVM),
     '  CGE(NMM,NGM,NEM),CIN(NMM,0:NGM,NNEPM),
     '  CP(NMM,NPM),RCQS_SPATIAL(NQRSVM,NQM),XIG(NIM,NGM,NBM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,i1,i2,i3,IBEG1,IBEG2,IBEG3,IBEG4,IBEG5,IBEG6,
     '  IBEG_T1,ICHAR,IEND1,IEND2,IEND3,IEND4,IEND5,IEND6,IEND_T1,
     '  il,INFO,IT1,ITOT,K52T,nb,ne,ng,nm_axial,nm_extn,
     '  noelem,nonode,NOQUES,np,num_entered
      CHARACTER CHAR1*100,CHAR2*100,CHAR3*100,
     '  CHAR4*100,CHAR5*100,CHAR6*100,CHAR7*3,TITLE1*200
      LOGICAL FILEIP,SETCONSTIT

      ICHAR=0 !temporary (needs deleting later)

      CALL ENTERS('IPMAT5',*9999)
      CALL ASSERT(NMM.EQ.35,' >>Arrays ILP,NMB in b13.cmn '
     '  //'and RMAT5,CMAT5 in ipma50.cmn must have first '
     '  //'dimension=NMM (change also BLK50)',ERROR,*9999)

      FILEIP=.FALSE.
      NOQUES=0

      TITLE1='/''  (1) Constant spatially'''//
     '  '/''  (2) Piecewise constant (defined by elements)'''//
     '  '/''  (3) Piecewise linear   (defined by nodes) '''//
     '  '/''  (4) Defined by Gauss points'''//
     '  '/''  (5) Defined by Grid points'''
      CALL STRING_TRIM(TITLE1,IBEG_T1,IEND_T1)

      IF(ITYP2(nr,nx).EQ.2) THEN !finite elasticity

        FORMAT='('' Stresses in constitutive law are '//
     '    'referred to [2]: '''//
     '    '/''   (1) Reference (theta) coordinates'''//
     '    '/''   (2) Body (fibre/transverse) coordinates'''//
     '    '/''   (3) Body coordinates with active fibre stress'''//
     '    '/$,''    '',I1)'
        IDEFLT(1)=2
        IF(IOTYPE.EQ.3) IDATA(1)=KTYP53(nr)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) KTYP53(nr)=IDATA(1)
C VJ 9Dec2003 adding option to use precomputed gauss point stresses
        FORMAT='('' Specify whether the constitutive law is'//
     '    ' defined by [1]: '''//
     '    '/''   (1) a Green strain energy function (hyperelastic) '''//
     '    '/''   (2) a stress/strain-rate relation'''//
     '    '/''   (3) Gauss point stresses (grid coupling)'''//
     '    '/$,''    '',I1)'
C        FORMAT='('' Specify whether the constitutive law is'//
C     '    ' defined by [1]: '''//
C     '    '/''   (1) a Green strain energy function (hyperelastic) '''//
C     '    '/''   (2) a stress/strain-rate relation'''//
C     '    '/''  *(3) a quasi-static Creep law'''//
C     '    '/$,''    '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=KTYP54(nr)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,3,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) KTYP54(nr)=IDATA(1)

C ---------------- Define type of constitutive law -------------------

        IF(KTYP54(nr).EQ.1) THEN !hyperelastic
          FORMAT='('' Specify whether the strain energy W is given as a'
     '      //' function of [3]: '''//
     '      '/''   (1) the principal strain invariants '''//
     '      '/''   (2) the principal extension ratios '''//
     '      '/''   (3) the fibre & transverse strains '''//
     '      '/$,''    '',I1)'
          IDEFLT(1)=3
          IF(IOTYPE.EQ.3) IDATA(1)=KTYP55(nr)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) KTYP55(nr)=IDATA(1)

          IF(KTYP52(nr).EQ.1) THEN       !compressible
            K52T=1
          ELSE IF(KTYP52(nr).EQ.2) THEN  !incompressible
            K52T=2
          ELSE IF(KTYP52(nr).EQ.3) THEN  !incomp + fluid
            K52T=2
          ELSE IF(KTYP52(nr).EQ.4) THEN  !compr + fluid
            K52T=1
          ELSE IF(KTYP52(nr).EQ.5) THEN  !incomp+inext
            K52T=2
          ELSE IF(KTYP52(nr).EQ.6) THEN  !compr + fluid for lung
            K52T=1
            MEAN_VV_RATIO=1.d0
          ENDIF

          CHAR2=TITL56(1,K52T,KTYP55(nr))
          CHAR3=TITL56(2,K52T,KTYP55(nr))
          CHAR4=TITL56(3,K52T,KTYP55(nr))
          CHAR5=TITL56(4,K52T,KTYP55(nr))
          CHAR6=TITL56(5,K52T,KTYP55(nr))
          CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
          CALL STRING_TRIM(CHAR3,IBEG3,IEND3)
          CALL STRING_TRIM(CHAR4,IBEG4,IEND4)
          CALL STRING_TRIM(CHAR5,IBEG5,IEND5)
          CALL STRING_TRIM(CHAR6,IBEG6,IEND6)
          IF(KTYP55(nr).EQ.1) THEN !princ. strain invariants
            IDEFLT(1)=1
          ELSE IF(KTYP55(nr).EQ.2) THEN !extension ratios
            IDEFLT(1)=1
          ELSE IF(KTYP55(nr).EQ.3) THEN !fibre strains
            IDEFLT(1)=3
          ENDIF
          WRITE(CHAR1,'(I1)') IDEFLT(1)
          FORMAT='('' Specify the form of strain energy'//
     '      ' function W ['//CHAR1(1:1)//']: '''//
     '      '/''   (1) '//CHAR2(IBEG2:IEND2)//''''//
     '      '/''   (2) '//CHAR3(IBEG3:IEND3)//''''//
     '      '/''   (3) '//CHAR4(IBEG4:IEND4)//''''//
     '      '/''   (4) '//CHAR5(IBEG5:IEND5)//''''//
     '      '/''   (5) '//CHAR6(IBEG6:IEND6)//''''//
     '      '/$,''    '',I1)'
          IF(IOTYPE.EQ.3) IDATA(1)=KTYP56(nr)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,5,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) KTYP56(nr)=IDATA(1)

          IF((KTYP55(nr).EQ.3).AND.(KTYP56(nr).EQ.2)) THEN
            FORMAT='('' Specify the form of exponential '
     '        //'funtion[1]: '''//
     '        '/''   (1) W=C1*exp(Q) where'//
     '        ' Q=2*C2*(Ef+Es+En)+C3*(Ef^2)+C4*(Es^2+En^2+2*Esn^2)'//
     '        '+2*C5*(Efs^2+Enf^2)'''//
     '        '/''   (2) Dr J.W. Holmes distributed '
     '        //'fibre formulation'''//
     '        '/''   (3) Tong & Fung skin function'''//
     '        '/''   (4) W=C1*exp(Q) where'//
     '        ' Q=C2*Ef^2+C3*Es^2+C4*En^2+2*C5*Efs*Esf+2*C6*Efn*Enf'//
     '        '+2*C7*Esn*Ens'''//
     '        '/$,''    '',I1)'

            IDEFLT(1)=1
            IF(IOTYPE.EQ.3) IDATA(1)=KTYP55a(nr)
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,4,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) KTYP55a(nr)=IDATA(1)

            IF(KTYP55a(nr).EQ.3) THEN
              CALL ASSERT(
     '          KTYP51(nr).EQ.4.OR.KTYP51(nr).EQ.1.OR.KTYP51(nr).EQ.2,
     '          ' >>Tong and Fung only used for 2D problems',
     '          ERROR,*9999)
            ENDIF
          ENDIF
        ELSE IF(KTYP54(nr).EQ.2) THEN !stress/strain-rate relation
          KTYP55(nr)=3 !fibre & transverse strains
          KTYP56(nr)=6 !linear viscous relation

C news VJ 16Dec2003         
        ELSE IF(KTYP54(nr).EQ.3) THEN !Gauss point stresses (grid coupling)
          KTYP55(nr) = 3 !For CellML always use dW/dEij so must
C                        define wrt fibre and transverse coordinates
          KTYP56(nr)=0 !not needed for grid coupling 
C         no constitutive parameters are required here as all defined using CellML
          IT(1,nr)=-1
          IT(2,nr)=-1
          IT(3,nr)=-1
        ENDIF !KTYP54

C news VJ 6Jan2004 This if statement is obsolete - KTYP56 cannot be set to 0
C        IF(KTYP56(nr).EQ.0) THEN !linear elastic relation
C          IT(1,nr)=2
C          IT(2,nr)=0
C          IT(3,nr)=0
C newe VJ 6J
        IF(KTYP56(nr).EQ.1) THEN !polynomial function
C KAT 2001-08-09: Certainly don't want to reset IT if writing.
C         Can't see why it needs initialization any other time.
C          DO i=1,3
C            IT(i,nr)=0
C          ENDDO
          IF(KTYP55(nr).LE.2) THEN !princ strain invars or princ exten.s
            IT1=1
            IF(KTYP52(nr).EQ.1.OR.KTYP52(nr).EQ.4.OR.KTYP52(nr).EQ.6)
     '        THEN !compressible
              ITOT=3
            ELSE                                        !incompressible
              ITOT=2
            ENDIF
            DO i=IT1,ITOT
              WRITE(CHAR1,'(I1)') i
              CHAR2=TITL57(KTYP55(nr))
              FORMAT='($,'' Enter the maximum power of '//CHAR2(1:1)//
     '          CHAR1(1:1)//' [1]: '',I2)'
              IF(IOTYPE.EQ.3) IDATA(1)=IT(i,nr)
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,9,
     '          LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) IT(i,nr)=IDATA(1)
            ENDDO
          ELSE IF(KTYP55(nr).EQ.3) THEN !fibre strains
C           For W a fn of strains, 2 is the lowest meaningful power
            IF(KTYP52(nr).EQ.1) THEN
              IT1=1
            ELSE IF(KTYP52(nr).GE.2) THEN
              IT1=2
            ENDIF
            ITOT=3
            DO i=IT1,ITOT
              WRITE(CHAR1,'(I1)') i
              CHAR2=TITL57(KTYP55(nr))
              FORMAT='($,'' Enter the maximum power of '//CHAR2(1:1)//
     '          CHAR1(1:1)//' [2]: '',I2)'
              IF(IOTYPE.EQ.3) IDATA(1)=IT(i,nr)
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,9,
     '          LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
              IT(i,nr)=IDATA(1)
            ENDDO
          ENDIF

        ELSE IF(KTYP56(nr).GE.2.AND.KTYP56(nr).LE.4) THEN
          IF(KTYP56(nr).EQ.3)THEN
            IF(KTYP52(nr).NE.6) THEN
             IDEFLT(1)=1
             WRITE(CHAR1,'(I1)') IDEFLT(1)
             FORMAT='('' The shear terms of the pole-zero are [1]: '''//
     '         '/''   (1) Determined using the fibre '
     '         //'distribution model'''//
     '         '/''   (2) Input'''//
     '         '/$,''    '',I1)'
             IF(IOTYPE.EQ.3) IDATA(1)=KTYP5F(nr)
             CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '         FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,
     &         IDEFLT,1,2,
     '         LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
             IF(IOTYPE.NE.3) KTYP5F(nr)=IDATA(1)
             IF(KTYP5F(nr).EQ.2) THEN
               IDEFLT(1)=30
             ELSE
               IDEFLT(1)=IMAT5(KTYP56(nr),KTYP55(nr),KTYP52(nr))
             ENDIF
           ENDIF
          ELSE
C     exponential, pole-zero or fibre in fluid functions
             IDEFLT(1)=IMAT5(KTYP56(nr),KTYP55(nr),KTYP52(nr))
          ENDIF

! Alex 02-Dec-02: Added the below string for Blatz-Ko function

          IF(KTYP56(nr).EQ.2.AND.KTYP55(nr).EQ.1.AND.
     &      KTYP52(nr).EQ.1)THEN
           ! Blatz-Ko material
            IT(1,nr)=0
            IT(2,nr)=0
            IT(3,nr)=2
          ELSE IF(KTYP56(nr).EQ.3.AND.KTYP55(nr).EQ.1.AND.
     &      KTYP52(nr).EQ.6)THEN
           ! Fung exponential law
            IT(1,nr)=0
            IT(2,nr)=0
            IT(3,nr)=2
          ELSE
            WRITE(CHAR1,'(I2)') IDEFLT(1)
            CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
            FORMAT='($,'' Enter the number of constitutive law'//
     '      ' parameters ['//CHAR1(IBEG1:IEND1)//']: '',I2)'
            IF(IOTYPE.EQ.3) THEN
              IF(KTYP52(nr).EQ.1) THEN
                IT(1,nr)=0
              ELSE IF(KTYP52(nr).GE.2) THEN
                IT(1,nr)=IT(1,nr)+1
              ENDIF
              IDATA(1)=IT(1,nr)
            ENDIF
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     '        NMM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) IT(1,nr)=IDATA(1)
            IF(KTYP52(nr).EQ.1) THEN !compressible
              IT(1,nr)=0
            ELSE IF(KTYP52(nr).GE.2) THEN !incompressible
              IT(1,nr)=IT(1,nr)-1
            ENDIF
            IT(2,nr)=0
            IT(3,nr)=0
          ENDIF

        ELSE IF(KTYP56(nr).EQ.5) THEN !user defined fn
          IT(2,nr)=0
          IT(3,nr)=0
          FORMAT='($,'' Enter the number of user-defined material'//
     '      ' parameters [1]: '',I2)'
          IF(IOTYPE.EQ.3) IDATA(1)=IT(1,nr)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '      FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,9,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            IT(1,nr)=IDATA(1)
            IT(1,nr)=IT(1,nr)-1
          ENDIF

        ELSE IF(KTYP56(nr).EQ.6) THEN !stress/strain-rate relation
          IT(1,nr)=3
          IT(2,nr)=0
          IT(3,nr)=0

        ENDIF !ktyp56

C --------------- Define constitutive law parameters ------------------

        ILT(1,nr,nx)=(IT(1,nr)+1)*(IT(2,nr)+1)*(IT(3,nr)+1)
        CALL ASSERT(ILT(1,nr,nx).LE.NMM,' >>Increase NMM',
     '    ERROR,*9999)

        il=0
        DO i3=0,IT(3,nr)
          DO i2=0,IT(2,nr)
            DO i1=0,IT(1,nr)
              il=il+1

              IF(KTYP56(nr).EQ.0) THEN !constant
C               Set FORMAT for parameter variation question

              ELSE IF(KTYP56(nr).EQ.1) THEN !poly(fib/trans norm strains)
                WRITE(CHAR1,'(I1)') i1
                WRITE(CHAR2,'(I1)') i2
                WRITE(CHAR3,'(I1)') i3
!Compressable
                IF((KTYP52(nr).EQ.1.OR.KTYP52(nr).EQ.4.OR.KTYP52(nr)
     '            .EQ.6).AND.KTYP55(nr).EQ.1) THEN !invariants
                  FORMAT=
     '              '(/'' Specify whether the coefficient of (I1-3)^'//
     '              CHAR1(1:1)//'*(I2-3)^'//CHAR2(1:1)//'*(I3-1)^'//
     '              CHAR3(1:1)//' is [1]:'''//
     '              TITLE1(IBEG_T1:IEND_T1)//'/$,''   '',I1)'
                ELSE IF(KTYP52(nr).GE.2.AND.KTYP55(nr).EQ.1) THEN
                  FORMAT=
     '              '(/'' Specify whether the coefficient of (I1-3)^'//
     '              CHAR1(1:1)//'*(I2-3)^'//CHAR2(1:1)//' is [1]:'''//
     '              TITLE1(IBEG_T1:IEND_T1)//'/$,''   '',I1)'
                ELSE IF(KTYP52(nr).EQ.1.AND.KTYP55(nr).EQ.2) THEN
                  FORMAT=
     '              '(/'' Specify whether the coefficient of (L1-1)^'//
     '              CHAR1(1:1)//'*(L2-1)^'//
     '              CHAR2(1:1)//'*(L3-1)^'//CHAR3(1:1)//' is [1]:'''//
     '              TITLE1(IBEG_T1:IEND_T1)//'/$,''   '',I1)'
                ELSE IF(KTYP52(nr).GE.2.AND.KTYP55(nr).EQ.2) THEN
                  FORMAT=
     '              '(/'' Specify whether the coefficient of (L1-1)^'//
     '              CHAR1(1:1)//'*(L2-1)^'//CHAR2(1:1)//' is [1]:'''//
     '              TITLE1(IBEG_T1:IEND_T1)//'/$,''   '',I1)'
                ELSE IF(KTYP52(nr).EQ.1.AND.KTYP55(nr).EQ.3) THEN
                  FORMAT=
     '              '(/'' Specify whether the coefficient of E1^'//
     '              CHAR1(1:1)//'*E2^'//
     '              CHAR2(1:1)//'*E3^'//CHAR3(1:1)//' is [1]:'''//
     '              TITLE1(IBEG_T1:IEND_T1)//'/$,''   '',I1)'
                ELSE IF(KTYP52(nr).GE.2.AND.KTYP55(nr).EQ.3) THEN
                  FORMAT=
     '              '(/'' Specify whether the coefficient of E2^'//
     '              CHAR1(1:1)//'*E3^'//
     '              CHAR2(1:1)//' is [1]:'''//
     '              TITLE1(IBEG_T1:IEND_T1)//'/$,''   '',I1)'
                ENDIF

              ELSE IF(KTYP56(nr).EQ.2) THEN !exp(fib/trans strains)
                IF(KTYP52(nr).EQ.1.AND.KTYP55(nr).EQ.1) THEN
C                 compress.& p.strain invars

! Alex 20-N0v-02: set three parameters for Blatz-Ko function

                  CHAR1=CMAT5(il,ktyp56(nr),ktyp55(nr),ktyp52(nr))
                  CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
                  FORMAT='(/'' Specify whether the Blatz-Ko'//
     '              ' material parameter '//CHAR1(IBEG1:IEND1)//
     '              ' is [1]:'''//TITLE1(IBEG_T1:IEND_T1)//
     '              '/$,''   '',I1)'
                 ELSE IF(KTYP52(nr).GE.2.AND.KTYP55(nr).EQ.1) THEN
C                 incompressible
                  FORMAT='(/'' Specify whether the'//
     '              ' exponential material parameter is [1]:'''//
     '              TITLE1(IBEG_T1:IEND_T1)//'/$,''   '',I1)'
                ELSE IF(KTYP52(nr).GE.2.AND.KTYP55(nr).EQ.2) THEN
C                 incompressible
                  WRITE(CHAR1,'(I1)') I1+1
                  FORMAT='(/'' Specify whether the'//
     '              ' Ogden material parameter '//CHAR1(1:1)//
     '              ' is [1]:'''//TITLE1(IBEG_T1:IEND_T1)//
     '              '/$,''   '',I1)'
                ELSE IF(KTYP55(nr).EQ.3) THEN !fibre strains
                  IEND1=0
                  CALL APPENDI(IEND1,IL,CHAR1)
                  FORMAT='(/'' Specify whether the exponential'
     '              //'  material parameter '//CHAR1(1:IEND1)
     '              //' is [1]:'''//
     '              TITLE1(IBEG_T1:IEND_T1)//'/$,''   '',I1)'
                ENDIF

              ELSE IF(KTYP56(nr).EQ.3) THEN
C               pole zero in fibre & transv. strains
                IF(KTYP52(nr).EQ.1.AND.KTYP55(nr).EQ.1) THEN
C                 compress.& p. strain invars
                ELSE IF(KTYP52(nr).GE.2.AND.KTYP55(nr).EQ.1) THEN
                   IF(KTYP52(nr).EQ.6)THEN 
C                  comp+fluid for lung
                     FORMAT='(/'' Specify whether the'//
     '                ' exponential material parameter is [1]:'''//
     &               TITLE1(IBEG_T1:IEND_T1)//'/$,''   '',I1)'
                   ELSE
C                  incompressible
                   ENDIF
                ELSE IF(KTYP52(nr).GE.2.AND.KTYP55(nr).EQ.2) THEN
C                 incompressible
                ELSE IF(KTYP55(nr).EQ.3) THEN !fibre strains
                  IF(KTYP52(nr).EQ.1) K52T=1
                  IF(KTYP52(nr).GE.2) K52T=2
                  IF(KTYP5F(nr).EQ.1) THEN
                    CHAR1=CMAT5(IL,KTYP56(nr),KTYP55(nr),K52T)
                  ELSE
                    CHAR1=CMAT6(IL)
                  ENDIF
                  CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
                  FORMAT='(/'' Specify whether the '
     '              //CHAR1(IBEG1:IEND1)
     '              //' is [1]:'''//TITLE1(IBEG_T1:IEND_T1)
     '              //'/$,''   '',I1)'
                ENDIF

              ELSE IF(KTYP56(nr).EQ.4) THEN !fibre in fluid
                IF(KTYP52(nr).EQ.1.AND.KTYP55(nr).EQ.1) THEN
C                 compress.& p. strain invars
                ELSE IF(KTYP52(nr).GE.2.AND.KTYP55(nr).EQ.1) THEN
C                 incompressible
                ELSE IF(KTYP52(nr).GE.2.AND.KTYP55(nr).EQ.2) THEN
                ELSE IF(KTYP55(nr).EQ.3) THEN !fibre strains
C                 incompressible
                  IF(KTYP52(nr).EQ.1) K52T=1
                  IF(KTYP52(nr).GE.2) K52T=2
                  CHAR1=CMAT5(IL,KTYP56(nr),KTYP55(nr),K52T)
                  CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
                  FORMAT='(/'' Specify whether the '
     '              //CHAR1(IBEG1:IEND1)
     '              //' is [1]:'''//TITLE1(IBEG_T1:IEND_T1)
     '              //'/$,''   '',I1)'
                  ENDIF

              ELSE IF(KTYP56(nr).EQ.5) THEN !user defined in USER53
                WRITE(CHAR1,'(I1)') IL
                FORMAT='(/'' Specify whether the user-defined'
     '            //' material parameter '//CHAR1(1:1)//' is [1]:'''//
     '            TITLE1(IBEG_T1:IEND_T1)//'/$,''   '',I1)'

              ELSE IF(KTYP56(nr).EQ.6) THEN !stress/strain-rate relation
                IF(il.EQ.1) THEN
                  FORMAT='(/'' Specify whether fibre elastic'
     '              //' modulus is [1]:'''
     '              //TITLE1(IBEG_T1:IEND_T1)//'/$,''    '',I1)'
                ELSE IF(il.EQ.2) THEN
                  FORMAT='(/'' Specify whether fibre-direction'
     '              //' viscosity is [1]:'''
     '              //TITLE1(IBEG_T1:IEND_T1)//'/$,''    '',I1)'
                ELSE IF(il.EQ.3) THEN
                  FORMAT='(/'' Specify whether transverse dirn'
     '              //' viscosity is [1]:'''
     '              //TITLE1(IBEG_T1:IEND_T1)//'/$,''    '',I1)'
                ELSE IF(il.EQ.4) THEN
                  FORMAT='(/'' Specify whether velocity field'
     '              //' is [1]:'''
     '              //TITLE1(IBEG_T1:IEND_T1)//'/$,''    '',I1)'
                ENDIF !il

              ENDIF !ktyp56

C             Set default value for current parameter
              IF(KTYP56(nr).EQ.0) THEN
                RDEFLT(1)=1.0d0
              ELSE
                IF(KTYP52(nr).EQ.1) K52T=1
                IF(KTYP52(nr).GE.2) K52T=2
                IF(KTYP5F(nr).EQ.1) THEN
                  RDEFLT(1)=RMAT5(il,KTYP56(nr),KTYP55(nr),K52T)
                ELSE
                  RDEFLT(1)=RMAT6(il)
                ENDIF
              ENDIF

              IF(KTYP52(nr).GE.2.AND.KTYP55(nr).EQ.3.AND.
     '          KTYP56(nr).EQ.3) THEN
C               Incomp (+ fluid); fibre strains; pole-zero law.
C               Constit parameters are determined later
                SETCONSTIT=.FALSE.
              ELSE
C               Copy entered params to constit pars as they're entered
                SETCONSTIT=.TRUE.
              ENDIF

              CALL IPMAT5_SINGLE(GRNGLIST,ICQS_SPATIAL,
     '          il,ILPIN,IRCQS_SPATIAL,NBJ,
     '          NEELEM,NELIST,NGLIST,
     '          NMBIN,NPLIST,NPNODE,NQET,NQNE,NQS,NQXI,nr,nx,
     '          CE,CELL_RCQS_VALUE,CGE,
     '          CIN,CP,RCQS_SPATIAL,RDEFLT,
     '          XIG,FORMAT,SETCONSTIT,ERROR,*9999)

            ENDDO !i1
          ENDDO !i2
        ENDDO !i3

! Alex 20-Nov-02: set IT(3,nr) back to zero

        IF(KTYP56(nr).EQ.2.AND.KTYP52(nr).EQ.1) IT(3,nr)=0

        ILTIN=il !set number of entered constitutive law params
        ILT(1,nr,nx)=ILTIN !total number of constitutive parameters

C--------------- Addiional constitutive law params --------------------

        IF((KTYP52(nr).EQ.4.OR.KTYP52(nr).EQ.6).AND.KTYP55(nr).EQ.1.
     '    AND.KTYP56(nr).EQ.1) THEN
C         Compress + fluid with strain energy polynomial in strain invariants
          ILT(1,nr,nx)=ILT(1,nr,nx)+1 !is new total for material consts
          IL_porosity=ILT(1,nr,nx)    !is IL# for porosity
          FORMAT='(/'' Specify porosity [0.1]:'''//
     '      TITLE1(IBEG_T1:IEND_T1)//'/$,''   '',I1)'
C         Set default value for current parameter
          RDEFLT(1)=0.1d0
C         Copy entered params to constit pars as they are entered
          SETCONSTIT=.TRUE.
          CALL IPMAT5_SINGLE(GRNGLIST,ICQS_SPATIAL,IL_porosity,
     '      ILPIN,IRCQS_SPATIAL,NBJ,NEELEM,NELIST,
     '      NGLIST,NMBIN,NPLIST,NPNODE,NQET,NQNE,NQS,NQXI,nr,nx,
     '      CE,CELL_RCQS_VALUE,CGE,
     '      CIN,CP,RCQS_SPATIAL,RDEFLT,
     '      XIG,FORMAT,SETCONSTIT,ERROR,*9999)

          ILT(1,nr,nx)=ILT(1,nr,nx)+1 !is new total for material consts
          IL_compliance=ILT(1,nr,nx)  !is IL# for compliance
          FORMAT='(/'' Specify compliance [0.1]:'''//
     '      TITLE1(IBEG_T1:IEND_T1)//'/$,''   '',I1)'
C         Set default value for current parameter
          RDEFLT(1)=1.d0
C         Copy entered params to constit pars as they are entered
          SETCONSTIT=.TRUE.
          CALL IPMAT5_SINGLE(GRNGLIST,ICQS_SPATIAL,IL_compliance,
     '      ILPIN,IRCQS_SPATIAL,NBJ,NEELEM,NELIST,
     '      NGLIST,NMBIN,NPLIST,NPNODE,NQET,NQNE,NQS,NQXI,
     '      nr,nx,CE,CELL_RCQS_VALUE,CGE,
     '      CIN,CP,RCQS_SPATIAL,RDEFLT,
     '      XIG,FORMAT,SETCONSTIT,ERROR,*9999)
        ENDIF

C--------------- Calculate constitutive law params --------------------

C ***   Determine constitutive parameters from entered parameters
        IF(KTYP52(nr).GE.2.AND.KTYP55(nr).EQ.3.AND.KTYP56(nr).EQ.3) THEN
C         Incomp (+ fluid); fibre strains; pole-zero law
          IF(KTYP5F(nr).EQ.1) THEN
            num_entered=9
          ELSE
            num_entered=30
          ENDIF
          DO il=1,num_entered !Axial parameters are same as CIN(1..9,0,ne/np)
            ILP(il,1,nr,nx)=ILPIN(il)
            IF(ILP(il,1,nr,nx).EQ.1.OR. !constant spatially
     '        ILP(il,1,nr,nx).EQ.2) THEN !defined by elements
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                CE(il,ne)=CIN(il,0,ne)
              ENDDO !noelem (ne)
            ELSE IF(ILP(il,1,nr,nx).EQ.3) THEN !defined by nodes
              NMB(il,nr,nx)=NMBIN(il)
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                CP(il,np)=CIN(il,0,np)
              ENDDO !nonode (np)
            ELSE IF((ILP(il,1,nr,nx).EQ.4).OR.
     '          (ILP(il,1,nr,nx).EQ.5)) THEN !defined by Gauss pr Grid points
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                nb=NBJ(1,ne)
                DO ng=1,NGT(nb)
                  CGE(il,ng,ne)=CIN(il,ng,ne)
                ENDDO !ng
              ENDDO !noelem (ne)
            ENDIF
          ENDDO !il

          IF(KTYP5F(nr).EQ.1) THEN !fibre distribution model

C ***       Pole-zero shear constitutive parameters are determined from
C ***       functionalrelationships found using the
C ***       fibre distribution model
C!!!        NOTE: any changes made to shear constitutive law parameter
C!!!        values below need to reflected in the derivatives
C!!!        calculated in subroutine D_ENERGY (FE50.F), and the updating
C!!!        of material parameters in subroutine EVRESI (FE21.F)
            DO il=10,25,3 !shear coeffs - all 1.0d0 for now
              IF(il.EQ.10.OR.il.EQ.13) THEN
                nm_axial=1
              ELSE IF(il.EQ.16.OR.il.EQ.19) THEN
                nm_axial=4
              ELSE IF(il.EQ.22.OR.il.EQ.25) THEN
                nm_axial=7
              ENDIF
C!!! CS 3/3/2000 Don't think this is right, prevents CP or CGE setup
C            ILP(il,1,nr,nx)=1 !TEMPORARY: shear coeffs const spatially
              ILP(il,1,nr,nx)=ILPIN(nm_axial)
              IF(ILP(il,1,nr,nx).EQ.1.OR. !constant spatially
     '          ILP(il,1,nr,nx).EQ.2) THEN !defined by elements
                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  CE(il,ne)=1.0d0
                ENDDO !noelem (ne)
              ELSE IF(ILP(il,1,nr,nx).EQ.3) THEN !defined by nodes
                NMB(il,nr,nx)=NMBIN(nm_axial)
                DO nonode=1,NPNODE(0,nr)
                  np=NPNODE(nonode,nr)
                  CP(il,np)=1.0d0
                ENDDO !nonode (np)
              ELSE IF((ILP(il,1,nr,nx).EQ.4).OR.
     '          (ILP(il,1,nr,nx).EQ.5)) THEN !defined by Gauss or Grid points
                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  nb=NBJ(1,ne)
                  DO ng=1,NGT(nb)
                    CGE(il,ng,ne)=1.0d0
                  ENDDO !ng
                ENDDO !noelem (ne)
              ENDIF
            ENDDO !il

            DO il=11,26,3 !shear poles - calc'ed from axial poles
              IF(il.EQ.11.OR.il.EQ.14) THEN
                nm_axial=2
              ELSE IF(il.EQ.17.OR.il.EQ.20) THEN
                nm_axial=5
              ELSE IF(il.EQ.23.OR.il.EQ.26) THEN
                nm_axial=8
              ENDIF
              ILP(il,1,nr,nx)=ILPIN(nm_axial)
              IF(ILP(il,1,nr,nx).EQ.1.OR. !constant spatially
     '          ILP(il,1,nr,nx).EQ.2) THEN !defined by elements
                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  CE(il,ne)=2.0d0*CE(nm_axial,ne)/DSQRT(1.0d0+
     '              2.0d0*CE(nm_axial,ne))
                ENDDO !noelem (ne)
              ELSE IF(ILP(il,1,nr,nx).EQ.3) THEN !defined by nodes
                NMB(il,nr,nx)=NMBIN(nm_axial)
                DO nonode=1,NPNODE(0,nr)
                  np=NPNODE(nonode,nr)
                  CP(il,np)=2.0d0*CP(nm_axial,np)/DSQRT(1.0d0+
     '              2.0d0*CP(nm_axial,np))
                ENDDO !nonode (np)
              ELSE IF((ILP(il,1,nr,nx).EQ.4).OR.
     '            (ILP(il,1,nr,nx).EQ.5)) THEN !defined by Gauss or Grid points
                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  nb=NBJ(1,ne)
                  DO ng=1,NGT(nb)
                    CGE(il,ng,ne)=2.0d0*CGE(nm_axial,ng,ne)/DSQRT(1.0d0+
     '              2.0d0*CGE(nm_axial,ng,ne))
                  ENDDO !ng
                ENDDO !noelem (ne)
              ENDIF
            ENDDO !il

            DO il=12,27,3 !shear curvatures - all 2.0d0 for now
              IF(il.EQ.12.OR.il.EQ.15) THEN
                nm_axial=3
              ELSE IF(il.EQ.18.OR.il.EQ.21) THEN
                nm_axial=6
              ELSE IF(il.EQ.24.OR.il.EQ.27) THEN
                nm_axial=9
              ENDIF
C!!! CS 3/3/2000 Don't think this is right, prevents CP or CGE setup
!            ILP(il,1,nr,nx)=1 !TEMP: shear curvatures const spatially
              ILP(il,1,nr,nx)=ILPIN(nm_axial)
              IF(ILP(il,1,nr,nx).EQ.1.OR. !constant spatially
     '          ILP(il,1,nr,nx).EQ.2) THEN !defined by elements
                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  CE(il,ne)=2.0d0
                ENDDO !noelem (ne)
              ELSE IF(ILP(il,1,nr,nx).EQ.3) THEN !defined by nodes
                NMB(il,nr,nx)=NMBIN(nm_axial)
                DO nonode=1,NPNODE(0,nr)
                  np=NPNODE(nonode,nr)
                  CP(il,np)=2.0d0
                ENDDO !nonode (np)
              ELSE IF((ILP(il,1,nr,nx).EQ.4).OR.
     '          (ILP(il,1,nr,nx).EQ.5)) THEN !defined by Gauss or Grid points
                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  nb=NBJ(1,ne)
                  DO ng=1,NGT(nb)
                    CGE(il,ng,ne)=2.0d0
                  ENDDO !ng
                ENDDO !noelem (ne)
              ENDIF
            ENDDO !il
  
            DO il=28,30 !initial strains same as CIN(14..16,0,ne/np)
C!!!          WARNING: the following could possibly change
              nm_extn=il-14
              ILP(il,1,nr,nx)=ILPIN(nm_extn)
              IF(ILP(il,1,nr,nx).EQ.1.OR. !constant spatially
     '          ILP(il,1,nr,nx).EQ.2) THEN !defined by elements
                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  CE(il,ne)=CIN(nm_extn,0,ne)
                ENDDO !noelem (ne)
              ELSE IF(ILP(il,1,nr,nx).EQ.3) THEN !defined by nodes
                NMB(il,nr,nx)=NMBIN(nm_extn)
                DO nonode=1,NPNODE(0,nr)
                  np=NPNODE(nonode,nr)
                  CP(il,np)=CIN(nm_extn,0,np)
                ENDDO !nonode (np)
              ELSE IF((ILP(il,1,nr,nx).EQ.4).OR.
     '          (ILP(il,1,nr,nx).EQ.5)) THEN !defined by Gauss or Grid points
                DO noelem=1,NEELEM(0,nr)
                  ne=NEELEM(noelem,nr)
                  nb=NBJ(1,ne)
                  DO ng=1,NGT(nb)
                    CGE(il,ng,ne)=CIN(nm_extn,ng,ne)
                  ENDDO !ng
                ENDDO !noelem (ne)
              ENDIF
            ENDDO !il
 
          ENDIF !fibre distribution model

          ILT(1,nr,nx)=30 !reset total number of constitutive parameters

        ENDIF !ktyp52,ktyp55,ktyp56

C       Membrane stress model: Membrane thickness
        IF(KTYP51(nr).EQ.4) THEN !membrane
          ILT(1,nr,nx)=ILT(1,nr,nx)+1 !is new total for material consts
          IL_thickness=ILT(1,nr,nx) !is IL# for membrane thickness
          FORMAT='(/'' Specify whether membrane thickness is [1]:'''//
     '      TITLE1(IBEG_T1:IEND_T1)//'/$,''   '',I1)'
C         Set default value for current parameter
          RDEFLT(1)=1.0d0
C         Copy entered params to constit pars as they are entered
          SETCONSTIT=.TRUE.
          CALL IPMAT5_SINGLE(GRNGLIST,ICQS_SPATIAL,IL_thickness,
     '      ILPIN,IRCQS_SPATIAL,NBJ,NEELEM,NELIST,
     '      NGLIST,NMBIN,NPLIST,NPNODE,NQET,NQNE,NQS,NQXI,nr,nx,CE,
     '      CELL_RCQS_VALUE,CGE,
     '      CIN,CP,RCQS_SPATIAL,RDEFLT,
     '      XIG,FORMAT,SETCONSTIT,ERROR,*9999)
        ENDIF !KTYP51(nr)=4 (membrane)

C       String stress model: String crosssection area
        IF(KTYP51(nr).EQ.5) THEN !string
          ILT(1,nr,nx)=ILT(1,nr,nx)+1 !is new total for material consts
          IL_thickness=ILT(1,nr,nx) !is IL# for string crosssection area
          FORMAT='(/'' Specify whether string Xsection area is [1]:'''//
     '      TITLE1(IBEG_T1:IEND_T1)//'/$,''   '',I1)'
C         Set default value for current parameter
          RDEFLT(1)=1.0d0
C         Copy entered params to constit pars as they are entered
          SETCONSTIT=.TRUE.
          CALL IPMAT5_SINGLE(GRNGLIST,ICQS_SPATIAL,IL_thickness,
     '      ILPIN,IRCQS_SPATIAL,NBJ,NEELEM,NELIST,
     '      NGLIST,NMBIN,NPLIST,NPNODE,NQET,NQNE,NQS,NQXI,nr,nx,CE,
     '      CELL_RCQS_VALUE,CGE,
     '      CIN,CP,RCQS_SPATIAL,RDEFLT,
     '      XIG,FORMAT,SETCONSTIT,ERROR,*9999)
        ENDIF !KTYP51(nr)=5 (string)

C       Active stress model: reference SL extension & time-delay
        IF(KTYP53(nr).EQ.3) THEN !active stress included
          ILT(1,nr,nx)=ILT(1,nr,nx)+1 !is new total for material consts
          IL_sarcomere=ILT(1,nr,nx) !is IL# for stress-free SL dist'n
C         Set FORMAT for parameter variation question
          FORMAT='(/'' Specify whether stress-free sarcomere extension'
     '      //' ratios are [1]:'''//TITLE1(IBEG_T1:IEND_T1)
     '      //'/$,''    '',I1)'
C         Set default value for current parameter
          RDEFLT(1)=1.0d0
C         Copy entered params to constit pars as they are entered
          SETCONSTIT=.TRUE.
          CALL IPMAT5_SINGLE(GRNGLIST,ICQS_SPATIAL,IL_sarcomere,
     '      ILPIN,IRCQS_SPATIAL,NBJ,NEELEM,NELIST,
     '      NGLIST,NMBIN,NPLIST,NPNODE,NQET,NQNE,NQS,NQXI,nr,nx,CE,
     '      CELL_RCQS_VALUE,CGE,
     '      CIN,CP,RCQS_SPATIAL,RDEFLT,
     '      XIG,FORMAT,SETCONSTIT,ERROR,*9999)

          ILT(1,nr,nx)=ILT(1,nr,nx)+1 !is new total for material consts
          IL_time_delay=ILT(1,nr,nx) !is IL# for time-delay variable
C         Set FORMAT for parameter variation question
          FORMAT='(/'' Specify whether time delay (in secs) to'
     '      //' activation is [1]:'''//TITLE1(IBEG_T1:IEND_T1)
     '      //'/$,''    '',I1)'
C         Set default value for current parameter
          RDEFLT(1)=0.0d0
C         Copy entered params to constit pars as they are entered
          SETCONSTIT=.TRUE.
          CALL IPMAT5_SINGLE(GRNGLIST,ICQS_SPATIAL,IL_time_delay,
     '      ILPIN,IRCQS_SPATIAL,NBJ,NEELEM,NELIST,
     '      NGLIST,NMBIN,NPLIST,NPNODE,NQET,NQNE,NQS,NQXI,nr,nx,CE,
     '      CELL_RCQS_VALUE,CGE,
     '      CIN,CP,RCQS_SPATIAL,RDEFLT,
     '      XIG,FORMAT,SETCONSTIT,ERROR,*9999)
        ENDIF

C       Fluid shift model: Fluid conductivity
        IF(KTYP52(nr).EQ.3) THEN !incompressible + fluid
          ILT(1,nr,nx)=ILT(1,nr,nx)+1 !is new total for material consts
          IL_fluid_conductivity=ILT(1,nr,nx) !is IL# of fluid conduct'y
C!!!     May need to extend this param to all 3 dirs (not just in Xi3).
C         Set FORMAT for parameter variation question
          FORMAT='(/'' Specify whether the through wall fluid '
     '      //'conductivity is [1]:'''//TITLE1(IBEG_T1:IEND_T1)
     '      //'/$,''    '',I1)'
C         Set default value for current parameter
          RDEFLT(1)=1.0d0
C         Copy entered params to constit pars as they are entered
          SETCONSTIT=.TRUE.
          CALL IPMAT5_SINGLE(GRNGLIST,ICQS_SPATIAL,
     '      IL_fluid_conductivity,
     '      ILPIN,IRCQS_SPATIAL,NBJ,
     '      NEELEM,NELIST,NGLIST,
     '      NMBIN,NPLIST,NPNODE,NQET,NQNE,NQS,NQXI,nr,nx,CE,
     '      CELL_RCQS_VALUE,CGE,
     '      CIN,CP,RCQS_SPATIAL,RDEFLT,
     '      XIG,FORMAT,SETCONSTIT,ERROR,*9999)
        ENDIF

C       density for gravity term
        ILT(1,nr,nx)=ILT(1,nr,nx)+1 !is new total for material consts
        IL_density=ILT(1,nr,nx) !is IL# for density
        FORMAT='(/'' Specify whether the density is [1]:'''//
     '    TITLE1(IBEG_T1:IEND_T1)//'/$,''   '',I1)'
C       Set default value for current parameter
        RDEFLT(1)=0.0d0
C       Copy entered params to constit pars as they are entered
        SETCONSTIT=.TRUE.
        CALL IPMAT5_SINGLE(GRNGLIST,ICQS_SPATIAL,IL_density,
     '    ILPIN,IRCQS_SPATIAL,NBJ,NEELEM,NELIST,
     '    NGLIST,NMBIN,NPLIST,NPNODE,NQET,NQNE,NQS,NQXI,nr,nx,CE,
     '    CELL_RCQS_VALUE,CGE,
     '    CIN,CP,RCQS_SPATIAL,RDEFLT,XIG,FORMAT,SETCONSTIT,ERROR,*9999)

        IF(KTYP52(nr).GE.2.AND.KTYP55(nr).EQ.1.AND.KTYP54(nr).EQ.1.
     &    AND.KTYP52(nr).EQ.6) THEN  !compr + fluid for lung
          FORMAT='($,'' Specify the limit on ratio of '//
     '      'deformed to undeformed volume [1.0]: '',D12.4)'
          RDEFLT(1)=1.d0
          IF(IOTYPE.EQ.3) RDATA(1)=RATIO_LIMIT
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,0.d0,RMAX,INFO,ERROR,*9999)
          RATIO_LIMIT=RDATA(1)
          FORMAT='($,'' Specify the penalty stiffness [1000]:'',D12.4)'
          RDEFLT(1)=1.d4
          IF(IOTYPE.EQ.3) RDATA(1)=RATIO_STIFFNESS
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,0.d0,RMAX,INFO,ERROR,*9999)
          RATIO_STIFFNESS=RDATA(1)
        ENDIF

      ELSE IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN
C       const vol constraint for cavity regions
        ILT(1,nr,nx)=1 !total number of constitutive parameters
C       Set FORMAT for parameter variation question
        FORMAT='('' Specify whether the spring stiffness for the '
     '    //'cavity is [1]:'''//TITLE1(IBEG_T1:IEND_T1)
     '    //'/$,''    '',I1)'
C       Set default value for current parameter
        RDEFLT(1)=1.0d0
C       Copy entered params to constit pars as they are entered
        SETCONSTIT=.TRUE.
        CALL IPMAT5_SINGLE(GRNGLIST,ICQS_SPATIAL,ILT(1,nr,nx),
     '    ILPIN,IRCQS_SPATIAL,NBJ,NEELEM,NELIST,
     '    NGLIST,NMBIN,NPLIST,NPNODE,NQET,NQNE,NQS,NQXI,nr,nx,CE,
     '    CELL_RCQS_VALUE,CGE,
     '    CIN,CP,RCQS_SPATIAL,RDEFLT,
     '    XIG,FORMAT,SETCONSTIT,ERROR,*9999)
      ELSE !incorrect problem type
        ERROR=' >>Problem type is incorrect'
        GOTO 9999
      ENDIF

      WRITE(CHAR7,'(I3)') ILT(1,nr,nx)
      CALL ASSERT(ILT(1,nr,nx).LE.NMM,'>>NMM too small. '
     '  //'Must be at least'//CHAR7,ERROR,*9999)

      CALL EXITS('IPMAT5')
      RETURN
 9999 CALL ERRORS('IPMAT5',ERROR)
      CALL EXITS('IPMAT5')
      RETURN 1
      END


