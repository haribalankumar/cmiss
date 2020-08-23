      SUBROUTINE OPMAT5(ILPIN,NBJ,NEELEM,NMBIN,NPNODE,nr,nx,CE,CGE,
     '  CIN,CP,CONSTIT,ERROR,*)

C#### Subroutine: OPMAT5
C###  Description:
C###    OPMAT5 outputs passive model equation parameters.

      IMPLICIT NONE
      INCLUDE 'b13.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ipma50.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'titl50.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER ILPIN(NMM),NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NMBIN(NMM),
     '  NPNODE(0:NP_R_M,0:NRM),nr,nx
      REAL*8 CE(NMM,NEM),CGE(NMM,NGM,NEM),
     '  CIN(NMM,0:NGM,NNEPM),CP(NMM,NPM)
      CHARACTER ERROR*(*)
      LOGICAL CONSTIT
!     Local Variables
      INTEGER i1,i2,i3,IBEG1,IEND1,il,IL_TOTAL,K52T,nb,ne,ng,
     '  noelem,nonode
      CHARACTER CHAR1*40,CONSTIT_NAME(31)*31,FORMAT_local*132,
     '  TITLE5(5)*35

      DATA TITLE5 /'constant spatially                 ',
     '             'piecewise constant (element params)',
     '             'piecewise linear (nodal parameters)',
     '             'by Gauss points                    ',
     '             'by Grid points                    '/
      DATA CONSTIT_NAME /'E(1,1) term coefficient (kPa)',
     '                   'E(1,1) term pole             ',
     '                   'E(1,1) term curvature        ',
     '                   'E(2,2) term coefficient (kPa)',
     '                   'E(2,2) term pole             ',
     '                   'E(2,2) term curvature        ',
     '                   'E(3,3) term coefficient (kPa)',
     '                   'E(3,3) term pole             ',
     '                   'E(3,3) term curvature        ',
     '                   'E(1,2) term coefficient (kPa)',
     '                   'E(1,2) term pole             ',
     '                   'E(1,2) term curvature        ',
     '                   'E(1,3) term coefficient (kPa)',
     '                   'E(1,3) term pole             ',
     '                   'E(1,3) term curvature        ',
     '                   'E(2,3) term coefficient (kPa)',
     '                   'E(2,3) term pole             ',
     '                   'E(2,3) term curvature        ',
     '                   'E(2,1) term coefficient (kPa)',
     '                   'E(2,1) term pole             ',
     '                   'E(2,1) term curvature        ',
     '                   'E(3,1) term coefficient (kPa)',
     '                   'E(3,1) term pole             ',
     '                   'E(3,1) term curvature        ',
     '                   'E(3,2) term coefficient (kPa)',
     '                   'E(3,2) term pole             ',
     '                   'E(3,2) term curvature        ',
     '                   'Initial fibre extension ratio',
     '                   'Initial sheet extension ratio',
     '                   'Initial sheet-norm extn ratio',
     &                   'Density                      '/

      CALL ENTERS('OPMAT5',*9999)
      CALL ASSERT(KTYP52(nr).GT.0,'>>Material parameters not defined',
     '  ERROR,*9999)

      IF(ITYP2(nr,nx).EQ.2) THEN !finite elasticity

        WRITE(OP_STRING,'(''               Material is: '',A/'
     '    //'''  Stresses are referred to: '',A)')
     '    TITL52(KTYP52(nr)),TITL53(KTYP53(nr))
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(KTYP54(nr).EQ.1) THEN
          IF(KTYP52(nr).EQ.1) K52T=1
          IF(KTYP52(nr).GE.2) K52T=2
          WRITE(OP_STRING,
     '      '(''  Strain energy: is a function of '',A,/15X,'
     '      //'''  and is '',A)')
     '      TITL55(KTYP55(nr)),TITL56(KTYP56(nr),K52T,KTYP55(nr))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE IF(KTYP54(nr).EQ.2) THEN
          WRITE(OP_STRING,
     '      '(''  Constitutive law is a function of strains and'
     '      //' strain rates'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE IF(KTYP54(nr).EQ.3) THEN
C         creep function output
        ENDIF !ktyp54

        IF(KTYP56(nr).EQ.0) THEN !constants

        ELSE IF(KTYP56(nr).EQ.1) THEN !polynomial function
          IF(KTYP55(nr).EQ.1) THEN
            FORMAT_local='(/I3,'') The coefficient of (I1-3)^'',I1,'
     '        //'''*(I2-3)^'',I1,''*(I3-1)^'',I1,'' is '',/4X,A)'
          ELSE IF(KTYP55(nr).EQ.2) THEN
            FORMAT_local='(/I3,'') The coefficient of (L1-1)^'',I1,'
     '        //'''*(L2-1)^'',I1,''*(L3-1)^'',I1,'' is '',/4X,A)'
          ELSE IF(KTYP55(nr).EQ.3) THEN
            FORMAT_local='(/I3,'') Parameter '',I1,'' in term '',I1,'
     '        //''' is:'',1X,A)'
          ENDIF

        ELSE IF(KTYP56(nr).EQ.2) THEN !exponential
          IF(KTYP55(nr).EQ.1) THEN
            IF(KTYP52(nr).EQ.1) THEN
              FORMAT_local='(/I3,'') The Blatz-Ko material '
     '          //'parameter ('',I1,'') is:''/4X,A)'
            ELSE IF(KTYP52(nr).GE.2) THEN
              FORMAT_local='(/I3,'') The exponential material '
     '          //'parameter ('',I1,'') is:''/4X,A)'
            ENDIF
          ELSE IF(KTYP55(nr).EQ.2) THEN
            FORMAT_local='(/I3,'') The Ogden material parameter ('',I1,'
     '        //''') is:''/4X,A)'
          ELSE IF(KTYP55(nr).EQ.3) THEN
            FORMAT_local='(/I3,'') The exponential material parameter '
     '        //'('',I1,'') is:''/4X,A)'
          ENDIF

        ELSE IF(KTYP56(nr).EQ.3) THEN !pole-zero
          FORMAT_local='(/I3,'') The user-defined material '
     '      //'parameter ('',I1,'') is '',A)'

        ELSE IF(KTYP56(nr).EQ.6) THEN !stress/strain-rate relation
          FORMAT_local='(/I3,'') Fibre direction viscosity is '',A)'
          FORMAT_local='(/I3,'') Transverse dirn viscosity is '',A)'

        ENDIF
        IF(KTYP52(nr).EQ.1) K52T=1
        IF(KTYP52(nr).GE.2) K52T=2

        IF(.NOT.CONSTIT) THEN !print material params that were entered
C         Entered material parameters are stored in CIN
          IF(KTYP52(nr).EQ.2.OR.KTYP52(nr).EQ.3.
     '      AND.KTYP55(nr).EQ.3) THEN
C           Incomp (+ fluid); fibre strains
            WRITE(OP_STRING,'(/'' Fibre Family parameters:'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE !all other materials
            WRITE(OP_STRING,'(/'' Entered material parameters:'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
          il=0
          DO i3=0,IT(3,nr)
            DO i2=0,IT(2,nr)
              DO i1=0,IT(1,nr)
                il=il+1
                IF(KTYP56(nr).EQ.0) THEN
                  WRITE(OP_STRING,FMT=FORMAT_local)
     '              il,TITLE5(ILPIN(il))
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                ELSE IF(KTYP56(nr).EQ.1.AND.KTYP55(nr).LE.2) THEN
                  WRITE(OP_STRING,FMT=FORMAT_local)
     '              il,i1,i2,i3,TITLE5(ILPIN(il))
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                ELSE IF(KTYP56(nr).EQ.1.AND.KTYP55(nr).EQ.3) THEN
                  WRITE(OP_STRING,FMT=FORMAT_local)
     '              il,i2+1,i3+1,TITLE5(ILPIN(il))
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                ELSE
                  CHAR1=CMAT5(il,KTYP56(nr),KTYP55(nr),K52T)
                  CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
                  FORMAT_local='(/I3,'') '//CHAR1(IBEG1:IEND1)
     '              //' is:''/4X,A)'
                  WRITE(OP_STRING,FMT=FORMAT_local) il,TITLE5(ILPIN(il))
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                ENDIF
                IF(ILPIN(il).EQ.1) THEN !constant spatially
                  WRITE(OP_STRING,'(4X,''Value='',D11.4)')
     '              CIN(il,0,NEELEM(1,nr))
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                ELSE IF(ILPIN(il).EQ.2) THEN !defined by elements
                  WRITE(OP_STRING,
     '              '(4X,''Element values:'',:/(5(1X,I3,'': '','
     '              //'D10.3)))') (NEELEM(noelem,nr),
     '              CIN(il,0,NEELEM(noelem,nr)),noelem=1,NEELEM(0,nr))
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                ELSE IF(ILPIN(il).EQ.3) THEN !defined by nodes
                  WRITE(OP_STRING,'(4X,''Using basis '',I3,'
     '              //''' for nodal values:'','
     '              //':/(5(1X,I3,'': '',D10.3)))')
     '              NMBIN(il),(NPNODE(nonode,nr),
     '              CIN(il,0,NPNODE(nonode,nr)),nonode=1,NPNODE(0,nr))
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                ELSE IF((ILPIN(il).EQ.4).OR.
     '             (ILPIN(il).EQ.5)) THEN !defined by Gauss or Grid pts
C new CS Fri May 28 1999 now implemented
C                  CALL ASSERT(.FALSE.,' >>> Gauss pt mat param '
C     '              //'variation not implemented',ERROR,*9999)
                  DO noelem=1,NEELEM(0,nr)
                    ne=NEELEM(noelem,nr)
                    nb=NBJ(1,ne)
                    WRITE(OP_STRING,'(4X,'' Element '',I5,'' Basis '',
     '                I2)') ne,nb
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                    DO ng=1,NGT(nb)
                      WRITE(OP_STRING,'(6X,''Guass Point '',I2,
     '                  '' value :'',D10.3)')
     '                  ng,CIN(il,ng,ne)
                      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                    ENDDO !ng
                  ENDDO !noelem
                ENDIF
              ENDDO !i1
            ENDDO !i2
          ENDDO !i3

        ELSE IF(CONSTIT) THEN !print constitutive parameters
C         CE/CP/CGE store constitutive parameters that have been
C         determined from entered parameters CIN
          WRITE(OP_STRING,'(/'' Constitutive parameters:'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IL_TOTAL=ILT(1,nr,nx)
C         some params output separately below so subtract from total
          IF(KTYP51(nr).EQ.4) IL_TOTAL=IL_TOTAL-1 !membrane thickness
          IF(KTYP53(nr).EQ.3) IL_TOTAL=IL_TOTAL-2 !active params
          IF(KTYP52(nr).EQ.3) IL_TOTAL=IL_TOTAL-1 !fluid conductivity
          DO il=1,IL_TOTAL !loop over all constitutive parameters
            IF(KTYP56(nr).EQ.1.AND.KTYP55(nr).LE.2.OR.
     '        KTYP56(nr).EQ.1.AND.KTYP55(nr).EQ.3) THEN
              CALL ASSERT(.FALSE.,' >>Omit "constitutive" '
     '          //'option to list the material parameters',ERROR,*9999)
            ELSE IF(KTYP52(nr).GE.2.AND.KTYP55(nr).EQ.3.AND.
     '          KTYP56(nr).EQ.3) THEN
C             Incomp (+ fluid); fibre strains; pole-zero law
              CHAR1=CONSTIT_NAME(il)
            ELSE
              CHAR1=CMAT5(il,KTYP56(nr),KTYP55(nr),K52T)
            ENDIF
            CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
            FORMAT_local='(/I3,'') '//CHAR1(IBEG1:IEND1)//' is:''/4X,A)'
            WRITE(OP_STRING,FMT=FORMAT_local) il,TITLE5(ILP(il,1,nr,nx))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            IF(ILP(il,1,nr,nx).EQ.1) THEN !constant spatially
              WRITE(OP_STRING,'(4X,''Value='',D11.4)')
     '          CE(il,NEELEM(1,nr))
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ELSE IF(ILP(il,1,nr,nx).EQ.2) THEN !defined by elements
              WRITE(OP_STRING,
     '          '(4X,''Element values:'',:/(5(1X,I3,'': '','
     '          //'D10.3)))') (NEELEM(noelem,nr),
     '          CE(il,NEELEM(noelem,nr)),noelem=1,NEELEM(0,nr))
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ELSE IF(ILP(il,1,nr,nx).EQ.3) THEN !defined by nodes
              WRITE(OP_STRING,'(4X,''Using basis '',I3,'' for nodal '
     '          //'values:'',:/(5(1X,I3,'': '',D10.3)))')
     '          NMB(il,nr,nx),(NPNODE(nonode,nr),
     '          CP(il,NPNODE(nonode,nr)),nonode=1,NPNODE(0,nr))
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ELSE IF((ILP(il,1,nr,nx).EQ.4).OR.
     '          (ILP(il,1,nr,nx).EQ.5)) THEN !by Gauss or Grid points
C new CS Fri May 28 1999 now implemented
C              CALL ASSERT(.FALSE.,' >>> Gauss pt mat param variation '
C     '          //'not implemented',ERROR,*9999)
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                nb=NBJ(1,ne)
                WRITE(OP_STRING,'(4X,'' Element '',I5,'' Basis '',
     '            I2)') ne,nb
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                DO ng=1,NGT(nb)
                  WRITE(OP_STRING,'(6X,''Guass Point '',I2,
     '              '' value :'',D10.3)')
     '              ng,CGE(il,ng,ne)
                  CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                ENDDO !ng
              ENDDO !noelem
            ENDIF
          ENDDO !il
        ENDIF !ktyp56

C       Membrane/string model: Membrane thickness or string Xsection area
        IF(KTYP51(nr).EQ.4.OR.KTYP51(nr).EQ.5) THEN !membrane or string
          il=IL_thickness
          IF(KTYP51(nr).EQ.4) THEN      !membrane
            WRITE(OP_STRING,'(/I3,'') Membrane thickness '
     '        //'is:''/4X,A)') il,TITLE5(ILP(il,1,nr,nx))
          ELSE IF(KTYP51(nr).EQ.5) THEN !string
            WRITE(OP_STRING,'(/I3,'') String Xsection area '
     '        //'is:''/4X,A)') il,TITLE5(ILP(il,1,nr,nx))
          ENDIF
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IF(ILP(il,1,nr,nx).EQ.1) THEN !constant spatially
            WRITE(OP_STRING,'(4X,''Value='',D11.4)') CE(il,NEELEM(1,nr))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE IF(ILP(il,1,nr,nx).EQ.2) THEN !defined by elements
            WRITE(OP_STRING,'(4X,''Element values:'','
     '        //'(/5(1X,I3,'': '',D10.3)))')
     '        (NEELEM(noelem,nr),CE(il,NEELEM(noelem,nr)),
     '        noelem=1,NEELEM(0,nr))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE IF(ILP(il,1,nr,nx).EQ.3) THEN !defined by nodes
            WRITE(OP_STRING,'(4X,''Using basis '',I3,'' for nodal '
     '        //'values:'',(/5(1X,I3,'': '',D10.3)))')
     '        NMB(il,nr,nx),(NPNODE(nonode,nr),
     '        CP(il,NPNODE(nonode,nr)),nonode=1,NPNODE(0,nr))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE IF((ILP(il,1,nr,nx).EQ.4).OR.
     '        (ILP(il,1,nr,nx).EQ.5)) THEN !by Gauss or Grid points
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              nb=NBJ(1,ne)
              WRITE(OP_STRING,'(4X,'' Element '',I5,'' Basis '',I2)')
     '          ne,nb
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              DO ng=1,NGT(nb)
                WRITE(OP_STRING,'(6X,''Guass Point '',I2,'' value :'','
     '            //'D10.3)') ng,CGE(il,ng,ne)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDDO !ng
            ENDDO !noelem
          ENDIF !ILP
        ENDIF !ktyp51

C       Compress + fluid with strain energy polynomial in strain invariants
        IF((KTYP52(nr).EQ.4.OR.KTYP52(nr).EQ.6).AND.KTYP55(nr).EQ.1.
     '    AND.KTYP56(nr).EQ.1) THEN
          il=IL_porosity
          WRITE(OP_STRING,'(/I3,'') Porosity '
     '        //'is:''/4X,A)') il,TITLE5(ILP(il,1,nr,nx))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IF(ILP(il,1,nr,nx).EQ.1) THEN !constant spatially
            WRITE(OP_STRING,'(4X,''Value='',D11.4)') CE(il,NEELEM(1,nr))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE IF(ILP(il,1,nr,nx).EQ.2) THEN !defined by elements
            WRITE(OP_STRING,'(4X,''Element values:'','
     '        //'(/5(1X,I3,'': '',D10.3)))')
     '        (NEELEM(noelem,nr),CE(il,NEELEM(noelem,nr)),
     '        noelem=1,NEELEM(0,nr))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE IF(ILP(il,1,nr,nx).EQ.3) THEN !defined by nodes
            WRITE(OP_STRING,'(4X,''Using basis '',I3,'' for nodal '
     '        //'values:'',(/5(1X,I3,'': '',D10.3)))')
     '        NMB(il,nr,nx),(NPNODE(nonode,nr),
     '        CP(il,NPNODE(nonode,nr)),nonode=1,NPNODE(0,nr))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE IF((ILP(il,1,nr,nx).EQ.4).OR.
     '        (ILP(il,1,nr,nx).EQ.5)) THEN !by Gauss or Grid points
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              nb=NBJ(1,ne)
              WRITE(OP_STRING,'(4X,'' Element '',I5,'' Basis '',I2)')
     '          ne,nb
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              DO ng=1,NGT(nb)
                WRITE(OP_STRING,'(6X,''Guass Point '',I2,'' value :'','
     '            //'D10.3)') ng,CGE(il,ng,ne)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDDO !ng
            ENDDO !noelem
          ENDIF !ILP
          il=IL_compliance
          WRITE(OP_STRING,'(/I3,'') Compliance '
     '        //'is:''/4X,A)') il,TITLE5(ILP(il,1,nr,nx))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IF(ILP(il,1,nr,nx).EQ.1) THEN !constant spatially
            WRITE(OP_STRING,'(4X,''Value='',D11.4)') CE(il,NEELEM(1,nr))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE IF(ILP(il,1,nr,nx).EQ.2) THEN !defined by elements
            WRITE(OP_STRING,'(4X,''Element values:'','
     '        //'(/5(1X,I3,'': '',D10.3)))')
     '        (NEELEM(noelem,nr),CE(il,NEELEM(noelem,nr)),
     '        noelem=1,NEELEM(0,nr))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE IF(ILP(il,1,nr,nx).EQ.3) THEN !defined by nodes
            WRITE(OP_STRING,'(4X,''Using basis '',I3,'' for nodal '
     '        //'values:'',(/5(1X,I3,'': '',D10.3)))')
     '        NMB(il,nr,nx),(NPNODE(nonode,nr),
     '        CP(il,NPNODE(nonode,nr)),nonode=1,NPNODE(0,nr))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE IF((ILP(il,1,nr,nx).EQ.4).OR.
     '        (ILP(il,1,nr,nx).EQ.5)) THEN !by Gauss or Grid points
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              nb=NBJ(1,ne)
              WRITE(OP_STRING,'(4X,'' Element '',I5,'' Basis '',I2)')
     '          ne,nb
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              DO ng=1,NGT(nb)
                WRITE(OP_STRING,'(6X,''Guass Point '',I2,'' value :'','
     '            //'D10.3)') ng,CGE(il,ng,ne)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDDO !ng
            ENDDO !noelem
          ENDIF !ILP
        ENDIF !Compress + fluid

C       Active stress model: reference SL extension & time-delay
        IF(KTYP53(nr).EQ.3) THEN !active stress included
          il=IL_sarcomere
          WRITE(OP_STRING,'(/I3,'') Stress-free SL extensions '
     '      //'is:''/4X,A)') il,TITLE5(ILP(il,1,nr,nx))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IF(ILP(il,1,nr,nx).EQ.1) THEN !constant spatially
            WRITE(OP_STRING,'(4X,''Value='',D11.4)') CE(il,NEELEM(1,nr))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE IF(ILP(il,1,nr,nx).EQ.2) THEN !defined by elements
            WRITE(OP_STRING,'(4X,''Element values:'','
     '        //'(/5(1X,I3,'': '',D10.3)))')
     '        (NEELEM(noelem,nr),CE(il,NEELEM(noelem,nr)),
     '        noelem=1,NEELEM(0,nr))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE IF(ILP(il,1,nr,nx).EQ.3) THEN !defined by nodes
            WRITE(OP_STRING,'(4X,''Using basis '',I3,'' for nodal '
     '        //'values:'',(/5(1X,I3,'': '',D10.3)))')
     '        NMB(il,nr,nx),(NPNODE(nonode,nr),
     '        CP(il,NPNODE(nonode,nr)),nonode=1,NPNODE(0,nr))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE IF((ILP(il,1,nr,nx).EQ.4).OR.
     '        (ILP(il,1,nr,nx).EQ.5)) THEN !by Gauss or Grid  points
C new CS Fri May 28 1999 now implemented
C            CALL ASSERT(.FALSE.,' >>> Gauss pt mat param variation '
C     '        //'not implemented',ERROR,*9999)
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              nb=NBJ(1,ne)
              WRITE(OP_STRING,'(4X,'' Element '',I5,'' Basis '',
     '          I2)') ne,nb
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              DO ng=1,NGT(nb)
                WRITE(OP_STRING,'(6X,''Guass Point '',I2,
     '            '' value :'',D10.3)')
     '            ng,CGE(il,ng,ne)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDDO !ng
            ENDDO !noelem
          ENDIF

          il=IL_time_delay
          WRITE(OP_STRING,'(/I3,'') Time-delay to activation '
     '      //'is:''/4X,A)') il,TITLE5(ILP(il,1,nr,nx))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IF(ILP(il,1,nr,nx).EQ.1) THEN !constant spatially
            WRITE(OP_STRING,'(4X,''Value='',D11.4)') CE(il,NEELEM(1,nr))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE IF(ILP(il,1,nr,nx).EQ.2) THEN !defined by elements
            WRITE(OP_STRING,'(4X,''Element values:'','
     '        //'(/5(1X,I3,'': '',D10.3)))')
     '        (NEELEM(noelem,nr),CE(il,NEELEM(noelem,nr)),
     '        noelem=1,NEELEM(0,nr))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE IF(ILP(il,1,nr,nx).EQ.3) THEN !defined by nodes
            WRITE(OP_STRING,'(4X,''Using basis '',I3,'' for nodal '
     '        //'values:'',(/5(1X,I3,'': '',D10.3)))')
     '        NMB(il,nr,nx),(NPNODE(nonode,nr),
     '        CP(il,NPNODE(nonode,nr)),nonode=1,NPNODE(0,nr))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE IF((ILP(il,1,nr,nx).EQ.4).OR.
     '        (ILP(il,1,nr,nx).EQ.5)) THEN !by Gauss or Grid  points
C new CS Fri May 28 1999 now implemented
C            CALL ASSERT(.FALSE.,' >>> Gauss pt mat param variation '
C     '        //'not implemented',ERROR,*9999)
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              nb=NBJ(1,ne)
              WRITE(OP_STRING,'(4X,'' Element '',I5,'' Basis '',
     '          I2)') ne,nb
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              DO ng=1,NGT(nb)
                WRITE(OP_STRING,'(6X,''Guass Point '',I2,
     '            '' value :'',D10.3)')
     '            ng,CGE(il,ng,ne)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDDO !ng
            ENDDO !noelem
          ENDIF
        ENDIF !KTYP53(nr)=3 (active stress)

C       Fluid shift model: Fluid conductivity
        IF(KTYP52(nr).EQ.3) THEN !incompressible + fluid
          il=IL_fluid_conductivity
          WRITE(OP_STRING,'(/I3,'') Through wall fluid conductivity '
     '      //'is:''/4X,A)') il,TITLE5(ILP(il,1,nr,nx))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IF(ILP(il,1,nr,nx).EQ.1) THEN !constant spatially
            WRITE(OP_STRING,'(4X,''Value='',D11.4)') CE(il,NEELEM(1,nr))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE IF(ILP(il,1,nr,nx).EQ.2) THEN !defined by elements
            WRITE(OP_STRING,'(4X,''Element values:'','
     '        //'(/5(1X,I3,'': '',D10.3)))')
     '        (NEELEM(noelem,nr),CE(il,NEELEM(noelem,nr)),
     '        noelem=1,NEELEM(0,nr))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE IF(ILP(il,1,nr,nx).EQ.3) THEN !defined by nodes
            WRITE(OP_STRING,'(4X,''Using basis '',I3,'' for nodal '
     '        //'values:'',(/5(1X,I3,'': '',D10.3)))')
     '        NMB(il,nr,nx),(NPNODE(nonode,nr),
     '        CP(il,NPNODE(nonode,nr)),nonode=1,NPNODE(0,nr))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE IF((ILP(il,1,nr,nx).EQ.4).OR.
     '        (ILP(il,1,nr,nx).EQ.5)) THEN !by Gauss or Grid points
C new CS Fri May 28 1999 now implemented
C            CALL ASSERT(.FALSE.,' >>> Gauss pt mat param variation '
C     '        //'not implemented',ERROR,*9999)
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              nb=NBJ(1,ne)
              WRITE(OP_STRING,'(4X,'' Element '',I5,'' Basis '',
     '          I2)') ne,nb
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              DO ng=1,NGT(nb)
                WRITE(OP_STRING,'(6X,''Guass Point '',I2,
     '            '' value :'',D10.3)')
     '            ng,CGE(il,ng,ne)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDDO !ng
            ENDDO !noelem
          ENDIF
        ENDIF !KTYP52(nr)=3 (incomp+fluid)

      ELSE IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).EQ.3) THEN
C       const vol constraint for cavity regions
        WRITE(OP_STRING,'('' Constant volume cavity:'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        il=ILT(1,nr,nx)
        WRITE(OP_STRING,'(/I3,'') Spring stiffness for cavity '
     '    //'is:''/4X,A)') il,TITLE5(ILP(il,1,nr,nx))
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(ILP(il,1,nr,nx).EQ.1) THEN !constant spatially
          WRITE(OP_STRING,'(4X,''Value='',D11.4)') CE(il,NEELEM(1,nr))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE IF(ILP(il,1,nr,nx).EQ.2) THEN !defined by elements
          WRITE(OP_STRING,'(4X,''Element values:'','
     '      //'(/5(1X,I3,'': '',D10.3)))')
     '      (NEELEM(noelem,nr),CE(il,NEELEM(noelem,nr)),
     '      noelem=1,NEELEM(0,nr))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE IF(ILP(il,1,nr,nx).EQ.3) THEN !defined by nodes
          WRITE(OP_STRING,'(4X,''Using basis '',I3,'' for nodal '
     '      //'values:'',(/5(1X,I3,'': '',D10.3)))')
     '      NMB(il,nr,nx),(NPNODE(nonode,nr),
     '      CP(il,NPNODE(nonode,nr)),nonode=1,NPNODE(0,nr))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE IF((ILP(il,1,nr,nx).EQ.4).OR.
     '      (ILP(il,1,nr,nx).EQ.5)) THEN !by Gauss or Grid points
C new CS Fri May 28 1999 now implemented
C            CALL ASSERT(.FALSE.,' >>> Gauss pt mat param variation '
C     '        //'not implemented',ERROR,*9999)
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              nb=NBJ(1,ne)
              WRITE(OP_STRING,'(4X,'' Element '',I5,'' Basis '',
     '          I2)') ne,nb
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              DO ng=1,NGT(nb)
                WRITE(OP_STRING,'(6X,''Guass Point '',I2,
     '            '' value :'',D10.3)')
     '            ng,CGE(il,ng,ne)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDDO !ng
            ENDDO !noelem
        ENDIF

      ELSE !incorrect problem type
        ERROR=' >>Problem type is incorrect'
        GOTO 9999
      ENDIF

      CALL EXITS('OPMAT5')
      RETURN
 9999 CALL ERRORS('OPMAT5',ERROR)
      CALL EXITS('OPMAT5')
      RETURN 1
      END


