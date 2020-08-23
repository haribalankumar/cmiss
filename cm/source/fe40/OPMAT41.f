      SUBROUTINE OPMAT41(ie,LIST,L1,L2,NPNODE,nr,nx,CE,CP,ERROR,*)

C#### Subroutine: OPMAT41
C###  Description:
C###    OPMAT41 outputs material parameters.

      IMPLICIT NONE
      INCLUDE 'b13.cmn'
      INCLUDE 'b40.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp40.cmn'
      INCLUDE 'titl40.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER ie,LIST(0:NLISTM),l1,l2,NPNODE(0:NP_R_M,0:NRM),nr,nx
      REAL*8 CE(NMM,NEM),CP(NMM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER il,l,nonode
      CHARACTER TITLE*40,TITLE5(4)*40

      DATA TITLE5 /'constant spatially                 ',
     '             'piecewise constant (element params)',
     '             'piecewise linear (nodal parameters)',
     '             'defined by Gauss points            '/

      CALL ENTERS('OPMAT41',*9999)
      IF(ie.EQ.3) THEN
        WRITE(OP_STRING,'(4X,''Cross-section is '',A)')
     '    TITL45(KTYP45)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ELSE
        WRITE(OP_STRING,'(4X,''Material is '',A)') TITL44(IMT(ie))
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF
      DO il=1,ILT(ie,nr,nx)
        IF(il.LE.NMPP(ie,IMT(ie))) THEN     !elastic material params
          TITLE=TITL46(il,IMT(ie))
        ELSE IF(il.LE.NMP(ie)) THEN         !thermal material params
          TITLE=TITL41(il-NMPP(ie,IMT(ie)),IMT(ie))
        ELSE IF(il.LE.NMP(ie)+NGP(ie)) THEN !geometric params
          TITLE=TITL47(il-NMP(ie),ie)
        ELSE IF(il.LE.NMP(ie)+NGP(ie)+NLP(ie)) THEN !element loads
          TITLE=TITL48(il-NMP(ie)-NGP(ie),ie)
        ELSE                                !mass density
          TITLE='Density (kg/m^3)'
        ENDIF
        WRITE(OP_STRING,
     '    '(2X,I2,'')'',A,'' is '',A) ') il,TITLE,
     '    TITLE5(ILP(il,ie,nr,nx))
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(ILP(il,ie,nr,nx).EQ.1) THEN
          IF(KTYP14.EQ.0) THEN        !constant spatially
            WRITE(OP_STRING,'(5X,''value='',E11.4)') CE(il,LIST(l1))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE IF(KTYP14.EQ.1) THEN
            WRITE(OP_STRING,
     '        '(5X,''value='',E11.4,'' increment='',E11.4)')
     '        CE(il,LIST(l1)),CE(il+ILT(ie,nr,nx),LIST(l1))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ELSE IF(ILP(il,ie,nr,nx).EQ.2) THEN !defined by elements
          WRITE(OP_STRING,'(5X,''element values:'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(5(3X,I3,'': '',E10.3))')
     '      (LIST(l),CE(IL,LIST(l)),l=l1,l2)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IF(KTYP14.EQ.1) THEN
            WRITE(OP_STRING,'(5X,''element increments:'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(5(3X,I3,'': '',E10.3))')
     '        (LIST(l),CE(il+ILT(ie,nr,nx),LIST(l)),l=l1,l2)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ELSE IF(ILP(il,ie,nr,nx).EQ.3) THEN !defined by nodes
          WRITE(OP_STRING,'(5X,''nodal values:'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(5(3X,I3,'': '',E10.3))')
     '      (NPNODE(nonode,nr),CP(il,NPNODE(nonode,nr)),
     '       nonode=1,NPNODE(0,nr))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IF(KTYP14.EQ.1) THEN
            WRITE(OP_STRING,'(5X,''nodal increments:'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(5(3X,I3,'': '',E10.3))')
     '        (NPNODE(nonode,nr),CP(il+ILT(ie,nr,nx),
     '         NPNODE(nonode,nr)),
     '         nonode=1,NPNODE(0,nr))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ELSE IF(ILP(il,1,nr,nx).EQ.4) THEN !defined by Gauss points
        ENDIF
      ENDDO

      CALL EXITS('OPMAT41')
      RETURN
 9999 CALL ERRORS('OPMAT41',ERROR)
      CALL EXITS('OPMAT41')
      RETURN 1
      END


