      SUBROUTINE OPCELL_PROMPT(NEELEM,NPNODE,nr,nx,CE,CP,CQ,ERROR,*)

C#### Subroutine: OPCELL_PROMPT
C###  Description:
C###    OPCELL outputs a summary of cellular parameter numbers for
C###    simple ionic current models.

      IMPLICIT NONE

      INCLUDE 'b13.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'titl30.cmn'

!     Parameter List
      INTEGER NEELEM(0:NE_R_M,0:NRM),NPNODE(0:NP_R_M,0:NRM),nr,nx
      REAL*8 CE(NMM,NEM),CP(NMM,NPM),CQ(NMM,NQM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IBEG1,IBEG2,IBEG3,IEND1,IEND2,IEND3,il,ILT_LOCAL(0:1),J,
     '  JLP,ne,noelem,nonode,nq
      CHARACTER CHAR1*34,TITL24(2)*17,TITLE5(5)*36

      DATA TITLE5 /'constant spatially                 ',
     '             'piecewise constant (element params)',
     '             'piecewise linear (nodal parameters)',
     '             'defined by grid points (CQ)        ',
     '             'defined by Gauss point array (YG)  '/
      DATA TITL24 /'constant wrt time','varying wrt time '/

      CALL ENTERS('OPCELL_PROMPT',*9999)

      ILT_LOCAL(0)=9
      IF(ITYP3(nr,nx).EQ.1) THEN
        ILT_LOCAL(1)=12
      ELSE IF(ITYP3(nr,nx).EQ.2) THEN
        ILT_LOCAL(1)=18
      ELSE IF(ITYP3(nr,nx).EQ.3) THEN
        ILT_LOCAL(1)=14
      ELSE IF(ITYP3(nr,nx).EQ.4) THEN
        ILT_LOCAL(1)=13
      ENDIF

      DO il=1,2
        JLP=IABS(ILP(il,1,nr,nx))
        J=1
        CHAR1=TITL32(il,ITYP2(nr,nx),ITYP3(nr,nx))

        CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
        CALL STRING_TRIM(TITLE5(JLP),IBEG2,IEND2)
        CALL STRING_TRIM(TITL24(J),IBEG3,IEND3)
        WRITE(OP_STRING,'(/I3,'') The '',A,'' is '',A/,'
     '    //'4X,'' and '',A)') il,CHAR1(IBEG1:IEND1),
     '    TITLE5(JLP)(IBEG2:IEND2),TITL24(J)(IBEG3:IEND3)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        IF(ILP(il,1,nr,nx).EQ.1) THEN      !constant spatially
          WRITE(OP_STRING,'('' Value='',E11.4)') CQ(il,1)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        ELSE IF(ILP(il,1,nr,nx).EQ.2) THEN !defined by elements
          WRITE(OP_STRING,'('' Element values:'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            WRITE(OP_STRING,'('' Element '',I5,'': '',E12.5)')
     '        ne,CE(il,ne)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO !noelem

        ELSE IF(ILP(il,1,nr,nx).EQ.3) THEN !defined by nodes
          WRITE(OP_STRING,
     '      '('' Basis function type number is '',I2)') NMB(il,1,nx)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Nodal values:'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(5(1X,I5,'': '',E10.3))')
     '      (NPNODE(nonode,nr),CP(il,NPNODE(nonode,nr)),
     '      nonode=1,NPNODE(0,nr))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        ELSE IF(ILP(il,1,nr,nx).EQ.4) THEN !defined by grid points
          WRITE(OP_STRING,'('' Grid point values:'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO nq=1,NQT
            WRITE(OP_STRING,'(1X,I6,'': '',E10.3)') nq,CQ(il,nq)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO

        ELSE IF(ILP(il,1,nr,nx).EQ.5) THEN !defined by YG

        ENDIF !type of interpolation
      ENDDO

      DO il=ILT_LOCAL(0),ILT_LOCAL(1)
        JLP=IABS(ILP(il,1,nr,nx))
        J=1
        CHAR1=TITL32(il,ITYP2(nr,nx),ITYP3(nr,nx))

        CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
        CALL STRING_TRIM(TITLE5(JLP),IBEG2,IEND2)
        CALL STRING_TRIM(TITL24(J),IBEG3,IEND3)
        WRITE(OP_STRING,'(/I3,'') The '',A,'' is '',A/,'
     '    //'4X,'' and '',A)') il,CHAR1(IBEG1:IEND1),
     '    TITLE5(JLP)(IBEG2:IEND2),TITL24(J)(IBEG3:IEND3)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        IF(ILP(il,1,nr,nx).EQ.1) THEN      !constant spatially
          WRITE(OP_STRING,'('' Value='',E11.4)') CQ(il,1)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        ELSE IF(ILP(il,1,nr,nx).EQ.2) THEN !defined by elements
          WRITE(OP_STRING,'('' Element values:'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            WRITE(OP_STRING,'('' Element '',I5,'': '',E12.5)')
     '        ne,CE(il,ne)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO !noelem

        ELSE IF(ILP(il,1,nr,nx).EQ.3) THEN !defined by nodes
          WRITE(OP_STRING,
     '      '('' Basis function type number is '',I2)') NMB(il,1,nx)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Nodal values:'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(5(1X,I5,'': '',E10.3))')
     '      (NPNODE(nonode,nr),CP(il,NPNODE(nonode,nr)),
     '      nonode=1,NPNODE(0,nr))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        ELSE IF(ILP(il,1,nr,nx).EQ.4) THEN !defined by grid points
          WRITE(OP_STRING,'('' Grid point values:'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO nq=1,NQT
            WRITE(OP_STRING,'(1X,I6,'': '',E10.3)') nq,CQ(il,nq)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO

        ELSE IF(ILP(il,1,nr,nx).EQ.5) THEN !defined by YG

        ENDIF !type of interpolation
      ENDDO

      CALL EXITS('OPCELL_PROMPT')
      RETURN
 9999 CALL ERRORS('OPCELL_PROMPT',ERROR)
      CALL EXITS('OPCELL_PROMPT')
      RETURN 1
      END


