      SUBROUTINE OPMAT3(NBJ,NEELEM,NPNODE,nr,nx,CE,CP,CQ,YG,ERROR,*)

C#### Subroutine: OPMAT3
C###  Description:
C###    OPMAT3 outputs material parameters for FE30 problems.

      IMPLICIT NONE
      INCLUDE 'b13.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'titl30.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NPNODE(0:NP_R_M,0:NRM),nr
      REAL*8 CE(NMM,NEM),CP(NMM,NPM),CQ(NMM,NQM),YG(NIYGM,NGM,NEM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IB2,IBEG1,IBEG2,IBEG3,IE2,IEND1,IEND2,IEND3,il,ILT_TOT,
     '  J,JLP,nb,ne,ng,nh,NH_GROUP_TOT,nm,noelem,nonode,np,nq,nx
      CHARACTER CHAR1*34,CHAR2*8,TITL24(2)*17,TITLE5(5)*36
      DATA TITLE5 /'constant spatially                 ',
     '             'piecewise constant (element params)',
     '             'piecewise linear (nodal parameters)',
     '             'defined by grid points (CQ)        ',
     '             'defined by Gauss point array (YG)  '/
      DATA TITL24 /'constant wrt time','varying wrt time '/

      CALL ENTERS('OPMAT3',*9999)

C SGM 26Oct2000 grid-based Finite element also
C MLT 29Nov02 Grid Finite Volume also
      IF((ITYP2(nr,nx).EQ.9).AND.
     '  (ITYP4(nr,nx).EQ.4.OR.ITYP4(nr,nx).EQ.6.OR.
     '   ITYP4(nr,nx).EQ.7).AND.
     '  (ITYP5(nr,nx).EQ.2)) THEN
        ILT(0,nr,nx)=3
        nm=2
      ELSE
        ILT(0,nr,nx)=1
        nm=0
      ENDIF

C DMAL 22-OCT-2002 Adding loop over dependent variables
      IF(ITYP5(nr,nx).EQ.2) THEN !time-dependent
        IF(ITYP2(nr,nx).EQ.3) THEN !advection-diffusion
          NH_GROUP_TOT=KTYP3A(nx)
          ILT_TOT=ILT(1,nr,nx)/NH_GROUP_TOT
        ELSE
          NH_GROUP_TOT=1
          ILT_TOT=ILT(1,nr,nx)
        ENDIF
      ELSE
        NH_GROUP_TOT=1
        ILT_TOT=ILT(1,nr,nx)
      ENDIF

      DO nh=1,NH_GROUP_TOT
        DO il=ILT(0,nr,nx),ILT_TOT
          nm=nm+1
          JLP=IABS(ILP(il,1,nr,nx))
          IF(ILP(nm,1,nr,nx).GT.0) THEN
            J=1
          ELSE
            J=2
          ENDIF

          IF(ITYP5(nr,nx).EQ.1) THEN        !static analysis
            CHAR1=TITL31(il,ITYP2(nr,nx),ITYP3(nr,nx))
          ELSE IF(ITYP5(nr,nx).EQ.2) THEN   !time integration
            CHAR1=TITL32(il,ITYP2(nr,nx),ITYP3(nr,nx))
          ELSE IF(ITYP5(nr,nx).EQ.3) THEN   !modal analysis
            CHAR1=TITL33(il,ITYP2(nr,nx),ITYP3(nr,nx))
          ELSE IF(ITYP5(nr,nx).EQ.4) THEN   !Quasi-static analysis
            CHAR1=TITL31(il,ITYP2(nr,nx),ITYP3(nr,nx))
          ELSE IF(ITYP5(nr,nx).EQ.5) THEN   !Wavefront path analysis
            CHAR1=TITL35(il,ITYP2(nr,nx),ITYP3(nr,nx))
          ELSE IF(ITYP5(nr,nx).EQ.6) THEN   !buckling analysis
            CHAR1=TITL36(il,ITYP2(nr,nx),ITYP3(nr,nx))
          ENDIF
          WRITE(CHAR2,'(I2)') nh
          CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
          CALL STRING_TRIM(CHAR2,IB2,IE2)
          CALL STRING_TRIM(TITLE5(JLP),IBEG2,IEND2)
          CALL STRING_TRIM(TITL24(J),IBEG3,IEND3)
          IF(NH_GROUP_TOT.GT.1)THEN
            WRITE(OP_STRING,'(/I3,'') The '',A,'' for dependent'
     '        //' variable '',A/,'' is '',A,'
     '        //'4X,'' and '',A)') il,CHAR1(IBEG1:IEND1),CHAR2(IB2:IE2),
     '        TITLE5(JLP)(IBEG2:IEND2),TITL24(J)(IBEG3:IEND3)
          ELSE
            WRITE(OP_STRING,'(/I3,'') The '',A,'' is '',A/,'
     '        //'4X,'' and '',A)') il,CHAR1(IBEG1:IEND1),
     '        TITLE5(JLP)(IBEG2:IEND2),TITL24(J)(IBEG3:IEND3)
          ENDIF
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

          IF(ILP(nm,1,nr,nx).EQ.1) THEN      !constant spatially
            WRITE(OP_STRING,'('' Value='',E11.4)') CE(nm,NEELEM(1,nr))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

          ELSE IF(ILP(nm,1,nr,nx).EQ.2) THEN !defined by elements
            WRITE(OP_STRING,'('' Element values:'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              WRITE(OP_STRING,'('' Element '',I5,'': '',E12.5)')
     '          ne,CE(nm,ne)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDDO !noelem

          ELSE IF(ILP(nm,1,nr,nx).EQ.3) THEN !defined by nodes
            WRITE(OP_STRING,
     '        '('' Basis function type number is '',I2)') NMB(nm,1,nx)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            IF(NKT(0,NMB(nm,1,nx)).EQ.1) THEN      !linear interpolation
              WRITE(OP_STRING,'('' Nodal values:'')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'(5(1X,I5,'': '',E10.3))')
     '          (NPNODE(nonode,nr),CP(nm,NPNODE(nonode,nr)),
     '           nonode=1,NPNODE(0,nr))
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ELSE IF(NKT(0,NMB(nm,1,nx)).EQ.2) THEN !cubic Hermite interp.
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                WRITE(OP_STRING,
     '            '('' Node '',I5,'': Value= '',E12.3,'' Deriv=''
     '            ,E12.3)')
     '            np,CP(nm,2*np-1),CP(nm,2*np)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDDO
            ENDIF

          ELSE IF(ILP(nm,1,nr,nx).EQ.4) THEN !defined by grid points
            WRITE(OP_STRING,'('' Grid point values:'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            DO nq=1,NQT
              WRITE(OP_STRING,'(1X,I6,'': '',E10.3)') nq,CQ(nm,nq)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDDO

          ELSE IF(ILP(nm,1,nr,nx).EQ.5) THEN !defined by YG
            WRITE(OP_STRING,'('' Element Gauss point values:'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              nb=NBJ(1,ne) !basis for first geometric variable
              WRITE(OP_STRING,'('' Element '',I5,'':'')') ne
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'((10E12.4))') (YG(nm,ng,ne),ng=1,NGT(nb))
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDDO !noelem
          ENDIF !type of interpolation
        ENDDO !il
      ENDDO !nh

      CALL EXITS('OPMAT3')
      RETURN
 9999 CALL ERRORS('OPMAT3',ERROR)
      CALL EXITS('OPMAT3')
      RETURN 1
      END



