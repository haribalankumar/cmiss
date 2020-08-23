      SUBROUTINE OPMAT4(LIST,NEELEM,NPNODE,nr,NW,nx,CE,CP,ERROR,*)

C#### Subroutine: OPMAT4
C###  Description:
C###    OPMAT4 outputs material parameters.

      IMPLICIT NONE
      INCLUDE 'b10.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'titl40.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER LIST(0:NLISTM),NEELEM(0:NE_R_M,0:NRM),
     '  NPNODE(0:NP_R_M,0:NRM),nr,NW(NEM,3),nx
      REAL*8 CE(NMM,NEM),CP(NMM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,ie,l,l1,l2,ne,nj,noelem
      REAL*8 GACN,RHOL

      CALL ENTERS('OPMAT4',*9999)
      DO ie=1,12
        IF(ETYP(ie)) THEN
          l1=1
          l2=0
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            IF(NW(ne,1).EQ.ie) THEN
              l2=l2+1
              LIST(l2)=ne
            ENDIF
          ENDDO
          WRITE(OP_STRING,'(2X,A,9I5/,(31X,9I5))')
     '      TITL42(ie),(LIST(l),l=l1,l2)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IF(ie.NE.10) THEN
            CALL OPMAT41(ie,LIST,l1,l2,NPNODE,nr,nx,
     '        CE,CP,ERROR,*9999)
          ENDIF
        ENDIF
      ENDDO

      IF(ETYP(8).OR.ETYP(10)) THEN
        WRITE(OP_STRING,
     '     '(''  Additional parameters (for fluid-tank system):''/'
     '    //'''             Acceleration due to gravity='',F8.2,/'
     '    //'''             Fluid mass density         ='',F8.2,/)')
     '    GACN,RHOL
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF

      IF(ETYP(3).AND.NJT.EQ.3) THEN
        DO nj=1,2
          WRITE(OP_STRING,
     '      '(4X,''The '',A,'' of beam elements are :'')')
     '      TITL49(nj)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            IF(NW(ne,1).EQ.3) THEN
              WRITE(OP_STRING,'(5X,''Element '',I5,'' : '',3E13.5)')
     '          ne,(CE(ILT(3,nr,nx)+3*(nj-1)+I,ne),I=1,3)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO
        ENDDO
      ENDIF

      CALL EXITS('OPMAT4')
      RETURN
 9999 CALL ERRORS('OPMAT4',ERROR)
      CALL EXITS('OPMAT4')
      RETURN 1
      END


