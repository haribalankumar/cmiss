      SUBROUTINE OPMATE(ILPIN,LIST,NBJ,NEELEM,NMBIN,NPNODE,
     '  nr,NW,nx,nxc,CE,CGE,CIN,CP,CQ,YG,CONSTIT,ERROR,*)

C#### Subroutine: OPMATE
C###  Description:
C###    OPMATE outputs materials for class nx, region nr.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
!     Parameter List
      INTEGER ILPIN(NMM),LIST(0:NLISTM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NMBIN(NMM),
     '  NPNODE(0:NP_R_M,0:NRM),nr,NW(NEM,3),nx,nxc
      REAL*8 CE(NMM,NEM),CGE(NMM,NGM,NEM),
     '  CIN(NMM,0:NGM,NNEPM),CP(NMM,NPM),CQ(NMM,NQM),YG(NIYGM,NGM,NEM)
      CHARACTER ERROR*(*)
      LOGICAL CONSTIT
!     Local Variables
      INTEGER IBEG,IEND,nm,nonode,np

      CALL ENTERS('OPMATE',*9999)

      WRITE(OP_STRING,'(/'' Class '',I1,'' (nx='',I1,'
     '  //''') Region '',I1,'':'')') nxc,nx,nr
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

! Time-varying material parameters
      IF(ITYP5(nr,nx).EQ.2) THEN !Time integration problem
        IF(KTYP3_mate(nx).EQ.1) THEN
          WRITE(OP_STRING,'('' Material params are constant in time'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE IF(KTYP3_mate(nx).EQ.2) THEN
          WRITE(OP_STRING,'('' Time-varying material params are defined'
     '      //' in USER_IPMATE subroutine'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE IF(KTYP3_mate(nx).EQ.3) THEN
          CALL STRING_TRIM(FILE03,IBEG,IEND)
          WRITE(OP_STRING,'('' Time-varying material params are defined'
     '      //' in '//FILE03(IBEG:IEND)//'.IPMATE_time'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF !KTYP3_mate
      ENDIF !ityp5(nr)

      IF(ITYP1(nr,nx).EQ.3) THEN      !FE30 problems
        CALL OPMAT3(NBJ,NEELEM,NPNODE,nr,nx,CE,CP,CQ,YG,ERROR,*9999)
C LKC 25-MAR-1998 split off into own routine - opcell
C LKC 16-MAR-98 output cell parameters
C        IF(.NOT.CELL) THEN
C          CALL OPMAT3(NBJ,NEELEM,NPNODE,nr,nx,CE,CP,CQ,YG,ERROR,*9999)
C        ELSE
C          CALL OPMAT3_CELL(NBJ,NEELEM,NPNODE,nr,nx,CE,CP,CQ,YG,
C     '      ERROR,*9999)
C        ENDIF
      ELSE IF(ITYP1(nr,nx).EQ.4) THEN !FE40 problems
        CALL OPMAT4(LIST,NEELEM,NPNODE,nr,NW,nx,CE,CP,ERROR,*9999)
      ELSE IF(ITYP1(nr,nx).EQ.5) THEN !FE50 problems
        CALL OPMAT5(ILPIN,NBJ,NEELEM,NMBIN,NPNODE,nr,nx,CE,CGE,CIN,CP,
     '    CONSTIT,ERROR,*9999)
      ELSE IF(ITYP1(nr,nx).EQ.6) THEN !FE60 problems-OPMAT6 not done
        CALL OPMAT3(NBJ,NEELEM,NPNODE,nr,nx,CE,CP,CQ,YG,ERROR,*9999)
      ELSE IF(ITYP1(nr,nx).EQ.9) THEN !FE90 problems
C news AJP 12/4/95
        IF(ITYP5(nr,nx).EQ.1.AND.ITYP2(nr,nx).EQ.1) THEN
          !Static linear elasticity
          CALL OPMAT4(LIST,NEELEM,NPNODE,nr,NW,nx,CE,CP,ERROR,*9999)
        ELSE !Rest of BEM
          CALL OPMAT3(NBJ,NEELEM,NPNODE,nr,nx,CE,CP,CQ,YG,ERROR,*9999)
        ENDIF
C newe
      ELSE !used for material parameter fits
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          WRITE(OP_STRING,'(1P,'' CP(nm,'',I5,''): '',5E14.6)')
     '      np,(CP(nm,np),nm=1,6)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO !nonode
      ENDIF

      CALL EXITS('OPMATE')
      RETURN
 9999 CALL ERRORS('OPMATE',ERROR)
      CALL EXITS('OPMATE')
      RETURN 1
      END


