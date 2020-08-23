      SUBROUTINE UPGROW(NBJ,NEELEM,NRLIST,NXLIST,XP,YG,STRING,ERROR,*)

C#### Subroutine: UPGROW
C###  Description:
C###    UPGROW updates growth by:
C###    (ktyp60=1,2) updating material density from strain energy
C###      density;
C###    (ktyp60=5)   updating position of growing tip.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'coup00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grow00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp60.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NRLIST(0:NRM),
     '  NXLIST(0:NXM)
      REAL*8 XP(NKM,NVM,NJM,NPM),YG(NIYGM,NGM,NEM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!    Local Variables
      INTEGER IBEG,IEND,nb,ne,ng,noelem,np,no_tip,
     '  no_nrlist,nr,nx,nxc
      LOGICAL ALL_REGIONS

      CALL ENTERS('UPGROW',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM update growth
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description: updates growth by:
C###    (ktyp60=1,2) updating material density from strain energy
C###      density;
C###    (ktyp60=5)   updating position of growing tip.
C###

        OP_STRING(1)=STRING(1:IEND)
        OP_STRING(2)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(3)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe23','doc','UPGROW',ERROR,*9999)
      ELSE
        CALL ASSERT(KTYP60.GT.0,'>>Growth law not defined',ERROR,*9999)

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)

        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '     ERROR,*9999)

        IF(KTYP60.LE.2) THEN !update material density from S.E. density
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            IF(ITYP2(nr,nx).EQ.1) THEN      !linear elasticity
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                nb=NBJ(1,ne)
                DO ng=1,NGT(nb)
                  YG(5,ng,ne)=GROW1+GROW2*YG(4,ng,ne) !YG(4)=SE,
                ENDDO                                 !YG(5)=density
              ENDDO

            ELSE IF(ITYP2(nr,nx).EQ.5) THEN !Poisson eqtn

            ENDIF
          ENDDO !no_nrlist
        ELSE IF(KTYP60.EQ.5) THEN !update position of growing tip
          ne=1
          DO no_tip=1,NP_GROWING_TIP(0,ne)
            np=NP_GROWING_TIP(no_tip,ne)
C GMH 8/1/97 Update cmgui link
            CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
            XP(1,1,2,np)=XP(1,1,2,np)+GROW1
          ENDDO
        ENDIF
      ENDIF

      CALL EXITS('UPGROW')
      RETURN
 9999 CALL ERRORS('UPGROW',ERROR)
      CALL EXITS('UPGROW')
      RETURN 1
      END


