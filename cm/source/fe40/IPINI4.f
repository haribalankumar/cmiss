      SUBROUTINE IPINI4(IDO,INP,NBH,NBHF,NBJ,NEELEM,NEL,NENP,NFF,
     '  NGAP,NHE,NHP,NKEF,NKH,NKHE,NLL,NNF,NNL,NPF,NPL,NPLIST,NPNE,
     '  NPNODE,nr,NVHE,NVHP,NW,nx,NYNE,NYNP,DF,DL,PG,SE,
     '  WG,XP,YP,ZA,ZP,FIX,NOFIX,ERROR,*)

C#### Subroutine: IPINI4
C###  Description:
C###    Inputs initial conditions and boundary conditions for linear
C###    elasticity problems.

C**** Note: Load params defined by elements or nodes are assumed to be
C****       aligned with the element or global coord system,respect.ly,
C****       (in the latter case a rotation is made for beams & plates).
C**** Note: Direction cosines of normals to principal axes of beam
C****       elements in 3D space are carried in CE(ILT(2,nr,nx)+1...,ne)
C****       but  are not interpolated at the element Gauss points.

      IMPLICIT NONE
      CHARACTER*(*) ROUTINENAME
      PARAMETER(ROUTINENAME='IPINI4')
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b10.cmn'
      INCLUDE 'dx00.cmn'
      INCLUDE 'fluid00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBH(NHM,NCM,NEM),NBHF(NHM,NCM,NFM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NEL(0:NELM,NLM),NENP(NPM,0:NEPM,0:NRM),NFF(6,NEM),
     '  NGAP(NIM,NBM),NHE(NEM),
     '  NHP(NPM),NKHE(NKM,NNM,NHM,NEM),
     '  NKEF(0:4,16,6,NBFM),
     '  NKH(NHM,NPM,NCM),NLL(12,NEM),NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),
     '  NPF(9,NFM),NPL(5,0:3,NLM),NPLIST(0:NPM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),nr,
     '  NVHE(NNM,NBFM,NHM,NEM),
     '  NVHP(NHM,NPM,NCM),NW(NEM,3),nx,NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 DF(NFM),DL(3,NLM),PG(NSM,NUM,NGM,NBM),SE(NSM,NBFM,NEM),
     '  WG(NGM,NBM),XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM),
     '  ZA(NAM,NHM,NCM,NEM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM),NOFIX
!     Local Variables
      INTEGER INFO,mh,nb,
     '  ne,nh,nhx,noelem,NOQUES
      LOGICAL FILEIP

      CALL ENTERS(ROUTINENAME,*9999)

      FILEIP=.FALSE.
      NOQUES=0

C *** CE(14,ne) used for total pressure load on membranes
C     DO noelem=1,NEELEM(0,nr)
C       ne=NEELEM(noelem,nr)
C       CE(14,ne)=0.d0
C     ENDDO

      IF(ETYP(8).OR.ETYP(10)) THEN !fluid-coupled
        FORMAT='($,'' Specify the fluid mass density'//
     '    ' [1000.0]: '',E12.5)'
        RDEFLT(1)=1000.0d0
        IF(IOTYPE.EQ.3) RDATA(1)=RHOL
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,0,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) RHOL=RDATA(1)

        FORMAT='($,'' Specify the g-acceleration [9.81]: '',E12.5)'
        RDEFLT(1)=9.81d0
        IF(IOTYPE.EQ.3) RDATA(1)=GACN
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,0,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) GACN=RDATA(1)
      ENDIF

      IF(ETYP(7).OR.ETYP(8)) THEN !shell
        !To determine the 3rd basis functions (GEOM) derivs for shells.
        !Only the hermite basis is considered, since hermite is used for
        !shell geometry.
        nb=NBJ(1,NEELEM(1,nr))
        CALL ASSERT(NKT(0,nb).EQ.4,'>>Geometry must be cubic Hermite',
     '    ERROR,*9999)
        CALL GAUS20(IDO(1,1,0,nb),INP(1,1,nb),nb,NGAP(1,nb),
     '    D3PG(1,1,1,nb),ERROR,*9999)
      ENDIF

C     Check basis of beam element variables
      IF(ETYP(3)) THEN !some beam elements
        !Note: geom basis is assumed linear & axial displ. is cubic
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          IF(NW(ne,1).EQ.3) THEN !beam element
            CALL ASSERT(NKT(0,NBJ(1,ne)).EQ.1,
     '        '>>Geometry basis sb linear',ERROR,*9999)
            nb=NBH(NH_LOC(1,nx),1,ne)
            IF(nb.GT.0) CALL ASSERT(NKT(0,nb).EQ.2,
     '        '>>Axial displ. basis sb cubic',ERROR,*9999)
            nb=NBH(NH_LOC(2,nx),1,ne)
            IF(nb.GT.0) CALL ASSERT(NKT(0,nb).EQ.2,
     '        '>>Transverse displ. basis sb cubic',ERROR,*9999)
            nb=NBH(Nh_LOC(3,nx),1,ne)
            IF(nb.GT.0) CALL ASSERT(NKT(0,nb).EQ.2,
     '        '>>Transverse displ. basis sb cubic',ERROR,*9999)
          ENDIF
        ENDDO !noelem (ne)
      ENDIF

      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
        DO nhx=1,NHE(ne)
          nh=NH_LOC(nhx,nx) !GMH 20-Apr-95
          IF(NBH(nh,1,ne).EQ.0) THEN
            DO mh=1,NHE(ne)
              IF(NBH(mh,1,ne).NE.0) THEN
                NBH(nh,1,ne)=NBH(mh,1,ne)
                GO TO 39
              ENDIF
            ENDDO !mh
 39         CONTINUE
          ENDIF
        ENDDO !nhx (nh)
      ENDDO ! noelem (ne)

c     DO nonode=1,NPNODE(0,nr)
c       np=NPNODE(nonode,nr)
c       NHP(np)=0
c       DO nh=1,NHM
c         NKH(nh,np,1)=0
c       ENDDO
c       DO noelem=1,NEELEM(0,nr)
c         ne=NEELEM(noelem,nr)
c         nb=NBH(NH_LOC(1,nx),1,ne)
c         ie=NW(ne,1)
c         DO nn=1,NNT(nb)
c           IF(NPNE(nn,nb,ne).EQ.NP) THEN
c             IF(NHE(ne).GT.NHP(np)) NHP(np)=NHE(ne)
c             DO nvar=1,NVE(ie)
c               nh=NHV(nvar,ie)
c               nb=NBH(nh,1,ne)
c               IF(NKT(0,nb).GT.NKH(nh,np,1)) NKH(nh,np,1)=NKT(0,nb)
c             ENDDO
c             GO TO 44
c           ENDIF
c         ENDDO
c44       CONTINUE
c       ENDDO
c       DO nh=1,NHP(np)
c         IF(NKH(nh,np,1).EQ.0) THEN
c           DO mh=1,NHP(np)
c             IF(NKH(mh,np,1).NE.0) THEN
c               NKH(nh,np,1)=NKH(mh,np,1)
c               GO TO 452
c             ENDIF
c           ENDDO
c452        CONTINUE
c         ENDIF
c       ENDDO
c     ENDDO

      CALL IPINIT_ELAS(%VAL(0),NBH,NBHF,NEELEM,NEL,NENP,NFF,
     '  NHE,NHP,NKEF,NKH,NKHE,NLL,NNF,NNL,NPF,NPL,NPLIST,NPNE,
     '  NPNODE,nr,NVHE,NVHP,NW,nx,NYNE,NYNP,DF,DL,PG,SE,WG,XP,
     '  YP,ZA,ZP,FIX,NOFIX,ERROR,*9999)

      CALL EXITS(ROUTINENAME)
      RETURN
 9999 CALL ERRORS(ROUTINENAME,ERROR)
      CALL EXITS(ROUTINENAME)
      RETURN 1
      END


