      SUBROUTINE UPCOUP(NBH,NEELEM,NHE,NHP,NKH,NKJ,
     '  NPNODE,NP_INTERFACE,NVHP,NVJP,NYNE,NYNP,
     '  PF,XP,YP,ZA,ZP,STRING,ERROR,*)

C#### Subroutine: UPCOUP
C###  Description:
C###    UPCOUP pdates coupling between regions as setup by 'define
C###    coupling'.
C###    For coupled aerofoil flow & stress (KTYP90=5):
C###      NP_INTERFACE(0,i),i=1,3 are #node triplets on aerofoil;
C###      NP_INTERFACE(n,i),i=1,3 are node #s on aerofoil;
C###      i=1 is upper surface aerofoil flow field;
C###      i=2 is lower surface aerofoil flow field;
C###      i=3 is aerofoil stress field.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'aero00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBH(NHM,NCM,NEM),NEELEM(0:NE_R_M,0:NRM),NHE(NEM,NXM),
     '  NHP(NPM,0:NRM,NXM),NKH(NHM,NPM,NCM,0:NRM),NKJ(NJM,NPM),
     '  NPNODE(0:NP_R_M,0:NRM),NP_INTERFACE(0:NPM,0:3),
     '  NVHP(NHM,NPM,NCM,0:NRM),NVJP(NJM,NPM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 PF(2,NEM),XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM),
     '  ZA(NAM,NHM,NCM,NEM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,N,nc,ne,nj,nk,noelem,nr,NP1,NP2,NP3,nv,nx

      CALL ENTERS('UPCOUP',*9999)

      nc=1 !temporary cpb 22/11/94
      nx=1 !temporary cpb 22/11/94

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM update coupling
C###  Description: Udates coupling between regions.
C###    For coupled aerofoil flow & stress (KTYP90=5):
C###      NP_INTERFACE(0,i),i=1,3 are #node triplets on aerofoil;
C###      NP_INTERFACE(n,i),i=1,3 are node #s on aerofoil;
C###      i=1 is upper surface aerofoil flow field;
C###      i=2 is lower surface aerofoil flow field;
C###      i=3 is aerofoil stress field.
C###

        OP_STRING(1)=STRING(1:IEND)
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe23','doc','UPCOUP',ERROR,*9999)
      ELSE

        IF(KTYP90.EQ.5) THEN !Coupled aerofoil flow & stress
          !Put aero pressure difference into sail element pressure
          nr=2 !for sail
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            PF(1,ne)=PRESS_DIFF_AERO(2,noelem) !middle gauss point
          ENDDO

          !Update flow field nodes from deformed sail nodes
!new MPN 6-Jan-95: current soln now stored in YP(ny,1) for nonlin probs
          CALL YPZP(1,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),NKH(1,1,1,nr),
     '      NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,
     '      YP(1,1,nx),ZA,ZP,ERROR,*9999)
!old
c          CALL YPZP(4,NBH,NEELEM,NHE(1,nx),NHP(1,nr,nx),NKH(1,1,1,nr),
c     '      NPNODE,nr,NVHP(1,1,1,nr),nx,NYNE,NYNP,
c     '      YP(1,1,nx),ZA,ZP,ERROR,*9999)
          DO N=1,NP_INTERFACE(0,3)
            NP1=NP_INTERFACE(N,1) !upper surface flow field nodes
            NP2=NP_INTERFACE(N,2) !lower surface flow field nodes
            NP3=NP_INTERFACE(N,3) !sail stress field nodes
C GMH 8/1/97 Update cmgui link
            CALL NODE_CHANGE(NP1,.FALSE.,ERROR,*9999)
C GMH 8/1/97 Update cmgui link
            CALL NODE_CHANGE(NP2,.FALSE.,ERROR,*9999)
            DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
              DO nk=1,NKJ(nj,NP3)
                DO nv=1,NVJP(nj,NP1)
                  XP(nk,nv,nj,NP1)=ZP(nk,nv,nj,NP3,nc)
                ENDDO !nv
                DO nv=1,NVJP(nj,NP2)
                  XP(nk,nv,nj,NP2)=ZP(nk,nv,nj,NP3,nc)
                ENDDO !nv
              ENDDO !nk
            ENDDO !nj
          ENDDO !N
        ENDIF

      ENDIF

      CALL EXITS('UPCOUP')
      RETURN
 9999 CALL ERRORS('UPCOUP',ERROR)
      CALL EXITS('UPCOUP')
      RETURN 1
      END


