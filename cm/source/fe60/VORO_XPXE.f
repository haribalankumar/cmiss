      SUBROUTINE VORO_XPXE(NBJ,ne,NKJE,NPNE,NVJE,SE,XE,XP,ERROR,*)

C#### Subroutine: VORO_XPXE
C###  Description:
C###    VORO_XPXE puts XP geometry information into XE.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),ne,NKJE(NKM,NNM,NJM,NEM),
     '  NPNE(NNM,NBFM,NEM),NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 SE(NSM,NBFM,NEM),XE(NSM,NJM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER mk,nb,nj,nk,nn,np,ns,nv

      CALL ENTERS('VORO_XPXE',*9999)

      IF(DOP)THEN
        WRITE(OP_STRING,'('' ne='',I5,'' nb='',I4)') ne,NBJ(1,ne)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      DO nj=1,3 !3d coordinates
        nb=NBJ(nj,ne) !basis fn
        ns=0
        DO nn=1,NNT(nb) !for each element node for basis nb
          np=NPNE(nn,nb,ne) !global node #
          nv=NVJE(nn,nb,nj,ne) !global version # for ne
          DO mk=1,NKT(nn,nb) !for # of nodal derivatives
            ns=ns+1 !element dof
            nk=NKJE(mk,nn,nj,ne) !global derivative #
            IF(DOP)THEN
              WRITE(OP_STRING,'('' nk,nv,nj,np,ns='',5(I4))')
     '          nk,nv,nj,np,ns
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            XE(ns,nj)=XP(nk,nv,nj,np)*SE(ns,nb,ne)
          ENDDO !mk
        ENDDO !nn
      ENDDO !nj

      CALL EXITS('VORO_XPXE')
      RETURN
 9999 CALL ERRORS('VORO_XPXE',ERROR)
      CALL EXITS('VORO_XPXE')
      RETURN 1
      END


