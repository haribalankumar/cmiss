      SUBROUTINE UPMASS(NEELEM,NRLIST,NXLIST,BBM,STRING,ERROR,*)

C#### Subroutine: UPMASS
C###  Description:
C###   AJS - include update of class-dependent CURRENT_VOLUME

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
!     Parameter List
      INTEGER NEELEM(0:NE_R_M,0:NRM),NRLIST(0:NRM),NXLIST(0:NXM)
      REAL*8 BBM(2,NEM,NXM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,ne,noelem,nr,nx,nxc,nx_from,N3CO
      INTEGER IFROMC
      LOGICAL ALL_REGIONS,CBBREV

      CALL ENTERS('UPMASS',*9999)


      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe23','doc','UPMASS',ERROR,*9999)
      ELSE

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        nr=NRLIST(1)
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        IF(CBBREV(CO,'FROM',2,noco+1,NTCO,N3CO)) 
     &    nx_from=IFROMC(CO(N3CO+1))
        MASS_CURRENT(nx)=MASS_CURRENT(nx_from)
        MASS_PREVIOUS(nx)=MASS_PREVIOUS(nx_from)
        CURRENT_VOLUME(nx)=CURRENT_VOLUME(nx_from)
        IDEAL_MASS(nx)=IDEAL_MASS(nx_from)

! AJS 29/3/2011 BBM is not class-dependent in FEM_DYNAM (or anywhere else in the code)
! Therefore remove updates for BBM below (do not need to copy volume and concentration 
! from inspiration to expiration class)
!         DO noelem=1,NEELEM(0,nr)
!           ne=NEELEM(noelem,nr)
!           BBM(1,ne,nx)=BBM(1,ne,nx_from)
!           BBM(2,ne,nx)=BBM(2,ne,nx_from)
!         ENDDO !noelem

      ENDIF

      CALL EXITS('UPMASS')
      RETURN
 9999 CALL ERRORS('UPMASS',ERROR)
      CALL EXITS('UPMASS')
      RETURN 1
      END


