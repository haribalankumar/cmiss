      SUBROUTINE UPMOTI(NXLIST,STRING,ERROR,*)

C#### Subroutine: UPMOTI
C###  Description:

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn'

!     Parameter List
      INTEGER NXLIST(0:NXM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,nx,nxc,nx_from,N3CO
      REAL*8 RFROMC
      LOGICAL CBBREV

      CALL ENTERS('UPMOTI',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe23','doc','UPMOTI',ERROR,*9999)
      ELSE

        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        IF(CBBREV(CO,'INLET_FLOW',2,noco+1,NTCO,N3CO)) 
     &    INLET_FLOW(nx)=RFROMC(CO(N3CO+1)) !*1.d6 !(L to mm^3)
        IF(CBBREV(CO,'DV_FROM',2,noco+1,NTCO,N3CO)) 
     &    DV_TOTAL(nx)=DV_TOTAL(IFROMC(CO(N3CO+1)))

      ENDIF

      CALL EXITS('UPMOTI')
      RETURN
 9999 CALL ERRORS('UPMOTI',ERROR)
      CALL EXITS('UPMOTI')
      RETURN 1
      END


